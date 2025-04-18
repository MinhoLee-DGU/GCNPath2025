import logging
import sys
from collections import OrderedDict

import pytoda
import torch
import torch.nn as nn
from pytoda.smiles.transforms import AugmentTensor

from ..utils.hyperparams import ACTIVATION_FN_FACTORY, LOSS_FN_FACTORY
from ..utils.interpret import monte_carlo_dropout, test_time_augmentation
from ..utils.layers import (
    ContextAttentionLayer, convolutional_layer, dense_layer
)
from ..utils.utils import get_device, get_log_molar

# setup logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger(__name__)


class PaccMannV2(nn.Module):
    """Based on the MCA model in Molecular Pharmaceutics:
        https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.9b00520.
        Updates:
            - Context instead of self attention on omic data
    """

    def __init__(self, params, *args, **kwargs):
        """Constructor.

        Args:
            params (dict): A dictionary containing the parameter to built the
                dense encoder.
                TODO params should become actual arguments (use **params).

        Items in params:
            smiles_padding_length (int): Padding length for SMILES.
            smiles_embedding_size (int): dimension of tokens' embedding.
            smiles_vocabulary_size (int): size of the tokens vocabulary.
            activation_fn (string, optional): Activation function used in all
                layers for specification in ACTIVATION_FN_FACTORY.
                Defaults to 'relu'.
            batch_norm (bool, optional): Whether batch normalization is
                applied. Defaults to True.
            dropout (float, optional): Dropout probability in all
                except parametric layer. Defaults to 0.5.
            filters (list[int], optional): Numbers of filters to learn per
                SMILES convolutional layer. Defaults to [64, 64, 64].
            kernel_sizes (list[list[int]], optional): Sizes of kernels per
                SMILES convolutional layer. Defaults to  [
                    [3, params['smiles_embedding_size']],
                    [5, params['smiles_embedding_size']],
                    [11, params['smiles_embedding_size']]
                ]
                NOTE: The kernel sizes should match the dimensionality of the
                smiles_embedding_size, so if the latter is 8, the images are
                t x 8, then treat the 8 embedding dimensions like channels
                in an RGB image.
            molecule_heads (list[int], optional): Amount of attentive molecule_heads
                per SMILES embedding. Should have len(filters)+1.
                Defaults to [4, 4, 4, 4].
            stacked_dense_hidden_sizes (list[int], optional): Sizes of the
                hidden dense layers. Defaults to [1024, 512].
            smiles_attention_size (int, optional): size of the attentive layer
                for the smiles sequence. Defaults to 64.
    """

        super(PaccMannV2, self).__init__(*args, **kwargs)

        # Model Parameter
        self.device = get_device()
        self.params = params
        self.loss_fn = LOSS_FN_FACTORY[params.get('loss_fn', 'mse')]
        self.min_max_scaling = True if params.get(
            'drug_sensitivity_processing_parameters', {}
        ) != {} else False
        if self.min_max_scaling:
            self.IC50_max = params[
                'drug_sensitivity_processing_parameters'
            ]['parameters']['max']  # yapf: disable
            self.IC50_min = params[
                'drug_sensitivity_processing_parameters'
            ]['parameters']['min']  # yapf: disable

        # Model inputs
        self.smiles_padding_length = params['smiles_padding_length']
        self.number_of_genes = params.get('number_of_genes', 2128)
        self.smiles_attention_size = params.get('smiles_attention_size', 64)
        self.gene_attention_size = params.get('gene_attention_size', 1)
        self.molecule_temperature = params.get('molecule_temperature', 1.)
        self.gene_temperature = params.get('gene_temperature', 1.)

        # Model architecture (hyperparameter)
        self.molecule_heads = params.get('molecule_heads', [4, 4, 4, 4])
        self.gene_heads = params.get('gene_heads', [2, 2, 2, 2])
        if len(self.gene_heads) != len(self.molecule_heads):
            raise ValueError('Length of gene and molecule_heads do not match.')

        self.filters = params.get('filters', [64, 64, 64])

        self.hidden_sizes = (
            [
                self.molecule_heads[0] * params['smiles_embedding_size'] + sum(
                    [
                        h * f
                        for h, f in zip(self.molecule_heads[1:], self.filters)
                    ]
                ) + sum(self.gene_heads) * self.number_of_genes
            ] + params.get('stacked_dense_hidden_sizes', [1024, 512])
        )

        self.dropout = params.get('dropout', 0.5)
        self.temperature = params.get('temperature', 1.)
        self.act_fn = ACTIVATION_FN_FACTORY[
            params.get('activation_fn', 'relu')]
        self.kernel_sizes = params.get(
            'kernel_sizes', [
                [3, params['smiles_embedding_size']],
                [5, params['smiles_embedding_size']],
                [11, params['smiles_embedding_size']]
            ]
        )
        if len(self.filters) != len(self.kernel_sizes):
            raise ValueError(
                'Length of filter and kernel size lists do not match.'
            )
        if len(self.filters) + 1 != len(self.molecule_heads):
            raise ValueError(
                'Length of filter and multihead lists do not match'
            )

        # Build the model
        self.smiles_embedding = nn.Embedding(
            self.params['smiles_vocabulary_size'],
            self.params['smiles_embedding_size'],
            scale_grad_by_freq=params.get('embed_scale_grad', False)
        )

        self.convolutional_layers = nn.Sequential(
            OrderedDict(
                [
                    (
                        f'convolutional_{index}',
                        convolutional_layer(
                            num_kernel,
                            kernel_size,
                            act_fn=self.act_fn,
                            batch_norm=params.get('batch_norm', False),
                            dropout=self.dropout
                        ).to(self.device)
                    ) for index, (num_kernel, kernel_size) in
                    enumerate(zip(self.filters, self.kernel_sizes))
                ]
            )
        )

        smiles_hidden_sizes = [params['smiles_embedding_size']] + self.filters

        self.molecule_attention_layers = nn.Sequential(OrderedDict([
            (
                f'molecule_attention_{layer}_head_{head}',
                ContextAttentionLayer(
                    reference_hidden_size=smiles_hidden_sizes[layer],
                    reference_sequence_length=self.smiles_padding_length,
                    context_hidden_size=1,
                    context_sequence_length=self.number_of_genes,
                    attention_size=self.smiles_attention_size,
                    individual_nonlinearity=params.get(
                        'context_nonlinearity', nn.Sequential()
                    ),
                    temperature=self.molecule_temperature
                )
            ) for layer in range(len(self.molecule_heads))
            for head in range(self.molecule_heads[layer])
            ]))  # yapf: disable

        # Gene attention stream
        self.gene_attention_layers = nn.Sequential(OrderedDict([
            (
                f'gene_attention_{layer}_head_{head}',
                ContextAttentionLayer(
                    reference_hidden_size=1,
                    reference_sequence_length=self.number_of_genes,
                    context_hidden_size=smiles_hidden_sizes[layer],
                    context_sequence_length=self.smiles_padding_length,
                    attention_size=self.gene_attention_size,
                    individual_nonlinearity=params.get(
                        'context_nonlinearity', nn.Sequential()
                    ),
                    temperature=self.gene_temperature
                )
            ) for layer in range(len(self.molecule_heads))
            for head in range(self.gene_heads[layer])
            ]))  # yapf: disable

        # Only applied if params['batch_norm'] = True
        self.batch_norm = nn.BatchNorm1d(self.hidden_sizes[0])
        self.dense_layers = nn.Sequential(
            OrderedDict(
                [
                    (
                        'dense_{}'.format(ind),
                        dense_layer(
                            self.hidden_sizes[ind],
                            self.hidden_sizes[ind + 1],
                            act_fn=self.act_fn,
                            dropout=self.dropout,
                            batch_norm=params.get('batch_norm', True)
                        ).to(self.device)
                    ) for ind in range(len(self.hidden_sizes) - 1)
                ]
            )
        )

        self.final_dense = (
            nn.Linear(self.hidden_sizes[-1], 1)
            if not params.get('final_activation', False) else nn.Sequential(
                OrderedDict(
                    [
                        ('projection', nn.Linear(self.hidden_sizes[-1], 1)),
                        ('sigmoidal', ACTIVATION_FN_FACTORY['sigmoid'])
                    ]
                )
            )
        )

    def forward(self, smiles, gep, confidence=False):
        """Forward pass through the PaccMannV2.

        Args:
            smiles (torch.Tensor): of type int and shape: [bs, smiles_padding_length]
            gep (torch.Tensor): of shape `[bs, number_of_genes]`.
            confidence (bool, optional) whether the confidence estimates are
                performed.

        Returns:
            (torch.Tensor, dict): predictions, prediction_dict
            predictions is IC50 drug sensitivity prediction of shape `[bs, 1]`.
            prediction_dict includes the prediction and attention weights.
        """

        gep = torch.unsqueeze(gep, dim=-1)
        # [bs, 2128] > [bs, 2128, 1]
        embedded_smiles = self.smiles_embedding(smiles.to(dtype=torch.int64))
        # [bs, 560] > [bs, 560, 128]
        
        # SMILES Convolutions. Unsqueeze has shape bs x 1 x T x H.
        encoded_smiles = [embedded_smiles] + [
            self.convolutional_layers[ind]
            (torch.unsqueeze(embedded_smiles, 1)).permute(0, 2, 1)
            for ind in range(len(self.convolutional_layers))
        ]
        ### Here is the location of error
        # [32, 560, 128], [32, 560, 64], [32, 560, 64], [32, 560, 64]

        # Molecule context attention
        encodings, smiles_alphas, gene_alphas = [], [], []
        for layer in range(len(self.molecule_heads)):
            for head in range(self.molecule_heads[layer]):

                ind = self.molecule_heads[0] * layer + head
                e, a = self.molecule_attention_layers[ind](
                    encoded_smiles[layer], gep
                )
                encodings.append(e)
                # [bs, 128] x 4
                # [bs, 64] x 12
                smiles_alphas.append(a)
                # [bs, 560] x 16

        # Gene context attention
        for layer in range(len(self.gene_heads)):
            for head in range(self.gene_heads[layer]):
                ind = self.gene_heads[0] * layer + head

                e, a = self.gene_attention_layers[ind](
                    gep, encoded_smiles[layer], average_seq=False
                )
                encodings.append(e)
                # [bs, 2128] x 8
                gene_alphas.append(a)
                # [bs, 2128] x 16

        encodings = torch.cat(encodings, dim=1)
        # [bs, 18304]
        
        # Apply batch normalization if specified
        inputs = self.batch_norm(encodings) if self.params.get(
            'batch_norm', False
        ) else encodings
        # [bs, 18304]
        # NOTE: stacking dense layers as a bottleneck
        for dl in self.dense_layers:
            inputs = dl(inputs)
        # [bs, 512]
        predictions = self.final_dense(inputs)
        # [bs, 1]
        prediction_dict = {}

        if not self.training:
            # The below is to ease postprocessing
            smiles_attention = torch.cat(
                [torch.unsqueeze(p, -1) for p in smiles_alphas], dim=-1
            )
            gene_attention = torch.cat(
                [torch.unsqueeze(p, -1) for p in gene_alphas], dim=-1
            )
            prediction_dict.update({
                'gene_attention': gene_attention,
                'smiles_attention': smiles_attention,
                'IC50': predictions,
                'log_micromolar_IC50':
                    get_log_molar(
                        predictions,
                        ic50_max=self.IC50_max,
                        ic50_min=self.IC50_min
                    ) if self.min_max_scaling else predictions
            })  # yapf: disable

            if confidence:
                augmenter = AugmentTensor(self.smiles_language)
                epi_conf, epi_pred = monte_carlo_dropout(
                    self,
                    regime='tensors',
                    tensors=(smiles, gep),
                    repetitions=5
                )
                ale_conf, ale_pred = test_time_augmentation(
                    self,
                    regime='tensors',
                    tensors=(smiles, gep),
                    repetitions=5,
                    augmenter=augmenter,
                    tensors_to_augment=0
                )

                prediction_dict.update({
                    'epistemic_confidence': epi_conf,
                    'epistemic_predictions': epi_pred,
                    'aleatoric_confidence': ale_conf,
                    'aleatoric_predictions': ale_pred
                })  # yapf: disable

        elif confidence:
            logger.info('Using confidence in training mode is not supported.')

        return predictions, prediction_dict

    def loss(self, yhat, y):
        return self.loss_fn(yhat, y)

    def _associate_language(self, smiles_language):
        """
        Bind a SMILES language object to the model. Is only used inside the
        confidence estimation.

        Arguments:
            smiles_language {[pytoda.smiles.smiles_language.SMILESLanguage]}
            -- [A SMILES language object]

        Raises:
            TypeError:
        """
        if not isinstance(
            smiles_language, pytoda.smiles.smiles_language.SMILESLanguage
        ):
            raise TypeError(
                'Please insert a smiles language (object of type '
                'pytoda.smiles.smiles_language.SMILESLanguage). Given was '
                f'{type(smiles_language)}'
            )
        self.smiles_language = smiles_language

    def load(self, path, *args, **kwargs):
        """Load model from path."""
        weights = torch.load(path, *args, **kwargs)
        self.load_state_dict(weights)

    def save(self, path, *args, **kwargs):
        """Save model to path."""
        torch.save(self.state_dict(), path, *args, **kwargs)

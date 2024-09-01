# import os
# os.environ['CUDA_LAUNCH_BLOCKING'] = "1"
# os.environ["CUDA_VISIBLE_DEVICES"] = "0"

#!/usr/bin/env python
"""Train PaccMann predictor."""
import argparse
import json
import logging
import os
import pickle
import sys
from copy import deepcopy
from time import time

import numpy as np
import torch
from paccmann_predictor.models import MODEL_FACTORY
from paccmann_predictor.utils.hyperparams import OPTIMIZER_FACTORY
from paccmann_predictor.utils.loss_functions import pearsonr
from paccmann_predictor.utils.utils import get_device
from pytoda.datasets import DrugSensitivityDataset
from pytoda.smiles.smiles_language import SMILESTokenizer

from utils_def import *
from utils_paccmann import *
# torch.multiprocessing.set_sharing_strategy('file_system')


# setup logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

# yapf: disable
parser = argparse.ArgumentParser()
# parser.add_argument(
#     'train_sensitivity_filepath', type=str,
#     help='Path to the drug sensitivity (IC50) data.'
# )
# parser.add_argument(
#     'test_sensitivity_filepath', type=str,
#     help='Path to the drug sensitivity (IC50) data.'
# )
# parser.add_argument(
#     'gep_filepath', type=str,
#     help='Path to the gene expression profile data.'
# )
# parser.add_argument(
#     'smi_filepath', type=str,
#     help='Path to the SMILES data.'
# )
# parser.add_argument(
#     'gene_filepath', type=str,
#     help='Path to a pickle object containing list of genes.'
# )
# parser.add_argument(
#     'smiles_language_filepath', type=str,
#     help='Path to a folder with SMILES language .json files.'
# )
# parser.add_argument(
#     'model_path', type=str,
#     help='Directory where the model will be stored.'
# )
# parser.add_argument(
#     'params_filepath', type=str,
#     help='Path to the parameter file.'
# )
# parser.add_argument(
#     'training_name', type=str,
#     help='Name for the training.'
# )
# yapf: enable

# Cell and Drug
parser.add_argument("-nth", type=int, default=0)
parser.add_argument("-cpu", type=int, default=4)
parser.add_argument("-choice", type=int, default=0)                  
parser.add_argument("-ic50", type=str, default="_data/IC50_GDSC2.csv")

parser.add_argument("-col_cell", type=str, default="Cell")
parser.add_argument("-col_drug", type=str, default="Drug")
parser.add_argument("-col_ic50", type=str, default="LN_IC50")

# def main(
#     train_sensitivity_filepath,
#     test_sensitivity_filepath,
#     gep_filepath,
#     smi_filepath,
#     gene_filepath,
#     smiles_language_filepath,
#     model_path,
#     params_filepath,
#     training_name,
# ):
  
def main(
    train_sensitivity_filepath,
    # valid_sensitivity_filepath,
    test_sensitivity_filepath,
    gep_filepath,
    smi_filepath,
    gene_filepath,
    smiles_language_filepath,
    params_filepath,
    args
):
    
    # logger = logging.getLogger(f"{training_name}")
    logger = logging.getLogger(args.dir_param)
    # Process parameter file:
    params = {}
    with open(params_filepath) as fp:
        params.update(json.load(fp))

    # Create model directory and dump files
    params["epochs"] = args.epochs
    params["device"] = args.device
    params["patience"] = args.patience
    # model_dir = os.path.join(model_path, training_name)
    # os.makedirs(os.path.join(model_dir, "weights"), exist_ok=True)
    # os.makedirs(os.path.join(model_dir, "results"), exist_ok=True)
    # with open(os.path.join(model_dir, "model_params.json"), "w") as fp:
    
    test_only = False
    if os.path.isfile(args.dir_hparam) and os.path.isfile(args.dir_param) :
        with open(args.dir_hparam, "r") as f : 
            params = json.load(f)
            params["device"] = args.device
        
        if params["early_stop_count"]>=params["patience"] :
            print("# The model is already trained to early stopping...")
            if os.path.isfile(args.dir_test) :
                print("# Test data is already predicted...")
                sys.exit()
            else :
                test_only = True
                print("# Test data is not predicted yet...")
        else :
            print("# Resume the model training...")
            
        print("# [epoch {}, patience : {}, count : {}]".format(params["epoch"], params["patience"], params["early_stop_count"]))
        checkpoint = torch.load(args.dir_param, map_location=params["device"])
        
    else :
        with open(args.dir_hparam, "w") as fp :
            json.dump(params, fp, indent=4)
    
    # Prepare the dataset
    logger.info("Start data preprocessing...")
    
    # Load SMILES language
    smiles_language = SMILESTokenizer.from_pretrained(smiles_language_filepath)
    smiles_language.set_encoding_transforms(
        add_start_and_stop=params.get("add_start_and_stop", True),
        padding=params.get("padding", True),
        padding_length=params.get("smiles_padding_length", None),
    )
    test_smiles_language = deepcopy(smiles_language)
    smiles_language.set_smiles_transforms(
        augment=params.get("augment_smiles", False),
        canonical=params.get("smiles_canonical", False),
        kekulize=params.get("smiles_kekulize", False),
        all_bonds_explicit=params.get("smiles_bonds_explicit", False),
        all_hs_explicit=params.get("smiles_all_hs_explicit", False),
        remove_bonddir=params.get("smiles_remove_bonddir", False),
        remove_chirality=params.get("smiles_remove_chirality", False),
        selfies=params.get("selfies", False),
        sanitize=params.get("selfies", False),
    )
    test_smiles_language.set_smiles_transforms(
        augment=False,
        canonical=params.get("test_smiles_canonical", True),
        kekulize=params.get("smiles_kekulize", False),
        all_bonds_explicit=params.get("smiles_bonds_explicit", False),
        all_hs_explicit=params.get("smiles_all_hs_explicit", False),
        remove_bonddir=params.get("smiles_remove_bonddir", False),
        remove_chirality=params.get("smiles_remove_chirality", False),
        selfies=params.get("selfies", False),
        sanitize=params.get("selfies", False),
    )
    
    # Load the gene list
    with open(gene_filepath, "rb") as f:
        gene_list = pickle.load(f)

    # Assemble datasets
    train_dataset = DrugSensitivityDataset(
        drug_sensitivity_filepath=train_sensitivity_filepath,
        column_names=[args.col_drug, args.col_cell, args.col_ic50],
        smi_filepath=smi_filepath,
        gene_expression_filepath=gep_filepath,
        smiles_language=smiles_language,
        gene_list=gene_list,
        drug_sensitivity_min_max=params.get("drug_sensitivity_min_max", True),
        drug_sensitivity_processing_parameters=params.get(
            "drug_sensitivity_processing_parameters", {}
        ),
        gene_expression_standardize=params.get("gene_expression_standardize", True),
        gene_expression_min_max=params.get("gene_expression_min_max", False),
        gene_expression_processing_parameters=params.get(
            "gene_expression_processing_parameters", {}
        ),
        # device=torch.device("cpu"),
        # device=torch.device(params.get("dataset_device", "cpu")),
        iterate_dataset=False,
    )
    train_loader = torch.utils.data.DataLoader(
        dataset=train_dataset,
        batch_size=params["batch_size"],
        shuffle=True,
        drop_last=True,
        # num_workers=params.get("num_workers", 0),
        num_workers=args.cpu,
    )
    
    # # Create valid dataset
    # valid_dataset = DrugSensitivityDataset(
    #     drug_sensitivity_filepath=valid_sensitivity_filepath,
    #     column_names=[args.col_drug, args.col_cell, args.col_ic50],
    #     smi_filepath=smi_filepath,
    #     gene_expression_filepath=gep_filepath,
    #     smiles_language=smiles_language,
    #     gene_list=gene_list,
    #     drug_sensitivity_min_max=params.get("drug_sensitivity_min_max", True),
    #     drug_sensitivity_processing_parameters=params.get(
    #         "drug_sensitivity_processing_parameters",
    #         train_dataset.drug_sensitivity_processing_parameters,
    #     ),
    #     gene_expression_standardize=params.get("gene_expression_standardize", True),
    #     gene_expression_min_max=params.get("gene_expression_min_max", False),
    #     gene_expression_processing_parameters=params.get(
    #         "gene_expression_processing_parameters",
    #         train_dataset.gene_expression_dataset.processing,
    #     ),
    #     # device=torch.device("cpu"),
    #     # device=torch.device(params.get("dataset_device", "cpu")),
    #     iterate_dataset=False,
    # )
    # valid_loader = torch.utils.data.DataLoader(
    #     dataset=valid_dataset,
    #     batch_size=params["batch_size"],
    #     shuffle=False,
    #     drop_last=True,
    #     # num_workers=params.get("num_workers", 0),
    #     num_workers=args.cpu,
    # )
    
    test_dataset = DrugSensitivityDataset(
        drug_sensitivity_filepath=test_sensitivity_filepath,
        column_names=[args.col_drug, args.col_cell, args.col_ic50],
        smi_filepath=smi_filepath,
        gene_expression_filepath=gep_filepath,
        smiles_language=smiles_language,
        gene_list=gene_list,
        drug_sensitivity_min_max=params.get("drug_sensitivity_min_max", True),
        drug_sensitivity_processing_parameters=params.get(
            "drug_sensitivity_processing_parameters",
            train_dataset.drug_sensitivity_processing_parameters,
        ),
        gene_expression_standardize=params.get("gene_expression_standardize", True),
        gene_expression_min_max=params.get("gene_expression_min_max", False),
        gene_expression_processing_parameters=params.get(
            "gene_expression_processing_parameters",
            train_dataset.gene_expression_dataset.processing,
        ),
        # device=torch.device("cpu"),
        # device=torch.device(params.get("dataset_device", "cpu")),
        iterate_dataset=False,
    )
    test_loader = torch.utils.data.DataLoader(
        dataset=test_dataset,
        # batch_size=params["batch_size"],
        batch_size=1024,
        shuffle=False,
        drop_last=False,
        # num_workers=params.get("num_workers", 0),
        num_workers=args.cpu,
    )
    
    logger.info(
        f"Training dataset has {len(train_dataset)} samples, test set has "
        f"{len(test_dataset)}."
    )

    # device = get_device()
    device = args.device
    
    logger.info(
        f"Device for data loader is {device} and for "
        f"model is {device}"
    )
    # save_top_model = os.path.join(model_dir, "weights/{}_{}_{}.pt")
    params.update(
        {  # yapf: disable
            "number_of_genes": len(gene_list),
            "smiles_vocabulary_size": smiles_language.number_of_tokens,
            "drug_sensitivity_processing_parameters": train_dataset.drug_sensitivity_processing_parameters,
            "gene_expression_processing_parameters": train_dataset.gene_expression_dataset.processing,
        }
    )
    model_name = params.get("model_fn", "paccmann_v2")
    model = MODEL_FACTORY[model_name](params).to(device)
    model._associate_language(smiles_language)

    # if os.path.isfile(os.path.join(model_dir, "weights", f"best_mse_{model_name}.pt")):
    if os.path.isfile(args.dir_param) :
        # logger.info("Found existing model, restoring now...")
        # model.load(os.path.join(model_dir, "weights", f"best_mse_{model_name}.pt"))
        model.load_state_dict(checkpoint['model_state_dict'])

        # with open(os.path.join(model_dir, "results", "mse.json"), "r") as f:
        #     info = json.load(f)
        #     
        #     min_rmse = info["best_rmse"]
        #     max_pearson = info["best_pearson"]
        #     min_loss = info["test_loss"]
        
        min_loss = params["test_loss"]
        min_rmse = params["best_rmse"]
        max_pearson = params["best_pearson"]
        early_stop_count = params["early_stop_count"]

    else:
        params["epoch"] = 0
        early_stop_count = 0
        min_loss, min_rmse, max_pearson = 100, 1000, 0
    
    # Define optimizer
    optimizer = OPTIMIZER_FACTORY[params.get("optimizer", "Adam")](
        model.parameters(), lr=params.get("lr", 0.01)
    )
    
    if os.path.isfile(args.dir_param) :
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    
    num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    params.update({"number_of_parameters": num_params})
    logger.info(f"Number of parameters {num_params}")
    # logger.info(model)

    # Overwrite params.json file with updated parameters.
    # with open(os.path.join(model_dir, "model_params.json"), "w") as fp:
    #     json.dump(params, fp)
    with open(args.dir_hparam, "w") as fp :
        json.dump(params, fp, indent=4)

    # Start training
    logger.info("Training about to start...\n")
    t = time()

    # model.save(save_top_model.format("epoch", "0", model_name))
    
    
    # for epoch in range(params["epochs"]):
    if not test_only :
        for epoch in range(params["epoch"]+1, params["epochs"]+1):
    
            model.train()
            logger.info(params_filepath.split("/")[-1])
            logger.info(f"== Epoch [{epoch}/{params['epochs']}] ==")
            train_loss = 0
            
            for ind, (smiles, gep, y) in enumerate(train_loader):
                y_hat, pred_dict = model(
                    torch.squeeze(smiles.to(device)), gep.to(device)
                )
                loss = model.loss(y_hat, y.to(device))
                optimizer.zero_grad()
                loss.backward()
                # Apply gradient clipping
                # torch.nn.utils.clip_grad_norm_(model.parameters(),1e-6)
                optimizer.step()
                train_loss += loss.item()
            
            logger.info(
                "\t **** TRAINING ****   "
                f"Epoch [{epoch + 1}/{params['epochs']}], "
                f"loss: {train_loss / len(train_loader):.5f}. "
                f"This took {time() - t:.1f} secs."
            )
            t = time()
            
            # Measure validation performance
            model.eval()
            with torch.no_grad():
                test_loss = 0
                predictions = []
                labels = []
                for ind, (smiles, gep, y) in enumerate(test_loader):
                # for ind, (smiles, gep, y) in enumerate(valid_loader):
                    y_hat, pred_dict = model(
                        torch.squeeze(smiles.to(device)), gep.to(device)
                    )
                    predictions.append(y_hat)
                    labels.append(y)
                    loss = model.loss(y_hat, y.to(device))
                    test_loss += loss.item()
            
            # predictions = np.array([p.cpu() for preds in predictions for p in preds])
            # labels = np.array([l.cpu() for label in labels for l in label])
            # test_pearson_a = pearsonr(torch.Tensor(predictions), torch.Tensor(labels))
            
            predictions = np.array([p.cpu().numpy() for preds in predictions for p in preds])
            labels = np.array([l.cpu().numpy() for label in labels for l in label])
            test_pearson_a = pearsonr(torch.Tensor(predictions).squeeze(), torch.Tensor(labels).squeeze())
            test_rmse_a = np.sqrt(np.mean((predictions - labels) ** 2))
            test_loss_a = test_loss / len(test_loader)
            # test_loss_a = test_loss / len(valid_loader)
            logger.info(
                f"\t **** VALIDATION **** Epoch [{epoch + 1}/{params['epochs']}], "
                f"loss: {test_loss_a:.5f}, "
                f"Pearson: {test_pearson_a:.3f}, "
                f"RMSE: {test_rmse_a:.3f}"
            )
    
            # def save(path, metric, typ, val=None):
            #     model.save(path.format(typ, metric, model_name))
            #     with open(os.path.join(model_dir, "results", metric + ".json"), "w") as f:
            #         json.dump(info, f)
            #     np.save(
            #         os.path.join(model_dir, "results", metric + "_preds.npy"),
            #         np.vstack([predictions, labels]),
            #     )
            #     if typ == "best":
            #         logger.info(
            #             f'\t New best performance in "{metric}"'
            #             f" with value : {val:.7f} in epoch: {epoch}"
            #         )
            # 
            # def update_info():
            #     return {
            #         "best_rmse": str(min_rmse),
            #         "best_pearson": str(float(max_pearson)),
            #         "test_loss": str(min_loss),
            #         "predictions": [float(p) for p in predictions],
            #     }
            
            def update_info(params, dir_hparam) :
                
                perf = {
                  	"epoch": epoch,
                  	"test_loss": float(min_loss),
                  	"best_rmse": float(min_rmse),
                  	"test_loss_a": float(test_loss_a),
                  	"best_pearson": float(max_pearson),
                  	"early_stop_count" : early_stop_count,
                }
                
                params.update(perf)
                with open(dir_hparam, "w") as fp :
                    json.dump(params, fp, indent=4)
                
                return params
    
            if test_loss_a < min_loss:
                min_rmse = test_rmse_a
                min_loss = test_loss_a
                min_loss_pearson = test_pearson_a
                # info = update_info()
                # save(save_top_model, "mse", "best", min_loss)
                
                ep_loss = epoch
                early_stop_count = 0
                params = update_info(params, args.dir_hparam)
                
                torch.save({
                    "epoch" : epoch,
                    "min_loss" : min_loss,
                    "test_loss_a" : test_loss_a,
                    "model_state_dict" : model.state_dict(), 
                    "optimizer_state_dict": optimizer.state_dict()
                }, args.dir_param)
                
            else :
                early_stop_count += 1
                print("Early stopping count {} at Epoch {}...".format(early_stop_count, epoch))
                params = update_info(params, args.dir_hparam)
                if early_stop_count>=params["patience"] :
                    print("Early stopping occured at Epoch {}...".format(epoch))
                    break
                
            if test_pearson_a > max_pearson:
                max_pearson = test_pearson_a
                max_pearson_loss = test_loss_a
                # info = update_info()
                # save(save_top_model, "pearson", "best", max_pearson)
                ep_pearson = epoch
                params = update_info(params, args.dir_hparam)
            # if (epoch + 1) % params.get("save_model", 100) == 0:
            #     save(save_top_model, "epoch", str(epoch))
    
    
    # Test Prediction
    model.eval()
    with torch.no_grad():
        test_loss = 0
        predictions = []
        labels = []
        for ind, (smiles, gep, y) in enumerate(test_loader):
            y_hat, pred_dict = model(
                torch.squeeze(smiles.to(device)), gep.to(device)
            )
            predictions.append(y_hat)
            labels.append(y)
            loss = model.loss(y_hat, y.to(device))
            test_loss += loss.item()
            
    
    # predictions = np.array([p.cpu() for preds in predictions for p in preds])
    # labels = np.array([l.cpu() for label in labels for l in label])
    predictions = np.array([p.cpu().numpy() for preds in predictions for p in preds])
    labels = np.array([l.cpu().numpy() for label in labels for l in label])
        
    # test_pearson_a = pearsonr(torch.Tensor(predictions), torch.Tensor(labels))
    test_pearson_a = pearsonr(torch.Tensor(predictions.squeeze(1)), torch.Tensor(labels.squeeze(1)))
    test_rmse_a = np.sqrt(np.mean((predictions - labels) ** 2))
    test_loss_a = test_loss / len(test_loader)
    logger.info(
        f"\t **** TESTING ****"
        f"loss: {test_loss_a:.5f}, "
        f"Pearson: {test_pearson_a:.3f}, "
        f"RMSE: {test_rmse_a:.3f}"
    )
    
    if not test_only :
        logger.info(
            "Overall best performances are: \n \t"
            f"Loss = {min_loss:.4f} in epoch {ep_loss} "
            f"\t (Pearson was {min_loss_pearson:4f}) \n \t"
            f"Pearson = {max_pearson:.4f} in epoch {ep_pearson} "
            f"\t (Loss was {max_pearson_loss:2f})"
        )
    # save(save_top_model, "training", "done")
    pred_to_csv(predictions, ic50_test, args.dir_test)
    logger.info("Done with training, models saved, shutting down.")


if __name__ == "__main__":
    # parse arguments
    args = parser.parse_args()
    # # run the training
    # main(
    #     args.train_sensitivity_filepath,
    #     args.test_sensitivity_filepath,
    #     args.gep_filepath,
    #     args.smi_filepath,
    #     args.gene_filepath,
    #     args.smiles_language_filepath,
    #     args.model_path,
    #     args.params_filepath,
    #     args.training_name,
    # )
    
    args.gep_filepath = "_data/EXP.csv"
    args.smi_filepath = "_data/SMILES_GDSC.smi"
    args.gene_filepath = "_data/2093_genes.pkl"
    args.smiles_language_filepath = "single_pytorch_model/smiles_language"
    args.params_filepath = "single_pytorch_model/model_params.json"

    nth = args.nth
    choice = args.choice
    args.device = "cuda" if torch.cuda.is_available() else "cpu"
    
    args.model_filepath = ""
    args = create_dir_out_(args, suf_param="pt", retrain=True)
    ic50_data = pd.read_csv(args.ic50, header=0)
    cell_data = pd.read_csv(args.gep_filepath, index_col=0)
    drug_data = pd.read_csv(args.smi_filepath, index_col=1, sep="\t", header=None)
    
    ic50_data = filt_ic50_(ic50_data, cell_data, drug_data, args)
    ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, retrain=True)
    args = ic50_split_path(args, ic50_train, None, ic50_test)
    # ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
    # args = ic50_split_path(args, ic50_train, ic50_valid, ic50_test)
    
    args.epochs = 300
    args.patience = 10
    args.train_sensitivity_filepath = args.dir_train_temp
    # args.valid_sensitivity_filepath = args.dir_valid_temp
    args.test_sensitivity_filepath = args.dir_test_temp
    
    main(
        args.train_sensitivity_filepath,
        # args.valid_sensitivity_filepath,
        args.test_sensitivity_filepath,
        args.gep_filepath,
        args.smi_filepath,
        args.gene_filepath,
        args.smiles_language_filepath,
        args.params_filepath,
        args
    )

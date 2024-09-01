import json
import torch
from pytoda.smiles.smiles_language import SMILESTokenizer
from pytoda.datasets import SMILESTokenizerDataset


def process_drug(smi_filepath) :
    padding: bool = True;
    padding_length: int = None;
    add_start_and_stop: bool = False;
    augment: bool = False;
    canonical: bool = False;
    kekulize: bool = False;
    all_bonds_explicit: bool = False;
    all_hs_explicit: bool = False;
    randomize: bool = False;
    remove_bonddir: bool = False;
    remove_chirality: bool = False;
    selfies: bool = False;
    sanitize: bool = True;
    vocab_file: str = None;
    iterate_dataset: bool = True;
    device: torch.device = torch.device(params.get("dataset_device", "cpu"));
    backend: str = 'eager';
    
    smiles_dataset = SMILESTokenizerDataset(
        smi_filepath,
        smiles_language=smiles_language,
        augment=augment,
        canonical=canonical,
        kekulize=kekulize,
        all_bonds_explicit=all_bonds_explicit,
        all_hs_explicit=all_hs_explicit,
        remove_bonddir=remove_bonddir,
        remove_chirality=remove_chirality,
        selfies=selfies,
        sanitize=sanitize,
        randomize=randomize,
        padding=padding,
        padding_length=padding_length,
        add_start_and_stop=add_start_and_stop,
        device=device,
        vocab_file=vocab_file,
        iterate_dataset=iterate_dataset,
        backend=backend,
    )
    
    return smiles_dataset

params_filepath = "single_pytorch_model/model_params.json"

params = {}
with open(params_filepath) as fp:
    params.update(json.load(fp))

smiles_language_filepath = "single_pytorch_model/smiles_language"
smiles_language = SMILESTokenizer(device=torch.device("cpu")).from_pretrained(smiles_language_filepath)
smiles_language.set_encoding_transforms(
    add_start_and_stop=params.get("add_start_and_stop", True),
    padding=params.get("padding", True),
    padding_length=params.get("smiles_padding_length", None),
)
# smiles_language.device = torch.device("cpu")

# smi_filepath: str = "data/smiles/gdsc.smi";
smi_filepath: str = "_data/SMILES_GDSC.smi";
smiles_dataset = process_drug(smi_filepath)

import pdb
pdb.set_trace()



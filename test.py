#!/usr/bin/env python

import os
import re
import sys
import pickle
import numpy as np
import pandas as pd

import torch
from torch.optim import Adam

import torch_geometric
from model.GCNPath_Plain import *
from model.GCNPath_Attn_v1 import *
from model.GCNPath_Attn_v2 import *
from model.GCNPath_SSI_v1 import *
from model.GCNPath_SSI_v2 import *

from utils.utils import *
from utils.utils_model import *


def parse_parameter() :
    import argparse
    parser = argparse.ArgumentParser()
    
    cell_ = "processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle"
    drug_ = "processed/drug_data/CCLE_Drug_Custom.pickle"
    ic50_ = "data/ic50_data/IC50_CCLE.txt"
    
    dir_model = "results/IC50_GDSC/Normal/RGCN"
    dir_param_ = "{}/param_0.pt".format(dir_model)
    dir_hparam_ = "{}/hyper_param_0.pickle".format(dir_model)
    out_file_ = "{}/pred_ccle_0.csv".format(dir_model)
    
    parser.add_argument("-dir_param", type=str, default=dir_param_, help="Model parameter (Default: {})".format(dir_param_))
    parser.add_argument("-dir_hparam", type=str, default=dir_hparam_, help="Model hyper-parameter (Default: {})".format(dir_hparam_))
    parser.add_argument("-out_file", type=str, default=out_file_, help="Output file in CSV format, created if non-existig directory (Default: {})".format(out_file_))
    parser.add_argument("-out_grad_cam", type=str, default="", help="Pathway importance scores calculated with Grad-CAM in CSV format, created if non-existig directory (Default: None)")
    
    parser.add_argument("-cell", type=str, default=cell_, help="Cell data in Pickle format (Default: {})".format(cell_))
    parser.add_argument("-drug", type=str, default=drug_, help="Drug data in Pickle format (Default: {})".format(drug_))
    parser.add_argument("-ic50", type=str, default=ic50_, help="IC50 data in CSV or TXT format (Default: {})".format(ic50_))
    
    parser.add_argument("-col_cell", type=str, default="Cell_BROAD_ID", help="Cell column in IC50 data (Default: Cell_BROAD_ID)")
    parser.add_argument("-col_drug", type=str, default="Drug_CID", help="Drug column in IC50 data (Default: Drug_CID)")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50", help="IC50 column in IC50 data (Default: LN_IC50)")
    
    parser.add_argument("-bs", type=int, default=256, help="Batch size (Default: 256)")
    parser.add_argument("-cpu", type=int, default=0, help="CPUs for subprocess workers (Default: 0)")
    parser.add_argument("-pin_memory", type=int, default=0, help="Pin memory for fast data loader in GPU training (Default: 0)")
    
    return parser.parse_args()


args = parse_parameter()
dir_param = args.dir_param     # Model Parameters
dir_hparam = args.dir_hparam   # Model Hyper-Parameters
out_file = args.out_file       # Output File
out_cam = args.out_grad_cam    # Output File

input_cell = args.cell         # Cell Data [CCLE]
input_drug = args.drug         # Drug Data [CCLE]
input_ic50 = args.ic50         # IC50 Data [CCLE]

batch_size = args.bs           # Batch Size [256]
num_workers = args.cpu         # Num of Subprocess [0]
do_grad_cam = out_cam not in ["", "None", None]


# Set output directory
breakpoint()
labels = args.col_ic50!="0"
out_dir = "/".join(out_file.split("/")[:-1])
os.makedirs(out_dir, exist_ok=True)
print("\n### Save a prediction file into {}...".format(out_file))

if do_grad_cam :
    out_dir_ = "/".join(out_cam.split("/")[:-1])
    os.makedirs(out_dir_, exist_ok=True)
    print("\n### Save a pathway importance file into {}...".format(out_cam))


# IC50 Filtering & DataLoader Preparation
print("\n### Data Processing")

ext_ic50 = input_ic50.split(".")[-1]
sep_ic50 = "," if ext_ic50=="csv" else "\t"
ic50_data = pd.read_csv(input_ic50, header=0, sep=sep_ic50)

with open(dir_hparam, "rb") as f :
    args_hparam = pickle.load(f)

with open(input_cell, "rb") as f :
    cell_data, path_dict, _, _ = pickle.load(f)
    cell_data = {str(k):cell_data[k] for k in cell_data.keys()}

with open(input_drug, "rb") as f :
    drug_data = pickle.load(f)
    drug_data = {str(k):drug_data[k] for k in drug_data.keys()}

device = "cuda" if torch.cuda.is_available() else "cpu"
ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)

args.device = device
no_labels = (labels==0)
args.pin_memory = args.pin_memory!=0 and device!="cpu"
 
args.gnn_cell = args_hparam["gnn_cell"]
args.gnn_drug = args_hparam["gnn_drug"]
args.attn_mode = args_hparam["attn_mode"]
 
test_loader = load_data(ic50_data, cell_data, drug_data, args, 
                        no_labels=no_labels, batch_size=batch_size, 
                        num_workers=num_workers, shuffle=False)


# Test using the pre-trained model
GCNPath = [GCNPath, GCNPath_Attn_v1, GCNPath_Attn_v2, GCNPath_SSI_v1, GCNPath_SSI_v2][args.attn_mode]
model = GCNPath(**args_hparam)
checkpoint = torch.load(dir_param, map_location=device)
model.load_state_dict(checkpoint['model_state_dict'])

verbose = True
model.summary(verbose=verbose)

print("\n### Data Prediction...")
test_ = test if labels else test_no_labels
pred_test, test_time = test_(model, test_loader, device=device, return_attn=False)
pred_to_csv(pred_test, ic50_data, out_file)

# Time Calculation
time_to_csv(test_time=test_time, dir_time=args.dir_time)

if do_grad_cam :
    print("\n### Get pathway importance scores [Grad-CAM]")
    col = [args.col_cell, args.col_drug]
    importance = grad_cam(model, cell_data, drug_data, ic50_data, args)
    importance = pd.DataFrame(importance, columns=list(path_dict.keys()))
    pred_to_csv(importance, ic50_data.loc[:, col], out_cam, grad_cam=True)

print("\n### Prediction Completed!!!")

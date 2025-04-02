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
    
    cell_ = "processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle"
    drug_ = "processed/drug_data/GDSC_Drug_Custom.pickle"
    ic50_ = "data/ic50_data/IC50_GDSC2.txt"
    out_dir_ = "results/example"
    
    parser.add_argument("-nth", type=int, default=0, help="Train fold (Default: 0, 0-24 in Strict-Blind tests, 0-9 in the rest tests)")
    parser.add_argument("-choice", type=int, default=0, help="Test type (Default: 0, 0: Normal, 1: Cell-Blind, 2: Drug-Blind, 3: Strict-Blind)")
    parser.add_argument("-out_dir", type=str, default=out_dir_, help="Output directory, created if non-existig directory (Default: {})".format(out_dir_))
    
    parser.add_argument("-cell", type=str, default=cell_, help="Cell data in Pickle format (Default: {})".format(cell_))
    parser.add_argument("-drug", type=str, default=drug_, help="Drug data in Pickle format (Default: {})".format(drug_))
    parser.add_argument("-ic50", type=str, default=ic50_, help="IC50 data in CSV or TXT format (Default: {})".format(ic50_))
    
    parser.add_argument("-col_cell", type=str, default="Cell", help="Cell column in IC50 data (Default: Cell)")
    parser.add_argument("-col_drug", type=str, default="Drug", help="Drug column in IC50 data (Default: Drug)")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50", help="IC50 column in IC50 data (Default: LN_IC50)")
    
    parser.add_argument("-gnn_cell", type=int, default=4, help="Cell GCN mode (Default: 4, 0: None, 1: GCN, 2: GAT, 3: GIN, 4: RGCN, 5: RGAT)")
    parser.add_argument("-gnn_drug", type=int, default=2, help="Drug GCN mode (Default: 2, 0: None, 1: GCN, 2: GAT, 3: GIN, 4: GINE)")
    parser.add_argument("-mode_cell", type=int, default=3, help="""Cell stacking mode (Default: 3, 0: plain, 1: res [ResNet], 2: res+ [ResNet+], 
                                                                3: dense [DenseNet], 4: jk_cat [Jumping Knowledge, Concat], 5: jk_max [Jumping Knowledge, Max Pooling])""")
    parser.add_argument("-mode_drug", type=int, default=3, help="""Drug stacking mode (Default: 3, 0: plain, 1: res [ResNet], 2: res+ [ResNet+], 
                                                                3: dense [DenseNet], 4: jk_cat [Jumping Knowledge, Concat], 5: jk_max [Jumping Knowledge, Max Pooling])""")
    
    parser.add_argument("-dim_cell", type=int, default=8, help="Cell hidden dimension (Default: 16)")
    parser.add_argument("-dim_drug", type=int, default=128, help="Drug hidden dimension (Default: 128)")
    parser.add_argument("-dim_pred", type=int, default=512, help="Prediction hidden dimension (Default: 512)")
    
    parser.add_argument("-dim_compact", type=int, default=32, help="Cell dim_compact dimension (Default: min(32, dim_cell_final))")
    parser.add_argument("-dim_embed_cell", type=int, default=256, help="Cell embedding dimension before Prediction module (Default: 256)")
    parser.add_argument("-dim_embed_drug", type=int, default=256, help="Drug embedding dimension before Prediction module (Default: 256)")
    
    parser.add_argument("-n_hid_cell", type=int, default=3, help="Cell layer number (Default: 3)")
    parser.add_argument("-n_hid_drug", type=int, default=3, help="Drug layer number (Default: 3)")
    parser.add_argument("-n_hid_pred", type=int, default=2, help="Prediction layer number (Default: 2)")
    
    parser.add_argument("-h_gat_cell", type=int, default=4, help="GAT heads in Cell layer (Default: 4)")
    parser.add_argument("-h_gat_drug", type=int, default=4, help="GAT heads in Drug layer (Default: 4)")
    parser.add_argument("-concat_gat", type=bool, default=True, help="The concat mode of GAT heads (Default: True, True: Concat, False: Average)")
    
    parser.add_argument("-attn_mode", type=int, default=0, help="Attention mode (Default: 0, 0: No attn, 1: Attn_v1, 2: Attn_v2, 3: SSI-DDI+Attn, 4: SA-DDI+Attn)")
    parser.add_argument("-dim_attn", type=int, default=32, help="Attention hidden dimension (Default: 32)")
    parser.add_argument("-coef_ffnn", type=int, default=4, help="Attention hidden dimension ratio during FFNN (Default: 4 [128=32x4])")
    parser.add_argument("-n_attn_layer", type=int, default=1, help="Attention layer number (Default: 1)")
    parser.add_argument("-h_attn", type=int, default=4, help="Attention heads (Default: 4)")
    
    parser.add_argument("-act", type=int, default=3, help="Activation (Default: 3, 0: ReLU, 1: ELU_0.1, 2: ELU_0.5, 3: ELU_1.0, 4: LeakyReLU_0.1, 5: LeakyReLU_0.2)")
    parser.add_argument("-norm", type=str, default="batch", help="Normalization (Default: batch, batch: BatchNorm1d, layer: LayerNorm, instance: InstanceNorm1d, None/False: None)")
    parser.add_argument("-drop", type=float, default=0.2, help="Dropout ratio (Default: 0.2)")
    
    parser.add_argument("-e", type=int, default=300, help="Epochs (Default: 300)")
    parser.add_argument("-p", type=int, default=10, help="Early stopping (Default: 10)")
    parser.add_argument("-sch", type=int, default=0, help="Training scheduler (Default: 0, 0: None, 1: Exponential, 2: Cosine)")
    
    parser.add_argument("-bs", type=int, default=256, help="Batch size (Default: 256)")
    parser.add_argument("-lr", type=float, default=0.001, help="Learning rate (Default: 0.001)")
    parser.add_argument("-cpu", type=int, default=0, help="CPUs for subprocess workers (Default: 0)")
        
    parser.add_argument("-seed_model", type=int, default=2021, help="Seed for model parameters (Default: 2021)")
    parser.add_argument("-seed_split", type=int, default=2021, help="Seed for split of IC50 Data (Default: 2021)")
    parser.add_argument("-pin_memory", type=int, default=0, help="Pin memory for fast data loader in GPU training (Default: 0)")
    
    return parser.parse_args()


args = parse_parameter()
nth = args.nth            # Train Fold [0-24 for S_Blind, 0-9 for the rest]
choice = args.choice      # Test Type [0-3, Normal, C_Blind, D_Blind, S_Blind]

input_cell = args.cell    # Cell Data [SANGER Pasports]
input_drug = args.drug    # Drug Data [GDSC]
input_ic50 = args.ic50    # IC50 Data [GDSC]
out_dir = args.out_dir    # Output Directory

epochs = args.e           # Epochs [300]
patience = args.p         # Patience [10]
scheduler = args.sch      # Scheduler [0-2]

batch_size = args.bs      # Batch Size [256]
num_workers = args.cpu    # Num of Subprocess [0]
learning_rate = args.lr   # Learning Rate [0.001]


# Set seed for the model
set_random_seed(args.seed_model)

print("\n### Input File List...")
print("# Cell : {}".format(input_cell))
print("# Drug : {}".format(input_drug))
print("# IC50 : {}".format(input_ic50))

# Set parameters in args
args = choose_params(args)   
# Activation, GNN_Cell, GNN_Drug, Mode_Cell, Mode_Drug
args = create_dir_out(args, dir_out=out_dir, suf_param="pt", retrain=False)


print("\n### Data Processing")
ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")

with open(input_cell, "rb") as f :
    cell_data, _, _, _ = pickle.load(f)
    cell_data = {str(k):cell_data[k] for k in cell_data.keys()}

with open(input_drug, "rb") as f :
    drug_data = pickle.load(f)
    drug_data = {str(k):drug_data[k] for k in drug_data.keys()}

device = "cuda" if torch.cuda.is_available() else "cpu"
ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, seed=args.seed_split, retrain=False)

no_labels = False
args.pin_memory = args.pin_memory!=0 and device!="cpu"

train_loader = load_data(ic50_train, cell_data, drug_data, args, no_labels, batch_size, num_workers, shuffle=True, fix_seed=True)
valid_loader = load_data(ic50_valid, cell_data, drug_data, args, no_labels, batch_size, num_workers, shuffle=False, fix_seed=True)
test_loader = load_data(ic50_test, cell_data, drug_data, args, no_labels, batch_size, num_workers, shuffle=False, fix_seed=True)

# Deep Learning Module & Dimension
args.dim_cell = args.dim_cell if args.gnn_cell!="linear" else 256
args.dim_drug = args.dim_drug if args.gnn_drug!="linear" else 256

sample_cell = list(cell_data.values())[0]
sample_drug = list(drug_data.values())[0]
args = get_cell_info(args, sample_cell)   # [1, 8, 8, 8]
args = get_drug_info(args, sample_drug)   # [85, 128, 128, 128]
args = get_pred_info(args)                # [512, 512]

# Train Model
args_hparam = vars(args)
param_names = ["gnn_cell", "gnn_drug", "mode_cell", "mode_drug", 
               "dim_cell", "dim_drug", "dim_pred", "act", "norm", "drop",
               "n_drug_edge", "concat_gat", "h_gat_cell", "h_gat_drug", "n_rel", 
               "attn_mode", "dim_attn", "h_attn", "dim_embed_cell", "dim_embed_drug", 
               "dim_compact", "n_path", "coef_ffnn", "n_attn_layer"]

args_hparam = {k:args_hparam[k] for k in param_names}
GCNPath = [GCNPath, GCNPath_Attn_v1, GCNPath_Attn_v2, GCNPath_SSI_v1, GCNPath_SSI_v2][args.attn_mode]
model = GCNPath(**args_hparam)

# Save hyper-parameters
with open(args.dir_hparam, "wb") as f :
    pickle.dump(args_hparam, f)
print("\n### Model hyper-parameters saved!!!")

# model = torch.compile(model)
# torch._dynamo.config.verbose=True
# torch._dynamo.config.suppress_errors = True

verbose = True
model.summary(verbose=verbose)
optimizer = Adam(model.parameters(), lr=learning_rate)
scheduler = choose_scheduler(scheduler)   # 0 [None]

print("\n### Train Start!!!")
model, train_time_list = train(model, train_loader, valid_loader, optimizer=optimizer, device=device, 
                               epochs=epochs, patience=patience, scheduler=scheduler, 
                               dir_param=args.dir_param, dir_log=args.dir_log)

# Valid Prediction
print("\n### Valid Performance...")
pred_valid, valid_time = test(model, valid_loader, device=device, return_attn=False)
pred_to_csv(pred_valid, ic50_valid, args.dir_valid)

# Test Prediction
print("\n### Test Performance...")
pred_test, test_time = test(model, test_loader, device=device, return_attn=False)
pred_to_csv(pred_test, ic50_test, args.dir_test)

# Time Calculation
time_to_csv(train_time_list, test_time, valid_time, dir_time=args.dir_time)

print("\n### Train & Test Completed!!!")

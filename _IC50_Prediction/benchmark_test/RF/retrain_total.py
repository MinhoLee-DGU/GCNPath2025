#!/usr/bin/env python

import os
import re
import sys
import pickle
import numpy as np
import pandas as pd

from utils.utils import *

import joblib
from sklearn.ensemble import RandomForestRegressor


def load_data_classic(ic50_data, cell_data, drug_data, args, no_labels=False):
    cell_name = ic50_data[args.col_cell]
    drug_name = ic50_data[args.col_drug]
    
    cell_data_ = np.vstack([cell_data[i] for i in cell_name])
    drug_data_ = np.vstack([drug_data[i] for i in drug_name])
    data_y = ic50_data[args.col_ic50] if not no_labels else None
    
    data_x = np.concatenate([cell_data_, drug_data_], axis=1)
    return data_x, data_y


def parse_parameter() :
    import argparse
    parser = argparse.ArgumentParser()
    
    cell_ = "_data/SANGER_RNA_Lin.pickle"
    drug_ = "_data/GDSC_Drug_Morgan.pickle"
    ic50_ = "_data/IC50_GDSC2.txt"
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
    
    parser.add_argument("-cpu", type=int, default=1, help="CPUs for training (Default: 1)")
    parser.add_argument("-seed_model", type=int, default=2021, help="Seed for model parameters (Default: 2021)")
    parser.add_argument("-seed_split", type=int, default=2021, help="Seed for split of IC50 Data (Default: 2021)")
    
    return parser.parse_args()


args = parse_parameter()
input_cell = args.cell    # Cell Data [SANGER Pasports]
input_drug = args.drug    # Drug Data [GDSC]
input_ic50 = args.ic50    # IC50 Data [GDSC]
out_dir = args.out_dir    # Output Directory


print("\n### Input File List...")
print("# Cell : {}".format(input_cell))
print("# Drug : {}".format(input_drug))
print("# IC50 : {}".format(input_ic50))
args = create_dir_out(args, dir_out=out_dir, suf_param="joblib", retrain=True)
args.dir_param = "{}/param_retrain_seed{}.joblib".format(args.dir_out, args.seed_model)


print("\n### Data Processing")
ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")

with open(input_cell, "rb") as f :
    cell_data = pickle.load(f)
    cell_data = {str(k):cell_data[k] for k in cell_data.keys()}

with open(input_drug, "rb") as f :
    drug_data = pickle.load(f)
    drug_data = {str(k):drug_data[k] for k in drug_data.keys()}

ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
train_x, train_y = load_data_classic(ic50_data, cell_data, drug_data, args, no_labels=False)


# Train Model
# min_samples_split ranges from [2, inf) (1 raised error...)
model = RandomForestRegressor(n_jobs=args.cpu, random_state=args.seed_model, 
                              n_estimators=100, criterion="squared_error", 
                              min_samples_split=2, min_samples_leaf=1, max_features=1.0)

print("\n### Train Start!!!")
model = model.fit(train_x, train_y)

# Save Model
print("\n### Save Model...")
joblib.dump(model, args.dir_param)

print("\n### Train Completed!!!")

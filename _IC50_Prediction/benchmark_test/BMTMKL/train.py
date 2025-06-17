#!/usr/bin/env python

import os
import re
import sys
import pickle
import numpy as np
import pandas as pd

from utils.utils import *


def parse_parameter() :
    import argparse
    parser = argparse.ArgumentParser()
    
    cell_ = "_data/cell_names.txt"
    ic50_ = "_data/IC50_GDSC2.txt"
    out_dir_ = "results/example"
    
    parser.add_argument("-nth", type=int, default=0, help="Train fold (Default: 0, 0-24 in Strict-Blind tests, 0-9 in the rest tests)")
    parser.add_argument("-choice", type=int, default=0, help="Test type (Default: 0, 0: Normal, 1: Cell-Blind, 2: Drug-Blind, 3: Strict-Blind)")
    parser.add_argument("-out_dir", type=str, default=out_dir_, help="Output directory, created if non-existig directory (Default: {})".format(out_dir_))
    
    parser.add_argument("-cell", type=str, default=cell_, help="Cell data in Pickle format (Default: {})".format(cell_))
    parser.add_argument("-ic50", type=str, default=ic50_, help="IC50 data in CSV or TXT format (Default: {})".format(ic50_))
    
    parser.add_argument("-col_cell", type=str, default="Cell", help="Cell column in IC50 data (Default: Cell)")
    parser.add_argument("-col_drug", type=str, default="Drug", help="Drug column in IC50 data (Default: Drug)")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50", help="IC50 column in IC50 data (Default: LN_IC50)")
    
    parser.add_argument("-cpu", type=int, default=1, help="CPUs for training (Default: 1)")
    parser.add_argument("-seed_model", type=int, default=2021, help="Seed for model parameters (Default: 2021)")
    parser.add_argument("-seed_split", type=int, default=2021, help="Seed for split of IC50 Data (Default: 2021)")
    
    return parser.parse_args()


args = parse_parameter()
nth = args.nth            # Train Fold [0-24 for S_Blind, 0-9 for the rest]
choice = args.choice      # Test Type [0-3, Normal, C_Blind, D_Blind, S_Blind]

input_cell = args.cell    # Cell Data [SANGER Pasports]
input_ic50 = args.ic50    # IC50 Data [GDSC]
out_dir = args.out_dir    # Output Directory


print("\n### Input File List...")
print("# Cell : {}".format(input_cell))
print("# IC50 : {}".format(input_ic50))
args = create_dir_out(args, dir_out=out_dir, suf_param="joblib", retrain=False)

print("\n### Data Processing")
sep = "," if ".csv" in input_ic50 else "\t"
ic50_data = pd.read_csv(input_ic50, header=0, sep=sep)
cell_data = pd.read_csv(input_cell, header=None, sep="\t")

cell_data = cell_data.values.squeeze()
cell_data = pd.DataFrame(cell_data, index=cell_data, columns=[args.col_cell])

drug_data = ic50_data[args.col_drug].unique()
drug_data = pd.DataFrame(drug_data, index=drug_data, columns=[args.col_drug])

cell_data.index = cell_data.index.astype(str)
drug_data.index = drug_data.index.astype(str)
ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, seed=args.seed_split, retrain=True)

file_train = f"{args.out_dir}/ic50_train_{nth}.txt"
file_test = f"{args.out_dir}/ic50_test_{nth}.txt"
ic50_train.to_csv(file_train, sep="\t", index=False)
ic50_test.to_csv(file_test, sep="\t", index=False)

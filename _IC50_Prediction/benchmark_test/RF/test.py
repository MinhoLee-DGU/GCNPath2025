#!/usr/bin/env python

import os
import re
import sys
import time
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
    ic50_ = "_data/IC50_ChEMBL.txt"
    
    parser.add_argument("-cell", type=str, default=cell_, help="Cell data in Pickle format (Default: {})".format(cell_))
    parser.add_argument("-drug", type=str, default=drug_, help="Drug data in Pickle format (Default: {})".format(drug_))
    parser.add_argument("-ic50", type=str, default=ic50_, help="IC50 data in CSV or TXT format (Default: {})".format(ic50_))
    
    parser.add_argument("-col_cell", type=str, default="Cell", help="Cell column in IC50 data (Default: Cell)")
    parser.add_argument("-col_drug", type=str, default="Drug", help="Drug column in IC50 data (Default: Drug)")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50", help="IC50 column in IC50 data (Default: LN_IC50)")
    
    parser.add_argument("-cpu", type=int, default=1, help="CPUs for training (Default: 1)")
    parser.add_argument("-dir_param", type=str, required=True, help="Model parameter")
    parser.add_argument("-out_file", type=str, required=True, help="Output file in CSV format, created if non-existig directory")
    parser.add_argument("-out_time", type=str, required=True, help="Log file for inference time")
    
    return parser.parse_args()


args = parse_parameter()
input_cell = args.cell    # Cell Data [SANGER Pasports]
input_drug = args.drug    # Drug Data [GDSC]
input_ic50 = args.ic50    # IC50 Data [GDSC]


print("\n### Input File List...")
print("# Cell : {}".format(input_cell))
print("# Drug : {}".format(input_drug))
print("# IC50 : {}".format(input_ic50))


print("\n### Data Processing")
ext_ic50 = input_ic50.split(".")[-1]
sep_ic50 = "," if ext_ic50=="csv" else "\t"
ic50_data = pd.read_csv(input_ic50, header=0, sep=sep_ic50)

with open(input_cell, "rb") as f :
    cell_data = pickle.load(f)
    cell_data = {str(k):cell_data[k] for k in cell_data.keys()}

with open(input_drug, "rb") as f :
    drug_data = pickle.load(f)
    drug_data = {str(k):drug_data[k] for k in drug_data.keys()}

ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
test_x, test_y = load_data_classic(ic50_data, cell_data, drug_data, args, no_labels=True)


# Test Model
model = joblib.load(args.dir_param)
model.n_jobs = args.cpu

start = time.perf_counter()
pred_test = model.predict(test_x)
end = time.perf_counter()

test_time = 1000 * (end - start)
pred_to_csv(pred_test, ic50_data, args.out_file)

# Time Calculation
if args.out_time not in [None, "None", ""] :
    time_to_csv(test_time=test_time, dir_time=args.out_time)

print("\n### Test Completed!!!")

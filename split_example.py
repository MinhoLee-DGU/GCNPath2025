#!/usr/bin/env python

import torch
import random
import pickle
from utils.utils import *
from utils.utils_model import *
import matplotlib.pyplot as plt

import argparse
parser = argparse.ArgumentParser()
args = parser.parse_args()

args.col_cell = "Cell"
args.col_drug = "Drug"
args.col_ic50 = "LN_IC50"


def split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col="Cell") :

    num_total = ic50_data[col].value_counts()
    num_train = ic50_train[col].value_counts()
    num_valid = ic50_valid[col].value_counts()
    num_test = ic50_test[col].value_counts()
    
    num_info = pd.concat([num_total, num_train, num_valid, num_test], axis=1, join="outer")
    num_info.columns = ["N_Total", "N_Train", "N_Valid", "N_Test"]
    num_info["Ratio_Train"] = num_info["N_Train"] / num_info["N_Total"]
    num_info["Ratio_Valid"] = num_info["N_Valid"] / num_info["N_Total"]
    num_info["Ratio_Test"] = num_info["N_Test"] / num_info["N_Total"]
    
    print("Mean [Train Ratio] : {}".format(round(num_info.Ratio_Train.mean(), 3)))
    print("Mean [Valid Ratio] : {}".format(round(num_info.Ratio_Valid.mean(), 3)))
    print("Mean [Test Ratio] : {}".format(round(num_info.Ratio_Test.mean(), 3)))
    
    plt.close()
    plt.hist(num_info.Ratio_Train)
    plt.hist(num_info.Ratio_Valid)
    plt.hist(num_info.Ratio_Test)
    plt.show()


def set_random_seed(seed, deterministic=True):
    """Set random seed."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


input_ic50 = "data/ic50_data/IC50_GDSC.txt"
input_cell = "processed/cell_data/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle"
input_drug = "processed/drug_data/GDSC_Drug_Custom.pickle"

ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")

with open(input_cell, "rb") as f :
    cell_data, _, _, _ = pickle.load(f)

with open(input_drug, "rb") as f :
    drug_data = pickle.load(f)

nth = 0
choice = 0
ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)

print("\n### Normal... [seed=2021 (default)]")
ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)

# Split ratio is kept stable
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_cell)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_drug)

# Seed is set (2021 in default) so split is reproducable
print("\n### Normal... [seed=2021]")
ic50_train_21, ic50_valid_21, ic50_test_21 = split_ic50(ic50_data, args, choice=choice, nth=nth)
print(ic50_train.compare(ic50_train_21))   # Empty DataFrame [No difference between two ic50_train_]
print(ic50_valid.compare(ic50_valid_21))   # Empty DataFrame [No difference between two ic50_train_]
print(ic50_test.compare(ic50_test_21))     # Empty DataFrame [No difference between two ic50_train_]

print("\n### Normal... [seed=2023]")
ic50_train_23, ic50_valid_23, ic50_test_23 = split_ic50(ic50_data, args, choice=choice, nth=nth, seed=2023)
print(ic50_train.compare(ic50_train_23))
print(ic50_valid.compare(ic50_valid_23))
print(ic50_test.compare(ic50_test_23))

print("\n### random.seed(2023)")
# import random
# random.seed(2023)
set_random_seed(2023)

ic50_train_, ic50_valid_, ic50_test_ = split_ic50(ic50_data, args, choice=choice, nth=nth)
print(ic50_train_.compare(ic50_train_23))
print(ic50_valid_.compare(ic50_valid_23))
print(ic50_test_.compare(ic50_test_23))

print(ic50_train_.compare(ic50_train_21))
print(ic50_valid_.compare(ic50_valid_21))
print(ic50_test_.compare(ic50_test_21))

### Cell Blind
nth = 0
choice = 1
print("\n### Cell-BLind...")

ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_cell)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_drug)


### Drug Blind
nth = 0
choice = 2
print("\n### Drug-BLind...")

ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_cell)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_drug)


### Strict Blind
nth = 0
choice = 3
print("\n### Strict-BLind...")

ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_cell)
split_info(ic50_data, ic50_train, ic50_valid, ic50_test, col=args.col_drug)

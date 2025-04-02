#!/usr/bin/env python

import numpy as np
import pandas as pd

def parse_drug():
    import argparse
    parser = argparse.ArgumentParser()
    input_dir_ = "_data/SMILES_GDSC.csv"
    output_dir_="_data/SMILES_GDSC.smi"

    parser.add_argument("-input_dir", type=str, default=input_dir_)
    parser.add_argument("-output_dir", type=str, default=output_dir_)
    parser.add_argument("-col_names", type=str, default="Drug_CID")
    parser.add_argument("-col_smiles", type=str, default="SMILES_CAN")
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_drug()
    col_names = args.col_names
    col_smiles = args.col_smiles
    
    sep = args.input_dir.split('/')[(-1)].split('.')[(-1)]
    sep = ',' if sep == 'csv' else '\t'
    drug_data = pd.read_csv(args.input_dir, sep=sep)
    drug_data[[col_smiles, col_names]].to_csv(args.output_dir, header=False, index=False, sep="\t")


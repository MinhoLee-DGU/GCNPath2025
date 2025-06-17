#!/usr/bin/env python

import os
import torch
import pickle
import numpy as np
import pandas as pd

def process_cell_linear(cell_data) :
    cell_data_dict = {}
    for cell in cell_data.index :
        cell_data_tp = cell_data.loc[[cell]].values[0]
        cell_data_dict[cell] = torch.from_numpy(cell_data_tp).float()
    return cell_data_dict

def add_rand_noise(cell_data, std=1.0, noise_seed=2021) :
    if noise_seed==-1 :
        pass
    else :
        rng = np.random.default_rng(seed=noise_seed)
        noise = rng.normal(0, std, size=cell_data.shape)
        cell_data = cell_data + noise
        print(f"# Adding noise N(0, {std}) to Cell Data...")
        print(f"{noise}\n")
    return cell_data

def parse_cell() :
    import argparse
    parser = argparse.ArgumentParser()
    input_omics_ = "_data/SANGER_RNA_TPM.csv"
    output_dir_ = "_data/SANGER_RNA_Lin.pickle"
    
    parser.add_argument("-omics", type=str, default=input_omics_)
    parser.add_argument("-out", type=str, default=output_dir_)
    parser.add_argument("-col1", type=str, default="Pathway1")
    parser.add_argument("-col2", type=str, default="Pathway2")
    return parser.parse_args()


if __name__ == '__main__' :
    args = parse_cell()
    input_omics = args.omics
    output_dir = args.out
    
    path_dict = None
    cell_data = pd.read_csv(input_omics, header=0, index_col=0, sep=",")
    cell_data.columns = cell_data.columns.astype(str)
    
    noise = "Noise" in output_dir
    if noise : cell_data = add_rand_noise(cell_data)
    cell_data = process_cell_linear(cell_data)
    
    with open(output_dir, "wb") as f :
        pickle.dump(cell_data, f)

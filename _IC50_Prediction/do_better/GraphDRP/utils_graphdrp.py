#!/usr/bin/env python

import os
import re
import torch
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
from torch_geometric import data as DATA


class DatasetDef(Dataset) :
    def __init__(self, ic50_data, cell_data, drug_data, args) :
        super(DatasetDef, self).__init__()
        
        self.ic50_data = ic50_data
        self.cell_data = cell_data
        self.drug_data = drug_data
        
        self.cell_name = ic50_data[args.col_cell]
        self.drug_name = ic50_data[args.col_drug]
        self.ic50_value = ic50_data[args.col_ic50]
    
    def __len__(self):
        return len(self.ic50_data)
    
    def __getitem__(self, idx) :
        cell = self.cell_name.iloc[idx]
        drug = self.drug_name.iloc[idx]
        
        cell_data_ = self.cell_data[cell]
        drug_data_ = self.drug_data[drug]
        ic50_data_ = self.ic50_value.iloc[idx]
        
        c_size, features, edge_index = drug_data_
        GCNData = DATA.Data(x=torch.Tensor(features),
                            edge_index=torch.LongTensor(edge_index).transpose(1, 0),
                            y=torch.FloatTensor([ic50_data_]))
        
        GCNData.target = torch.FloatTensor([cell_data_])
        
        return GCNData


def filt_ic50_(ic50_data, cell_data, drug_data, args) :
    
    ic50_data.dropna(subset=[args.col_cell, args.col_drug], inplace=True)
    ic50_data = ic50_data.astype({args.col_cell:"int"})
    ic50_data = ic50_data.astype({args.col_cell:"str"})
    ic50_data = ic50_data.astype({args.col_drug:"str"})
    
    cell_list = list(cell_data.keys())
    drug_list = list(drug_data.keys())

    cell_list = list(set(ic50_data[args.col_cell]) & set(cell_list))
    drug_list = list(set(ic50_data[args.col_drug]) & set(drug_list))
    ic50_data = ic50_data[ic50_data[args.col_cell].isin(cell_list)]
    ic50_data = ic50_data[ic50_data[args.col_drug].isin(drug_list)]

    print("# [Cell] {}".format(len(cell_list)))
    print("# [Drug] {}".format(len(drug_list)))
    print("# [IC50] {}".format(ic50_data.shape[0]))

    return ic50_data


def create_dir_out_(args, dir_out=None, suf_param="pt", retrain=False) :
    
    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        if args.mode=="test" and args.pretrain==1 : dir2 = "GINConv"
        if args.mode=="test" and args.pretrain==0 : dir1 = "IC50_GDSC"
        dir_out = "Results/{}/{}".format(dir1, dir2)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
    args.dir_hparam = "{}/hyper_param_{}.pickle".format(dir_out, args.nth)
    
    if args.mode=="test" and args.pretrain==1 : 
        args.dir_test = "{}/pred_test.csv".format(dir_out)
        args.dir_perf = "{}/perf_test.csv".format(dir_out)
    elif args.mode=="test" and args.pretrain==0 : 
        args.dir_param = "{}/param_retrain_{}.{}".format(dir_out, args.nth, suf_param)
        args.dir_test = "{}/pred_{}_{}.csv".format(dir_out, args.test_name, args.nth)
        args.dir_perf = "{}/perf_{}_{}.csv".format(dir_out, args.test_name, args.nth)
    else :
        if not retrain : 
            args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
            args.dir_perf = "{}/perf_test_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_{}.{}".format(dir_out, args.nth, suf_param)
        else :
            args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_perf = "{}/perf_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_retrain_{}.{}".format(dir_out, args.nth, suf_param)
    
    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args

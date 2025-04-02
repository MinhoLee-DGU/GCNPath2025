#!/usr/bin/env python

import os
import re
import math
import numpy as np
import pandas as pd

# import keras
# from keras.utils import Sequence


class Batch():
    
    def __init__(self, ic50_data, cell_data, drug_data, args, batch_size=100, shuffle=True):
        self.batch_size = batch_size
        self.ic50_data = ic50_data
        self.drug_data = drug_data
        self.cell_data = cell_data
        
        self.offset = 0
        self.shuffle = shuffle
        self.size = len(ic50_data)
        self.indices = np.arange(len(self.ic50_data))
        
        self.col_cell = args.col_cell
        self.col_drug = args.col_drug
        self.col_ic50 = args.col_ic50

    def mini_batch(self):
        try :
            if self.offset >= self.size :
                return None
            else :
                if self.offset + self.batch_size <= self.size:
                    sub_posi = range(self.offset, self.offset + self.batch_size)
                else :
                    sub_posi = range(self.offset, self.size)
                
                self.offset += self.batch_size
                sub_posi = self.indices[sub_posi]
                drug_list = self.ic50_data[self.col_drug].iloc[sub_posi]
                cell_list = self.ic50_data[self.col_cell].iloc[sub_posi]
                ic50_list = self.ic50_data[self.col_ic50].iloc[sub_posi]
                
                drug = [self.drug_data[_] for _ in drug_list]
                cell = [self.cell_data[_] for _ in cell_list]
                value = [np.array([_]) for _ in ic50_list]
                return np.array(value), np.array(drug), np.array(cell)
              
        except :
            import pdb
            pdb.set_trace()
    
    def whole_batch(self):
        drug_list = self.ic50_data[self.col_drug].iloc[self.indices]
        cell_list = self.ic50_data[self.col_cell].iloc[self.indices]
        ic50_list = self.ic50_data[self.col_ic50].iloc[self.indices]
        
        drug = [self.drug_data[_] for _ in drug_list]
        cell = [self.cell_data[_] for _ in cell_list]
        value = [np.array([_]) for _ in ic50_list]
        
        return np.array(value), np.array(drug), np.array(cell)
    
    def reset(self):
        self.offset = 0
        if self.shuffle :
            np.random.shuffle(self.indices)
    
    def available(self):
        if self.offset < self.size:
            return True
        else:
            return False


def create_dir_out_(args, dir_out=None, suf_param="pt", retrain=False) :
    
    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        dir_out = "Results/{}/{}".format(dir1, dir2)
    
    args.dir_out = dir_out
    # os.makedirs(dir_out, exist_ok=True)
    if not os.path.exists(dir_out) : os.makedirs(dir_out)
    args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
    
    if not retrain : 
        args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
        args.dir_param = "{}/param_{}.{}".format(dir_out, args.nth, suf_param)
        args.dir_hparam = "{}/hyper_param_{}.pickle".format(dir_out, args.nth)
    else :
        args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
        args.dir_param = "{}/param_retrain_{}.{}".format(dir_out, args.nth, suf_param)
        args.dir_hparam = "{}/hyper_param_retrain_{}.pickle".format(dir_out, args.nth)
    
    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def filt_ic50_(ic50_data, cell_data, drug_data, args) :
    
    ic50_data.dropna(subset=[args.col_cell, args.col_drug], inplace=True)
    ic50_data = ic50_data.astype({args.col_cell:"int"})   # COSMIC ID [Int > Str]
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


# class DataLoader(Sequence) :
#     
#     def __init__(self, ic50_data, cell_data, drug_data, args, batch_size=100, shuffle=True) :
#         
#         self.ic50_data = ic50_data
#         self.cell_data = cell_data
#         self.drug_data = drug_data
#         self.batch_size = batch_size
#         
#         self.shuffle = shuffle
#         self.col_cell = args.col_cell
#         self.col_drug = args.col_drug
#         self.col_ic50 = args.col_ic50
#         self.on_epoch_end()
# 
#     def __len__(self) :
#         return int(np.ceil(len(self.ic50_data) / self.batch_size))
# 
#     def __getitem__(self, idx) :
#         indices = self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
#         cells = self.ic50_data[self.col_cell].iloc[indices]
#         drugs = self.ic50_data[self.col_drug].iloc[indices]
#         
#         Batch_C = self.cell_data.loc[cells].values
#         Batch_D = np.array([self.drug_data[drug] for drug in drugs])
#         batch_y = self.ic50_data[self.col_ic50].iloc[indices]
#         batch_y = np.array(batch_y.astype("float32"))
#         return [*Batch_C, Batch_D], batch_y
#         
#     def on_epoch_end(self) :
#         self.indices = np.arange(len(self.ic50_data))
#         if self.shuffle :
#             np.random.shuffle(self.indices)


if __name__ == '__main__' :
    
    def parse_parameter() :
        import argparse
        parser = argparse.ArgumentParser()
        
        parser.add_argument("-nth", type=int, default=0)
        parser.add_argument("-choice", type=int, default=0)
        parser.add_argument("-ic50", type=str, default="IC50_GSDC2.txt")
        
        parser.add_argument("-col_cell", type=str, default="Cell_COSMIC")
        parser.add_argument("-col_drug", type=str, default="Drug")
        parser.add_argument("-col_ic50", type=str, default="LN_IC50")
        
        return parser.parse_args()
    
    # We confirmed that shuffle of train-loader works!!!
    from utils_def import *
    args = parse_parameter()
    ic50_data = pd.read_csv("_data/IC50_GDSC2.txt", header=0, sep="\t")

    cell_data = np.load("_data/cell_mut_matrix.npy", encoding="latin1", allow_pickle=True).item()
    drug_data = np.load("_data/drug_onehot_smiles.npy", encoding="latin1", allow_pickle=True).item()
    cell_data = {k:v for k,v in zip(cell_data["cell_names"], cell_data["cell_mut"])}
    drug_data = {k:v for k,v in zip(drug_data["drug_cids"], drug_data["canonical"])}
    
    ic50_data = ic50_data[:5000]
    ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
    ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=0, nth=0)
    
    batch_size = 100
    train = Batch(ic50_train, cell_data, drug_data, args, batch_size=batch_size, shuffle=True)
    valid = Batch(ic50_valid, cell_data, drug_data, args, batch_size=batch_size, shuffle=False)
    test = Batch(ic50_test, cell_data, drug_data, args, batch_size=batch_size, shuffle=False)
    
    indices_train = []
    indices_test = []
    print("Train Shuffle : {}".format(train.shuffle))
    print("Test Shuffle : {}".format(test.shuffle))
    
    for i in range(5) :
        train.reset()
        test.reset()
        print(train.indices[:10])
        print(test.indices[:10])
        
        while(train.available()):
            if train.offset==0 : indices_train.append(train.indices.copy())
            train_values, train_drugs, train_cells = train.mini_batch()
        
        while(test.available()) :
            if test.offset==0 : indices_test.append(test.indices.copy())
            test_values, test_drugs, test_cells = test.mini_batch()
        
        if i>0 :
            sim_train = np.mean(indices_train[0]==indices_train[i])
            sim_test = np.mean(indices_test[0]==indices_test[i])
            print("Index Similarity [{}th, Train] : {}".format(i, round(sim_train, 5)))
            print("Index Similarity [{}th, Test] : {}".format(i, round(sim_test, 5)))

#!/usr/bin/env python

import os
import re
import torch
import random
import numpy as np
import pandas as pd


def filt_ic50_(ic50_data, cell_data, drug_data, args) :

    ic50_data = ic50_data.astype({args.col_cell:"str"})
    ic50_data = ic50_data.astype({args.col_drug:"str"})
    # cell_list = list(cell_data.keys())
    # drug_list = list(drug_data.keys())
    cell_list = list(cell_data.index.astype(str))
    drug_list = list(drug_data.index.astype(str))

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
        if args.model_filepath=="" :
            dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        else :
            dir2 = args.model_filepath.split("/")[-1]
            dir2 = re.sub(".pt", "", dir2)
        dir_out = "Results/{}/{}".format(dir1, dir2)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    
    if args.model_filepath!="" :
        # Pre-trained model
        args.dir_test = "{}/pred_test.csv".format(dir_out)
        args.dir_param = args.model_filepath
        args.dir_hparam = "{}/hyper_param.json".format(dir_out)
    else :
        # Train & Test model from scratch
        args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
        
        if not retrain : 
            args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_{}.json".format(dir_out, args.nth)
        else :
            args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_retrain_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_retrain_{}.json".format(dir_out, args.nth)

    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def create_dir_out_lite(args, dir_out=None, suf_param="pt", retrain=False) :
    
    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        if args.model_filepath=="" :
            dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        else : 
            dir2 = args.model_filepath.split("/")[-1]
            dir2 = re.sub(".pt", "", dir2)
        dir3 = args.params_filepath.split("/")[-1]
        dir3 = re.sub(".json", "", dir3)
        dir_out = "Results/{}/{}/{}".format(dir1, dir2, dir3)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    
    if args.model_filepath!="" : 
        # Pre-trained model
        args.dir_test = "{}/pred_test.csv".format(dir_out)
        args.dir_param = args.model_filepath
        args.dir_hparam = "{}/hyper_param.json".format(dir_out)
    else :
        # Train & Test model from scratch
        args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
        
        if not retrain : 
            args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_{}.json".format(dir_out, args.nth)
        else :
            args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_retrain_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_retrain_{}.json".format(dir_out, args.nth)

    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def create_dir_out_seed(args, dir_out=None, suf_param="pt", retrain=False) :

    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        if args.model_filepath=="" :
            dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        else :
            dir2 = args.model_filepath.split("/")[-1]
            dir2 = re.sub(".pt", "", dir2)
        dir_out = "Results/{}/{}/Seed{}".format(dir1, dir2, args.seed)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    
    if args.model_filepath!="" :
        # Pre-trained model
        args.dir_test = "{}/pred_test.csv".format(dir_out)
        args.dir_param = args.model_filepath
        args.dir_hparam = "{}/hyper_param.json".format(dir_out)
    else :
        # Train & Test model from scratch
        args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
        
        if not retrain : 
            args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_{}.json".format(dir_out, args.nth)
        else :
            args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_retrain_{}.{}.tar".format(dir_out, args.nth, suf_param)
            args.dir_hparam = "{}/hyper_param_retrain_{}.json".format(dir_out, args.nth)

    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def ic50_split_path(args, ic50_train, ic50_valid, ic50_test, dir_ic50=None) :

    if dir_ic50 is None :
        dir_ic50 = "_data/Temp"
        os.makedirs(dir_ic50, exist_ok=True)

        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        
        if ic50_valid is not None :
            args.dir_train_temp = "{}/{}_{}_{}_Train.csv".format(dir_ic50, dir1, dir2, args.nth)
            args.dir_valid_temp = "{}/{}_{}_{}_Valid.csv".format(dir_ic50, dir1, dir2, args.nth)
            args.dir_test_temp = "{}/{}_{}_{}_Test.csv".format(dir_ic50, dir1, dir2, args.nth)
        else :
            args.dir_train_temp = "{}/{}_{}_{}_Train_Retrain.csv".format(dir_ic50, dir1, dir2, args.nth)
            args.dir_test_temp = "{}/{}_{}_{}_Test_Retrain.csv".format(dir_ic50, dir1, dir2, args.nth)
    
    else :
        os.makedirs(dir_ic50, exist_ok=True)
        if ic50_valid is not None :
            args.dir_train_temp = "{}/Train.csv".format(dir_ic50)
            args.dir_valid_temp = "{}/Valid.csv".format(dir_ic50)
            args.dir_test_temp = "{}/Test.csv".format(dir_ic50)
        else :
            args.dir_train_temp = "{}/Train_Retrain.csv".format(dir_ic50)
            args.dir_test_temp = "{}/Test_Retrain.csv".format(dir_ic50)
    
    ic50_train.to_csv(args.dir_train_temp)
    if ic50_test is not None :
        ic50_test.to_csv(args.dir_test_temp)
    if ic50_valid is not None :
        ic50_valid.to_csv(args.dir_valid_temp)

    return args


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


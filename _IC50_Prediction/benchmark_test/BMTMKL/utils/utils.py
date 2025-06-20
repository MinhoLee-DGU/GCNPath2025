#!/usr/bin/env python

import re
import os
import numpy as np
import pandas as pd


def create_dir_out(args, dir_out=None, suf_param="pt", retrain=False) :
    
    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        dir_out = "Results/{}/{}".format(dir1, dir2)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
    
    if not retrain : 
        args.dir_log = "{}/perf_log_{}.txt".format(dir_out, args.nth)
        args.dir_time = "{}/perf_time_{}.txt".format(dir_out, args.nth)
        args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
        args.dir_param = "{}/param_{}.{}".format(dir_out, args.nth, suf_param)
        args.dir_hparam = "{}/hyper_param_{}.pickle".format(dir_out, args.nth)
    else :
        args.dir_log = "{}/perf_log_retrain_{}.txt".format(dir_out, args.nth)
        args.dir_time = "{}/perf_time_retrain_{}.txt".format(dir_out, args.nth)
        args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
        args.dir_param = "{}/param_retrain_{}.{}".format(dir_out, args.nth, suf_param)
        args.dir_hparam = "{}/hyper_param_retrain_{}.pickle".format(dir_out, args.nth)
    
    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def filt_ic50(ic50_data, cell_data, drug_data, args) :
    
    ic50_data.dropna(subset=[args.col_cell, args.col_drug], inplace=True)
    if args.col_cell=="COSMIC_ID" :
        ic50_data = ic50_data.astype({args.col_cell:"int"})
    
    ic50_data = ic50_data.astype({args.col_cell:"str"})
    ic50_data = ic50_data.astype({args.col_drug:"str"})
    
    cell_list = list(str(_) for _ in cell_data.index)
    drug_list = list(str(_) for _ in drug_data.index)

    cell_list = list(set(ic50_data[args.col_cell]) & set(cell_list))
    drug_list = list(set(ic50_data[args.col_drug]) & set(drug_list))
    ic50_data = ic50_data[ic50_data[args.col_cell].isin(cell_list)]
    ic50_data = ic50_data[ic50_data[args.col_drug].isin(drug_list)]

    print("# [Cell] {}".format(len(cell_list)))
    print("# [Drug] {}".format(len(drug_list)))
    print("# [IC50] {}".format(ic50_data.shape[0]))

    return ic50_data


def split_even(ic50_data, args, nth=0, n_splits=10, n_repeats=1, seed=2021, choice=0) :
                 
    from sklearn.model_selection import RepeatedKFold
    split = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    idx = list(split.split(ic50_data, groups=(ic50_data[args.col_cell].astype(str) + '_' + ic50_data[args.col_drug].astype(str))))
    
    idx_tn = idx[nth][0]
    idx_tt = idx[nth][1]
    ic50_train = ic50_data.iloc[idx_tn]
    ic50_test = ic50_data.iloc[idx_tt]
    
    ic50_train.reset_index(drop=True, inplace=True)
    ic50_test.reset_index(drop=True, inplace=True)
    return ic50_train, ic50_test


def split_blind(ic50_data, args, nth=0, n_splits=10, n_repeats=1, seed=2021, choice=1) :
  
    from sklearn.model_selection import RepeatedKFold
    split = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    cells = ic50_data[args.col_cell].unique()
    drugs = ic50_data[args.col_drug].unique()
    idx_cell = list(split.split(cells))
    idx_drug = list(split.split(drugs))
    
    if choice == 1 :
        # Cell Blind Test
        cell_tn = cells[idx_cell[nth][0]]
        cell_tt = cells[idx_cell[nth][1]]
        ic50_train = ic50_data[ic50_data[args.col_cell].isin(cell_tn)]
        ic50_test = ic50_data[ic50_data[args.col_cell].isin(cell_tt)]
    elif choice == 2 :
        # Drug Blind Test
        drug_tn = drugs[idx_drug[nth][0]]
        drug_tt = drugs[idx_drug[nth][1]]
        ic50_train = ic50_data[ic50_data[args.col_drug].isin(drug_tn)]
        ic50_test = ic50_data[ic50_data[args.col_drug].isin(drug_tt)]
    else :
        # Cell & Drug Blind Test [Strict]
        # All combinations of cell x drug index
        idx_test = np.empty(shape=(0, 2), dtype=int)
        
        for i in range(n_repeats) :
            idx_test_temp = np.arange(n_splits * i, n_splits * (i + 1))
            idx_test_temp = np.array(np.meshgrid(idx_test_temp, idx_test_temp)).T.reshape(-1, 2)
            idx_test = np.concatenate([idx_test, idx_test_temp], axis=0)
            
        nth_cell = idx_test[nth][0]
        nth_drug = idx_test[nth][1]
        
        cell_tn = cells[idx_cell[nth_cell][0]]
        cell_tt = cells[idx_cell[nth_cell][1]]
        drug_tn = drugs[idx_drug[nth_drug][0]]
        drug_tt = drugs[idx_drug[nth_drug][1]]
        
        ic50_train = ic50_data[(ic50_data[args.col_cell].isin(cell_tn) & ic50_data[args.col_drug].isin(drug_tn))]
        ic50_test = ic50_data[(ic50_data[args.col_cell].isin(cell_tt) & ic50_data[args.col_drug].isin(drug_tt))]

    ic50_train.reset_index(drop=True, inplace=True)
    ic50_test.reset_index(drop=True, inplace=True)
    return ic50_train, ic50_test


def split_ic50(ic50_data, args, choice=0, nth=0, n_splits=10, 
               n_splits_strict=5, n_repeats=3, seed=2021, retrain=False) :
    
    # Use this function after the function filt_ic50
    split_def = split_even if choice==0 else split_blind
    n_splits = n_splits if choice in [0, 1, 2] else n_splits_strict
    ic50_train, ic50_test = split_def(ic50_data, args, nth=nth, choice=choice,
                                      n_splits=n_splits, n_repeats=n_repeats, seed=seed)
    
    if choice not in [0, 1, 2, 3] :
        raise Exception("\nchoice of IC50 split : 0 [Normal], 1 [Cell_Blind], 2 [Drug_Blind], 3 [Strict_Blind]")
    
    if not retrain :
        ic50_train, ic50_valid = split_def(ic50_train, args, nth=0, choice=choice, 
                                           n_splits=n_splits-1, n_repeats=n_repeats, seed=seed)
        nc_valid = ic50_valid[args.col_cell].nunique()
        nd_valid = ic50_valid[args.col_drug].nunique()

    nc_train = ic50_train[args.col_cell].nunique()
    nd_train = ic50_train[args.col_drug].nunique()
    nc_test = ic50_test[args.col_cell].nunique()
    nd_test = ic50_test[args.col_drug].nunique()
    
    if retrain :
        print("# [Cell] Train : Test = {} : {}".format(nc_train, nc_test))
        print("# [Drug] Train : Test = {} : {}".format(nd_train, nd_test))
        print("# [IC50] Train : Test = {} : {}".format(len(ic50_train), len(ic50_test)))
        return ic50_train, ic50_test
    else :
        print("# [Cell] Train : Valid : Test = {} : {} : {}".format(nc_train, nc_valid, nc_test))
        print("# [Drug] Train : Valid : Test = {} : {} : {}".format(nd_train, nd_valid, nd_test))
        print("# [IC50] Train : Valid : Test = {} : {} : {}".format(len(ic50_train), len(ic50_valid), len(ic50_test)))
        return ic50_train, ic50_valid, ic50_test

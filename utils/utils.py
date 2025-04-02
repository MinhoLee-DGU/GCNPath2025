#!/usr/bin/env python

import re
import os
import numpy as np
import pandas as pd


def calc_rmse(y_test, y_predict, digit=3) :
    from sklearn.metrics import mean_squared_error
    value = mean_squared_error(y_test, y_predict, squared=False)
    return round(float(value), digit)


def calc_r2(y_test, y_predict, digit=3) :
    from sklearn.metrics import r2_score
    value = r2_score(y_test, y_predict)
    return round(float(value), digit)


def calc_pcc(y_test, y_predict, digit=3) :
    from scipy.stats.stats import pearsonr
    value = pearsonr(y_test, y_predict)[0]
    return round(float(value), digit)


def calc_scc(y_test, y_predict, digit=3) :
    from scipy.stats.stats import spearmanr
    value = spearmanr(y_test, y_predict)[0]
    return round(float(value), digit)


def norm_ic50(ln_ic50, reverse=False) :
    if not reverse :
        ln_ic50 = np.exp(ln_ic50)
        ln_ic50 = 1 / (1 + ln_ic50**(-0.1))
    else :
        ln_ic50 = (1/ln_ic50 - 1)**(-10)
        ln_ic50 = np.log(ln_ic50)
    return ln_ic50


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
    if args.col_cell in ["COSMIC_ID"] :
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


def pred_to_csv(pred_y, ic50_test, dir_pred, grad_cam=False) :
    if grad_cam :
        ic50_pred = pred_y
    else :
        ic50_pred = pd.DataFrame(pred_y, columns=['Prediction'])
    
    ic50_pred = ic50_pred.reset_index(drop=True)
    ic50_test = ic50_test.reset_index(drop=True)
    ic50_pred = pd.concat([ic50_test, ic50_pred], axis=1)
    ic50_pred.to_csv(dir_pred, header=True, index=False)


def pred_to_csv_norm(pred_y, ic50_test, dir_pred) :
    ic50_pred = pd.DataFrame(pred_y, columns=["Prediction_Norm"])
    ic50_pred["Prediction"] = norm_ic50(ic50_pred["Prediction_Norm"], reverse=True)
    ic50_pred = ic50_pred.reset_index(drop=True)
    ic50_test = ic50_test.reset_index(drop=True)
    ic50_pred = pd.concat([ic50_test, ic50_pred], axis=1)
    ic50_pred.to_csv(dir_pred, header=True, index=False)


def time_to_csv(train_time_list=None, test_time=None, valid_time=None, dir_time=None) :
    if train_time_list is None :
        # Test only
        steps = ["Test"]
        times = [test_time]
    elif valid_time is None :
        # Retrain total dataset without valid dataset
        steps = ["Train_Sum", "Train_Mean"]
        times = [np.sum(train_time_list), np.mean(train_time_list)]
        steps_epoch = [f"Train_Epoch{epoch+1}" for epoch in range(len(train_time_list))]
        steps.extend(steps_epoch)
        times.extend(train_time_list)
    else :
        # Train with train, valid, test dataset
        steps = ["Train_Sum", "Train_Mean", "Valid", "Test"]
        times = [np.sum(train_time_list), np.mean(train_time_list), valid_time, test_time]
        steps_epoch = [f"Train_Epoch{epoch+1}" for epoch in range(len(train_time_list))]
        steps.extend(steps_epoch)
        times.extend(train_time_list)
    
    times_info = pd.DataFrame({"Step": steps, "Time": times})
    times_info.to_csv(dir_time, header=True, index=False)
    print(times_info)


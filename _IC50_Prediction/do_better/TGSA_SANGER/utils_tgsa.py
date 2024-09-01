#!/usr/bin/env python

import os
import re
import torch
import numpy as np
import pandas as pd


def create_dir_out_(args, dir_out=None, suf_param=None, tgsa=False, retrain=False) :
    
    if dir_out is None :
        dir1 = args.ic50.split("/")[-1]
        dir1 = re.sub(".txt|.csv|.xlsx", "", dir1)
        dir2 = ["Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind"][args.choice]
        dir3 = "TGDRP" if not tgsa else "TGSA"
    
    if args.pretrain==1 :
        dir2 = "SA" if tgsa else "TGDRP"
        args.weight_tgdrp = "./weights/TGDRP.pth"
        args.weight_tgsa = "./weights/SA.pth"
    elif args.pretrain==2 :
        dir2 = "SA_pre" if tgsa else "TGDRP_pre"
        args.weight_tgdrp = "./weights/TGDRP_pre.pth"
        args.weight_tgsa = "./weights/SA_pre.pth"
    else :
        if not retrain: 
            args.weight_tgdrp = "Results/{}/{}/TGDRP/param_{}.{}".format(dir1, dir2, args.nth, suf_param)
            args.weight_tgsa = "Results/{}/{}/TGSA/param_{}.{}".format(dir1, dir2, args.nth, suf_param)
        else :
            args.weight_tgdrp = "Results/{}/{}/TGDRP/param_retrain_{}.{}".format(dir1, dir2, args.nth, suf_param)
            args.weight_tgsa = "Results/{}/{}/TGSA/param_retrain_{}.{}".format(dir1, dir2, args.nth, suf_param)
    
    if args.pretrain==0 :
        dir_out = "Results/{}/{}/{}".format(dir1, dir2, dir3)
        args.dir_valid = "{}/pred_valid_{}.csv".format(dir_out, args.nth)
        
        if not retrain : 
            args.dir_test = "{}/pred_test_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_{}.{}".format(dir_out, args.nth, suf_param)
        else :
            args.dir_test = "{}/pred_test_retrain_{}.csv".format(dir_out, args.nth)
            args.dir_param = "{}/param_retrain_{}.{}".format(dir_out, args.nth, suf_param)
        
    else :
        dir_out = "Results/{}/{}".format(dir1, dir2)
        args.dir_valid = "{}/pred_valid.csv".format(dir_out)
        args.dir_test = "{}/pred_test.csv".format(dir_out)
    
    args.dir_out = dir_out
    os.makedirs(dir_out, exist_ok=True)
    args.dir_hparam = "{}/hyper_param_{}.pickle".format(dir_out, args.nth)
    
    print("\n### Save all files into {}...".format(dir_out))
    print("# Save parameters into {}...".format(args.dir_param))
    return args


def summary_model(model, verbose=True) :
    if verbose : 
        print("\n### Model Summary\n{}\n".format(model))

    param_total = sum(p.numel() for p in model.parameters())
    param_cell = sum(p.numel() for p in model.GNN_cell.parameters())
    param_drug = sum(p.numel() for p in model.GNN_drug.parameters())
    
    param_cell_ = sum(p.numel() for p in model.cell_emb.parameters())
    param_drug_ = sum(p.numel() for p in model.drug_emb.parameters())
    param_final = sum(p.numel() for p in model.regression.parameters())
    
    print("# Model Parameters : {}".format(param_total))
    print("# Model Parameters : {} [Cell]".format(param_cell))
    print("# Model Parameters : {} [Drug]".format(param_drug))
    print("# Model Parameters : {} [Cell_]".format(param_cell_))
    print("# Model Parameters : {} [Drug_]".format(param_drug_))
    print("# Model Parameters : {} [Final]".format(param_final))


def summary_model_sa(model, verbose=True) :
    if verbose :
        print("\n### Model Summary\n{}\n".format(model))

    param_total = sum(p.numel() for p in model.parameters())
    print("# Model Parameters [SA] : {}".format(param_total))

#!/usr/bin/env python

import os
import sys
import random
import logging
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

import torch
from torch.nn import MSELoss
from torch.utils.data import Dataset, DataLoader
from torch.nn import Linear, Tanh, ReLU, ELU, LeakyReLU, Softmax
from torch.optim.lr_scheduler import ExponentialLR, CosineAnnealingLR

from utils.utils import *
from torch_geometric.data import Batch


class DatasetDef(Dataset) :
    def __init__(self, ic50_data, cell_data, drug_data, args, no_labels=False) :
        super(DatasetDef, self).__init__()
        
        self.ic50_data = ic50_data
        self.cell_data = cell_data
        self.drug_data = drug_data
        self.no_labels = no_labels
        
        self.cell_name = ic50_data[args.col_cell]
        self.drug_name = ic50_data[args.col_drug]
        
        if not self.no_labels :
            self.ic50_value = ic50_data[args.col_ic50]

    def __len__(self):
        return len(self.ic50_data)

    def __getitem__(self, idx) :
        cell = self.cell_name[idx]
        drug = self.drug_name[idx]
        
        cell_data_ = self.cell_data[cell]
        drug_data_ = self.drug_data[drug]
        
        if not self.no_labels :
            ic50_data_ = self.ic50_value[idx]
            return cell_data_, drug_data_, ic50_data_
        else :
            return cell_data_, drug_data_


class CollateDef(object) :
    def __init__(self, gnn_cell=True, gnn_drug=True, no_labels=False, sa_ddi=False) :
        super().__init__()
        
        self.sa_ddi = sa_ddi
        self.gnn_cell = gnn_cell
        self.gnn_drug = gnn_drug
        self.no_labels = no_labels
    
    def collate(self, data_list, gnn=True, sa_ddi=False) :
        if not gnn :
            batch = torch.stack([data for data in data_list], dim=0)
        elif sa_ddi :
            batch = Batch.from_data_list(data_list, follow_batch=['edge_index'])
        else :
            batch = Batch.from_data_list(data_list)
        
        return batch
    
    def __call__(self, batch) :
        if self.no_labels :
            cell_data_, drug_data_ = map(list, zip(*batch))
        else :
            cell_data_, drug_data_, ic50_data_ = map(list, zip(*batch))
            ic50_data_batch = torch.tensor(ic50_data_)
        
        cell_data_batch = self.collate(cell_data_, gnn=self.gnn_cell, sa_ddi=False)
        drug_data_batch = self.collate(drug_data_, gnn=self.gnn_drug, sa_ddi=self.sa_ddi)
        
        if self.no_labels :
            return cell_data_batch, drug_data_batch
        else :
            return cell_data_batch, drug_data_batch, ic50_data_batch


class Logger_Perf :
    def __init__(self, file) :
        self.logger = logging.getLogger("Logger")
        self.logger.setLevel(logging.DEBUG)

        # stream_handler = logging.StreamHandler()
        # self.logger.addHandler(stream_handler)
            
        file_handler = logging.FileHandler(file)
        self.logger.addHandler(file_handler)
    
    def set_column(self, no_valid=False) :
        col_info = "Epoch\tEarly_Stop"
        metrics = ["MSE", "RMSE", "PCC", "SCC", "Time"]
        metrics_train = "\t".join(["Train_{}".format(m) for m in metrics])
        metrics_valid = "\t".join(["Valid_{}".format(m) for m in metrics])
        
        if no_valid :
            self.logger.info("{}\t{}".format(col_info, metrics_train))
        else :
            self.logger.info("{}\t{}\t{}".format(col_info, metrics_train, metrics_valid))
    
    def info(self, epoch, early_stop, perf_train, perf_valid=None) :
        perf_train = "\t".join([str(round(m, 3)) for m in perf_train])
        
        if perf_valid is None :
            self.logger.info("{}\t{}\t{}".format(epoch, False, perf_train))
        else :
            perf_valid = "\t".join([str(round(m, 3)) for m in perf_valid])
            self.logger.info("{}\t{}\t{}\t{}".format(epoch, early_stop, perf_train, perf_valid))


def set_random_seed(seed, deterministic=True) :
    
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    print("# Set random seed : {}".format(seed))
    
    if deterministic:
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        # torch.use_deterministic_algorithms(True)


def seed_worker(worker_id) :
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def load_data(ic50_data, cell_data, drug_data, args, no_labels=False,
              batch_size=256, num_workers=0, shuffle=False, fix_seed=False) :
    
    sa_ddi = args.attn_mode==4
    gnn_cell = args.gnn_cell!="linear"
    gnn_drug = args.gnn_drug!="linear"
    collate_fn = CollateDef(gnn_cell=gnn_cell, gnn_drug=gnn_drug, no_labels=no_labels, sa_ddi=sa_ddi)
    
    ic50_data.reset_index(drop=True, inplace=True)
    dataset = DatasetDef(ic50_data, cell_data, drug_data, args, no_labels=no_labels)
    
    _init_fn = seed_worker if fix_seed else None
    dataloader = DataLoader(dataset, collate_fn=collate_fn, batch_size=batch_size, 
                            num_workers=num_workers, shuffle=shuffle, worker_init_fn=_init_fn)
    
    return dataloader


def choose_params(args, verbose=False) :
  
    gnn_cell_list = ["linear", "gcn", "gat", "gin", "rgcn", "rgat"]
    gnn_drug_list = ["linear", "gcn", "gat", "gin", "gine"]
    mode_list = ["plain", "res", "res+", "dense", "jk_cat", "jk_max"]
    act_list = [ReLU(), ELU(0.1), ELU(0.5), ELU(1.0), LeakyReLU(0.1), LeakyReLU(0.2)]
    
    args.act = act_list[args.act]
    args.gnn_cell = gnn_cell_list[args.gnn_cell]
    args.gnn_drug = gnn_drug_list[args.gnn_drug]
    args.mode_cell = mode_list[args.mode_cell]
    args.mode_drug = mode_list[args.mode_drug]
    
    if args.gnn_cell=="linear" or args.gnn_drug=="linear" : args.attn_mode = False
    if verbose : print("# Get parameters in args : act, gnn_cell, gnn_drug, mode_cell, mode_drug")
    return args


def get_cell_info(args, sample_cell, dim_list=None, verbose=False) :
    
    if args.gnn_cell!="linear" :
        n_path = sample_cell.x.shape[0]
        n_rel = hasattr(sample_cell, "edge_type")
        n_rel = len(sample_cell.edge_type.unique()) if n_rel else None
        n_cell_feat = sample_cell.x.shape[1]
    else :
        n_path = None
        n_rel = None
        n_cell_feat = sample_cell.shape[0]
    
    dim = args.dim_cell
    n_hidden = args.n_hid_cell
    dim_def = [n_cell_feat] + [dim for _ in range(n_hidden)]
    dim_cell = dim_def if dim_list is None else dim_list
    # GNN [1, 8, 8, 8] or Linear [n_path, 256, 256, 256]
    
    args.n_rel = n_rel
    args.n_path = n_path
    args.dim_cell = dim_cell
    if verbose : print("# Get parameters : n_rel, n_path, dim_cell")
    
    return args


def get_drug_info(args, sample_drug, dim_list=None, verbose=False) :
    
    if args.gnn_drug!="linear" :
        n_drug_edge = hasattr(sample_drug, "edge_attr")
        n_drug_edge = sample_drug.edge_attr.shape[1] if n_drug_edge else None
        n_drug_feat = sample_drug.x.shape[1]
    
    else :
        n_drug_edge = None
        n_drug_feat = sample_drug.shape[0]
    
    dim = args.dim_drug
    n_hidden = args.n_hid_drug
    dim_def = [n_drug_feat] + [dim for _ in range(n_hidden)]
    dim_drug = dim_def if dim_list is None else dim_list
    # GNN [n_fp, 256, 256, 256] or Linear [192, 256, 256, 256]
    
    args.n_drug_edge = n_drug_edge
    args.dim_drug = dim_drug
    if verbose : print("# Get parameters : n_drug_edge, dim_drug")
    
    return args


def get_pred_info(args, dim_list=None, verbose=False) :
    
    dim = args.dim_pred
    n_hidden = args.n_hid_pred
    dim_def = [dim for _ in range(n_hidden)]
    dim_pred = dim_def if dim_list is None else dim_list
    
    args.dim_pred = dim_pred
    if verbose : print("# Get parameters : dim_pred")
    
    return args


def choose_scheduler(scheduler, gamma=0.99, T_max=10, lr=0.001) :
    if scheduler==0 :
        scheduler = None
    elif scheduler==1 :
        scheduler = ExponentialLR(optimizer, gamma=gamma)
    elif scheduler==2 :
        scheduler = CosineAnnealingLR(optimizer, T_max=T_max, eta_min=lr/2)
    else : raise Exception("\nchoice of scheduler : 0 [None], 1 [Exponential], 2 [Cosine]")
    
    return scheduler


def train_epoch(model, train_loader, device, optimizer) :
    
    model.train()
    model.to(device)
    loss_epoch = 0.0
    loss_fn = MSELoss()
    
    y_pred = torch.Tensor()
    y_true = torch.Tensor()
    
    if device!="cpu" :
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
    
    for batch, (cell, drug, ic50) in enumerate(tqdm(train_loader, desc="Train")) :
        cell = cell.to(device)
        drug = drug.to(device)
        ic50 = ic50.view(-1, 1).float().to(device)
        
        pred = model(cell, drug, return_attn=False)
        loss_batch = loss_fn(pred, ic50)
        loss_epoch += loss_batch.item()
        
        optimizer.zero_grad()
        loss_batch.backward()
        optimizer.step()
        
        y_pred = torch.cat((y_pred, pred.cpu()))
        y_true = torch.cat((y_true, ic50.cpu()))
    
    if device!="cpu" :
        end.record()
        torch.cuda.synchronize()
        time = start.elapsed_time(end)
    else :
        time = 0
    
    loss_epoch = loss_epoch / len(train_loader)
    y_pred = y_pred.detach().squeeze().numpy()
    y_true = y_true.detach().squeeze().numpy()
    
    rmse = calc_rmse(y_pred, y_true, digit=3)
    pcc = calc_pcc(y_pred, y_true, digit=3)
    scc = calc_scc(y_pred, y_true, digit=3)
    
    return loss_epoch, rmse, pcc, scc, time


def valid_epoch(model, valid_loader, device) :
    
    model.eval()
    model.to(device)
    loss_epoch = 0.0
    loss_fn = MSELoss()
    
    y_pred = torch.Tensor()
    y_true = torch.Tensor()
    
    if device!="cpu" :
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
    
    with torch.no_grad() :
        for batch, (cell, drug, ic50) in enumerate(tqdm(valid_loader, desc="Valid")) :
            cell = cell.to(device)
            drug = drug.to(device)
            ic50 = ic50.view(-1, 1).float().to(device)
            
            pred = model(cell, drug, return_attn=False)
            loss_batch = loss_fn(pred, ic50)
            loss_epoch += loss_batch.item()
            
            y_pred = torch.cat((y_pred, pred.cpu()))
            y_true = torch.cat((y_true, ic50.cpu()))
    
    if device!="cpu" :
        end.record()
        torch.cuda.synchronize()
        time = start.elapsed_time(end)
    else :
        time = 0
        
    loss_epoch = loss_epoch / len(valid_loader)
    y_pred = y_pred.detach().squeeze().numpy()
    y_true = y_true.detach().squeeze().numpy()
    
    rmse = calc_rmse(y_pred, y_true, digit=3)
    pcc = calc_pcc(y_pred, y_true, digit=3)
    scc = calc_scc(y_pred, y_true, digit=3)
    
    return loss_epoch, rmse, pcc, scc, time


def train(model, train_loader, valid_loader, device, optimizer,
          epochs=300, patience=10, dir_param=None, scheduler=None, dir_log=None) :
    
    if dir_log is not None : 
        Logger = Logger_Perf(dir_log)
        Logger.set_column()
    
    time_list = []
    early_stopping = EarlyStopping(patience=patience, path=dir_param, verbose=False)
    
    for epoch in range(1, epochs+1) :
        train_loss, train_rmse, train_pcc, train_scc, train_time = train_epoch(model, train_loader, device, optimizer)
        valid_loss, valid_rmse, valid_pcc, valid_scc, valid_time = valid_epoch(model, valid_loader, device)
        if scheduler is not None : scheduler.step()
        time_list.append(train_time+valid_time)
        
        epoch_str = str(epoch).zfill(3)
        early_stopping(valid_loss, model, optimizer, scheduler)
        print("# [Epoch {}] Train & Valid Loss : {} & {}".format(epoch_str, round(train_loss, 3), round(valid_loss, 3)))
        
        if dir_log is not None : 
            train_perf = [train_loss, train_rmse, train_pcc, train_scc, train_time]
            valid_perf = [valid_loss, valid_rmse, valid_pcc, valid_scc, valid_time]
            Logger.info(epoch, early_stopping.early_stop, train_perf, valid_perf)
            
        if early_stopping.early_stop :
            print("### Early Stopping Occurred...")
            break
    
    # Load the last checkpoint with the best model
    checkpoint = torch.load(dir_param, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    return model, time_list


def train_wo_valid(model, train_loader, device, optimizer,
                   epochs=100, dir_param=None, scheduler=None, dir_log=None) :
    
    if dir_log is not None : 
        Logger = Logger_Perf(dir_log)
        Logger.set_column(no_valid=True)
    
    time_list = []
    for epoch in range(1, epochs+1) :
        train_loss, train_rmse, train_pcc, train_scc, train_time = train_epoch(model, train_loader, device, optimizer)
        if scheduler is not None : scheduler.step()
        time_list.append(train_time)
        
        epoch_str = str(epoch).zfill(3)
        print("# [Epoch {}] Train Loss : {}".format(epoch_str, round(train_loss, 3)))
        
        if dir_log is not None : 
            train_perf = [train_loss, train_rmse, train_pcc, train_scc, train_time]
            Logger.info(epoch, False, train_perf)
        
        scheduler_ = None if scheduler is None else scheduler.state_dict()
        torch.save({"model_state_dict" : model.state_dict(), 
                    "optimizer_state_dict" : optimizer.state_dict(),
                    "scheduler_state_dict" : scheduler_}, dir_param)
    
    # Load the last checkpoint
    checkpoint = torch.load(dir_param, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    return model, time_list


def test(model, test_loader, device, return_attn=False) :
    
    model.eval()
    model.to(device)
    loss_epoch = 0.0
    loss_fn = MSELoss()
    
    y_pred = torch.Tensor()
    y_true = torch.Tensor()
    
    attn = torch.Tensor()
    attn_path = torch.Tensor()
    attn_subs = torch.Tensor()
    
    if device!="cpu" :
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
    
    with torch.no_grad() :
        for batch, (cell, drug, ic50) in enumerate(tqdm(test_loader, desc="Test")) :
            cell = cell.to(device)
            drug = drug.to(device)
            ic50 = ic50.view(-1, 1).float().to(device)
            
            if model.attn_mode in [1, 2] and return_attn :
                pred, attn_ = model(cell, drug, return_attn=return_attn)
                attn = torch.cat((attn, attn_.cpu()))
            elif model.attn_mode in [3, 4] and return_attn :
                pred, attn_path_, attn_subs_ = model(cell, drug, return_attn=return_attn)
                attn_path = torch.cat((attn_path, attn_path_.cpu()))
                attn_subs = torch.cat((attn_subs, attn_subs_.cpu()))
            else :
                pred = model(cell, drug, return_attn=return_attn)
            
            y_pred = torch.cat((y_pred, pred.cpu()))
            y_true = torch.cat((y_true, ic50.cpu()))
            loss_batch = loss_fn(pred, ic50)
            loss_epoch += loss_batch.item()
    
    if device!="cpu" :
        end.record()
        torch.cuda.synchronize()
        time = start.elapsed_time(end)
    else :
        time = 0
    
    loss_epoch = loss_epoch / len(test_loader)
    y_pred = y_pred.detach().squeeze().numpy()
    y_true = y_true.detach().squeeze().numpy()
    
    rmse = calc_rmse(y_pred, y_true, digit=3)
    pcc = calc_pcc(y_pred, y_true, digit=3)
    scc = calc_scc(y_pred, y_true, digit=3)
    
    print("# Performance (RMSE) : {}".format(rmse))
    print("# Performance (PCC) : {}".format(pcc))
    print("# Performance (SCC) : {}".format(scc))
    print("# Performance (SCC) : {}".format(time))
    
    if model.attn_mode in [1, 2] and return_attn :
        return y_pred, attn, time
    elif model.attn_mode in [3, 4] and return_attn :
        return y_pred, attn_path, attn_subs, time
    else :
        return y_pred, time


def test_no_labels(model, test_loader, device, return_attn=False) :
    
    model.eval()
    model.to(device)
    pred_test = torch.Tensor()
    
    attn = torch.Tensor()
    attn_path = torch.Tensor()
    attn_subs = torch.Tensor()
    
    if device!="cpu" :
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
    
    with torch.no_grad() :
        for batch, (cell, drug) in enumerate(tqdm(test_loader, desc="Test")) :
            cell = cell.to(device)
            drug = drug.to(device)
            
            if model.attn_mode in [1, 2] and return_attn :
                pred, attn_ = model(cell, drug, return_attn=return_attn)
                attn = torch.cat((attn, attn_.cpu()))
            elif model.attn_mode in [3, 4] and return_attn :
                pred, attn_path_, attn_subs_ = model(cell, drug, return_attn=return_attn)
                attn_path = torch.cat((attn_path, attn_path_.cpu()))
                attn_subs = torch.cat((attn_subs, attn_subs_.cpu()))
            else :
                pred = model(cell, drug, return_attn=return_attn)
            
            pred_test = torch.cat((pred_test, pred.cpu()))
    
    if device!="cpu" :
        end.record()
        torch.cuda.synchronize()
        time = start.elapsed_time(end)
    else :
        time = 0
    
    pred_test = pred_test.squeeze().numpy()
    
    if model.attn_mode in [1, 2] and return_attn :
        return pred_test, attn, time
    elif model.attn_mode in [3, 4] and return_attn :
        return pred_test, attn_path, attn_subs, time
    else :
        return pred_test, time


def grad_cam(model, cell_data, drug_data, ic50_data, args):
    model.to(args.device)
    model.eval()
    importance = []
    
    for i in range(ic50_data.shape[0]) :
        cell = ic50_data[args.col_cell][i]
        drug = ic50_data[args.col_drug][i]
        cell = Batch.from_data_list([cell_data[cell]]).to(args.device)
        drug = Batch.from_data_list([drug_data[drug]]).to(args.device)
        
        pred, cell_explain = model(cell, drug, grad_cam=True)
        pred.backward()
        
        imp_weight = torch.mean(cell_explain.grad, dim=0)
        imp = ReLU()((cell_explain*imp_weight).sum(dim=1))
        imp = imp / imp.sum()
        importance.append(imp.detach().cpu())
    
    importance = torch.stack(importance).squeeze().numpy()
    return importance


class EarlyStopping :

    def __init__(self, patience=10, verbose=False, delta=0, path='checkpoint.pt') :
        self.path = path
        self.delta = delta
        self.verbose = verbose
        self.patience = patience
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf

    def __call__(self, val_loss, model, optimizer, scheduler=None) :
        score = -val_loss
        if self.best_score is None :
            self.best_score = score
            self.save_checkpoint(val_loss, model, optimizer, scheduler)
        else :
            if score < self.best_score + self.delta :
                self.counter += 1
                print(f"EarlyStopping Counter: {self.counter} out of {self.patience}")
                if self.counter >= self.patience :
                    self.early_stop = True
            else:
                self.best_score = score
                self.save_checkpoint(val_loss, model, optimizer, scheduler)
                self.counter = 0

    def save_checkpoint(self, val_loss, model, optimizer, scheduler=None) :
        if self.verbose :
            print(f"Valid Loss Decreased ({self.val_loss_min:.6f} --> {val_loss:.6f})")
        scheduler_ = None if scheduler is None else scheduler.state_dict()
        
        self.val_loss_min = val_loss
        torch.save({"valid_loss_min" : self.val_loss_min,
                    "model_state_dict" : model.state_dict(), 
                    "optimizer_state_dict" : optimizer.state_dict(),
                    "scheduler_state_dict" : scheduler_}, self.path)
        

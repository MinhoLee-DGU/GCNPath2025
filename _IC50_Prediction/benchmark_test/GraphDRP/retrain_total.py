import numpy as np
import pandas as pd
import sys, os
from random import shuffle
import pdb
import torch
import torch.nn as nn
from models.gat import GATNet
from models.gat_gcn import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet
import datetime
import argparse
import pdb
import pickle
from torch_geometric.data import DataLoader
import random

from utils import *
from utils_def import *
from utils_graphdrp import *


# training function at each epoch
def train(model, device, train_loader, optimizer, epoch, log_interval):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    model.train()
    loss_fn = nn.MSELoss()
    avg_loss = []
    for batch_idx, data in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        output, _ = model(data)
        loss = loss_fn(output, data.y.view(-1, 1).float().to(device))
        loss.backward()
        optimizer.step()
        avg_loss.append(loss.item())
        if batch_idx % log_interval == 0:
            print('Train epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,
                                                                           batch_idx * len(data.x),
                                                                           len(train_loader.dataset),
                                                                           100. * batch_idx / len(train_loader),
                                                                           loss.item()))
    return sum(avg_loss)/len(avg_loss)


def predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for idx, data in enumerate(loader):
            data = data.to(device)
            output, _ = model(data)
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu()), 0)
    return total_labels.numpy().flatten(), total_preds.numpy().flatten()


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


# def main(modeling, train_batch, val_batch, test_batch, lr, num_epoch, log_interval, cuda_name):
def main(modeling, train_batch, lr, num_epoch, log_interval, args) :

    print('Learning rate: ', lr)
    print('Epochs: ', num_epoch)
    set_random_seed(args.seed)

    model_st = modeling.__name__
    dataset = 'GDSC'
    train_losses = []
    val_losses = []
    val_pearsons = []
    # print('\nrunning on ', model_st + '_' + dataset )
    
    # processed_data_file_train = 'data/processed/' + dataset + '_train_mix'+'.pt'
    # processed_data_file_val = 'data/processed/' + dataset + '_val_mix'+'.pt'
    # processed_data_file_test = 'data/processed/' + dataset + '_test_mix'+'.pt'

    # train_data = TestbedDataset(root='data', dataset=dataset+'_train_mix')
    # val_data = TestbedDataset(root='data', dataset=dataset+'_val_mix')
    # test_data = TestbedDataset(root='data', dataset=dataset+'_test_mix')
    
    # make data PyTorch mini-batch processing ready
    nth = args.nth
    ic50 = args.ic50
    choice = args.choice
    args = create_dir_out_(args, suf_param="pt", retrain=True)
    
    args.dir_time = "{}/log_time_seed_{}.csv".format(args.dir_out, args.seed)
    args.dir_pred = "{}/pred_total_seed{}.csv".format(args.dir_out, args.seed)
    args.dir_param = "{}/param_retrain_seed{}.pt".format(args.dir_out, args.seed)
    args.dir_test = "{}/pred_chembl_seed{}.csv".format(args.dir_out, args.seed)
    
    if "ChEMBL" in args.ic50:
        args.dir_time = "{}/log_time_chembl_seed_{}.csv".format(args.dir_out, args.seed)
    
    cell_data = pd.read_csv(args.dir_cell, index_col=0)
    cell_data.index = cell_data.index.astype(str)
    cell_data = {k:v for k, v in zip(cell_data.index, cell_data.values)}
    
    with open(args.dir_drug, "rb") as f :
        drug_data = pickle.load(f)
    
    ic50_data = pd.read_csv(ic50, sep="\t")
    ic50_data["Norm_IC50"] = norm_ic50(ic50_data[args.col_ic50])

    args.col_ic50 = "Norm_IC50"
    ic50_data = filt_ic50_(ic50_data, cell_data, drug_data, args)
    
    if args.mode=="train" :
        train_data = DatasetDef(ic50_data, cell_data, drug_data, args)
        # ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, retrain=True)
        # ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
        # train_data = DatasetDef(ic50_train, cell_data, drug_data, args)
        # valid_data = DatasetDef(ic50_valid, cell_data, drug_data, args)
        # test_data = DatasetDef(ic50_test, cell_data, drug_data, args)
        
        train_loader = DataLoader(train_data, batch_size=train_batch, shuffle=True, num_workers=args.cpu)
        train_loader_ = DataLoader(train_data, batch_size=train_batch, shuffle=False, num_workers=args.cpu)
        # val_loader = DataLoader(valid_data, batch_size=val_batch, shuffle=False, num_workers=args.cpu)
        # test_loader = DataLoader(test_data, batch_size=test_batch, shuffle=False, num_workers=args.cpu)
    else :
        test_data = DatasetDef(ic50_data, cell_data, drug_data, args)
        test_loader = DataLoader(test_data, batch_size=test_batch, shuffle=False, num_workers=args.cpu)
    
    print("CPU/GPU: ", torch.cuda.is_available())
    print("Save the prediction file in {}".format(args.dir_test))
    
    # training the model
    # device = torch.device(cuda_name if torch.cuda.is_available() else "cpu")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = modeling().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    
    if args.mode=="test" :
        if args.pretrain==1 :
            checkpoint = torch.load("model_GINConvNet_GDSC.model", map_location=device)
        else :
            checkpoint = torch.load(args.dir_param, map_location=device)
        model.load_state_dict(checkpoint)
    
    if args.mode=="train" :
        best_mse = 1000
        best_pearson = 1
        best_epoch = -1
        
        early_stop = 0
        early_stop_count = 10
        model_file_name = args.dir_param
        result_file_name = args.dir_perf
        
        # model_file_name = 'model_' + model_st + '_' + dataset +  '.model'
        # result_file_name = 'result_' + model_st + '_' + dataset +  '.csv'
        # loss_fig_name = 'model_' + model_st + '_' + dataset + '_loss'
        # pearson_fig_name = 'model_' + model_st + '_' + dataset + '_pearson'
        
        train_time_list = []
        args.device = device
        
        if args.device!="cpu" :
            start = torch.cuda.Event(enable_timing=True)
            end = torch.cuda.Event(enable_timing=True)
            start.record()
        
        for epoch in range(num_epoch):
            train_loss = train(model, device, train_loader, optimizer, epoch+1, log_interval)
            
            if args.device!="cpu" :
                end.record()
                torch.cuda.synchronize()
                train_time = start.elapsed_time(end)
            else :
                train_time = 0
            
            train_time_list.append(train_time)
            torch.save(model.state_dict(), model_file_name)
            
            # G,P = predicting(model, device, val_loader)
            # ret = [rmse(G,P),mse(G,P),pearson(G,P),spearman(G,P)]
            # G_test, P_test = predicting(model, device, test_loader)
            # ret_test = [rmse(G_test, P_test), mse(G_test, P_test), pearson(G_test, P_test), spearman(G_test, P_test)]
    
            # train_losses.append(train_loss)
            # val_losses.append(ret[1])
            # val_pearsons.append(ret[2])
    
            # if ret[1]<best_mse:
            # if ret_test[1]<best_mse:
            #     torch.save(model.state_dict(), model_file_name)
            #     with open(result_file_name,'w') as f:
            #         f.write(','.join(map(str,ret_test)))
            #     
            #     early_stop = 0
            #     # best_mse = ret[1]
            #     # best_pearson = ret[2]
            #     best_mse = ret_test[1]
            #     best_pearson = ret_test[2]
            #     best_epoch = epoch + 1
            #     print(' rmse improved at epoch ', best_epoch, '; best_mse:', best_mse, model_st, dataset)
            # else:
            #     early_stop += 1
            #     print(' no improvement since epoch ', best_epoch, '; best_mse, best pearson:', best_mse, best_pearson, model_st, dataset)
            #     
            #     if early_stop >= early_stop_count :
            #         print("Early stopping occured...")
            #         break
        
        if args.device!="cpu" :
            start = torch.cuda.Event(enable_timing=True)
            end = torch.cuda.Event(enable_timing=True)
            start.record()
        
        G_Total, P_Total = predicting(model, device, train_loader_)
        
        if args.device!="cpu" :
            end.record()
            torch.cuda.synchronize()
            test_time = start.elapsed_time(end)
        else :
            test_time = 0
        
        pred_to_csv_norm(P_Total, ic50_data, args.dir_pred)
        time_to_csv(train_time_list, test_time, dir_time=args.dir_time)
        
    else :
        args.device = device
        if args.device!="cpu" :
            start = torch.cuda.Event(enable_timing=True)
            end = torch.cuda.Event(enable_timing=True)
            start.record()
        else :
            import time
            start = time.perf_counter()
      
        G_test, P_test = predicting(model, device, test_loader)
        
        if args.device!="cpu" :
            end.record()
            torch.cuda.synchronize()
            test_time = start.elapsed_time(end)
        else :
            end = time.perf_counter()
            test_time = 1000*(end - start)
        
        pred_to_csv_norm(P_test, ic50_data, args.dir_test)
        time_to_csv(test_time=test_time, dir_time=args.dir_time)
        # draw_loss(train_losses, val_losses, loss_fig_name)
        # draw_pearson(val_pearsons, pearson_fig_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='train model')
    parser.add_argument('--model', type=int, required=False, default=0, help='0: GINConvNet, 1: GATNet, 2: GAT_GCN, 3: GCNNet')
    parser.add_argument('--train_batch', type=int, required=False, default=1024,  help='Batch size training set')
    parser.add_argument('--val_batch', type=int, required=False, default=1024, help='Batch size validation set')
    parser.add_argument('--test_batch', type=int, required=False, default=1024, help='Batch size test set')
    parser.add_argument('--lr', type=float, required=False, default=1e-4, help='Learning rate')
    parser.add_argument('--num_epoch', type=int, required=False, default=100, help='Number of epoch')
    parser.add_argument('--log_interval', type=int, required=False, default=20, help='Log interval')
    # parser.add_argument('--cuda_name', type=str, required=False, default="cuda:0", help='Cuda')
    
    parser.add_argument("-nth", type=int, default=0)
    parser.add_argument("-cpu", type=int, default=4)
    parser.add_argument("-seed", type=int, default=2021)
    parser.add_argument("-choice", type=int, default=0)
    parser.add_argument("-ic50", type=str, default="IC50_GDSC2.txt")
    
    parser.add_argument("-col_cell", type=str, default="Cell_COSMIC")
    parser.add_argument("-col_drug", type=str, default="Drug")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50")
    
    parser.add_argument("-mode", type=str, default="train")
    parser.add_argument("-pretrain", type=int, default=0)
    parser.add_argument("-test_name", type=str, default="chembl")
    
    parser.add_argument("-dir_cell", type=str, default="cell_feat.csv")
    parser.add_argument("-dir_drug", type=str, default="drug_dict")
    args = parser.parse_args()

    modeling = [GINConvNet, GATNet, GAT_GCN, GCNNet][args.model]
    train_batch = args.train_batch
    val_batch = args.val_batch
    test_batch = args.test_batch
    lr = args.lr
    num_epoch = args.num_epoch
    log_interval = args.log_interval
    # cuda_name = args.cuda_name
    
    main(modeling, train_batch, lr, num_epoch, log_interval, args)
    # main(modeling, train_batch, val_batch, test_batch, lr, num_epoch, log_interval, cuda_name)
    

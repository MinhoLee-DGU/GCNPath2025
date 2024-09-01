import os
import torch
import torch.nn as nn
import numpy as np
import pandas as pd

# from utils import load_data
from torch.utils.data import Dataset, DataLoader
from utils import EarlyStopping, set_random_seed
from utils import train, validate
from preprocess_gene import get_STRING_graph, get_predefine_cluster
from models.TGDRP import TGDRP

import argparse
import fitlog
from utils import *
from utils_def import *
from utils_tgsa import *


def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42, help='seed')
    parser.add_argument('--model', type=str, default='TGDRP', help='Name of the model')
    parser.add_argument('--batch_size', type=int, default=128, help='batch size (default: 128)')
    parser.add_argument('--lr', type=float, default=0.0001, help='learning rate')
    parser.add_argument('--layer_drug', type=int, default=3, help='layer for drug')
    parser.add_argument('--dim_drug', type=int, default=128, help='hidden dim for drug')
    parser.add_argument('--layer', type=int, default=3, help='number of GNN layer')
    parser.add_argument('--hidden_dim', type=int, default=8, help='hidden dim for cell')
    parser.add_argument('--weight_decay', type=float, default=0, help='weight decay')
    parser.add_argument('--dropout_ratio', type=float, default=0.2, help='dropout ratio')
    parser.add_argument('--epochs', type=int, default=300, help='maximum number of epochs (default: 300)')
    parser.add_argument('--patience', type=int, default=10, help='patience for earlystopping (default: 10)')
    parser.add_argument('--edge', type=float, default=0.95, help='threshold for cell line graph')
    parser.add_argument('--pretrain', type=int, default=0, help='whether use pre-trained weights (0 for False, 1 for True')
    parser.add_argument('--weight_path', type=str, default='', help='filepath for pretrained weights')
    parser.add_argument('--mode', type=str, default='test', help='train or test')
    
    parser.add_argument("-nth", type=int, default=0)
    parser.add_argument("-cpu", type=int, default=4)
    parser.add_argument("-choice", type=int, default=0) 
    
    parser.add_argument("-ic50", type=str, default="_data/IC50_GDSC2.txt")
    parser.add_argument("-cell", type=str, default="_data/cell_feature_all.npy")
    parser.add_argument("-drug", type=str, default="_data/drug_feature_graph.npy")
    
    parser.add_argument("-col_cell", type=str, default="Cell")
    parser.add_argument("-col_drug", type=str, default="Drug")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50")
    
    parser.add_argument("-dir_test", type=str, default=None)
    parser.add_argument("-dir_param", type=str, default=None)
    
    return parser.parse_args()


def main():
    args = arg_parse()
    set_random_seed(args.seed)
    
    nth = args.nth
    choice = args.choice
    dir_ic50 = "./_data"
    input_ic50 = os.path.join(dir_ic50, args.ic50)
    args.device = "cuda" if torch.cuda.is_available() else "cpu"
    
    if args.dir_test is not None and args.dir_param is not None :
        args.weight_tgdrp = args.dir_param
    else : 
        args = create_dir_out_(args, suf_param="pth", tgsa=False, retrain=True)
    
    ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")
    cell_data = np.load(args.cell, allow_pickle=True).item()
    drug_data = np.load(args.drug, allow_pickle=True).item()
    
    edge_index = np.load("./_data/edge_index_PPI_{}.npy".format(args.edge))
    cluster_predefine = np.load("./_data/cluster_predefine_PPI_{}.npy".format(args.edge), allow_pickle=True).item()
    cluster_predefine = {i: j.to(args.device) for i, j in cluster_predefine.items()}
    
    # genes_path = './data_custom/cell_data'
    # edge_index = get_STRING_graph(genes_path, args.edge)
    # cluster_predefine = get_predefine_cluster(edge_index, genes_path, args.edge, args.device)
    
    cell_data = {str(k) : cell_data[k] for k in cell_data.keys()}
    drug_data = {str(k) : drug_data[k] for k in drug_data.keys()}
    ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
    
    if args.mode=="train" :
        ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, retrain=True)
        # ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
        train_data = MyDataset(drug_data, cell_data, ic50_train, edge_index, args)
        # valid_data = MyDataset(drug_data, cell_data, ic50_valid, edge_index, args)
        test_data = MyDataset(drug_data, cell_data, ic50_test, edge_index, args)
        
        train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True, collate_fn=collate_fn, num_workers=args.cpu)
        # val_loader = DataLoader(valid_data, batch_size=args.batch_size, shuffle=False, collate_fn=collate_fn, num_workers=args.cpu)
        test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, collate_fn=collate_fn, num_workers=args.cpu)
    
    else :
        ic50_test = ic50_data
        test_data = MyDataset(drug_data, cell_data, ic50_test, edge_index, args)
        test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, collate_fn=collate_fn, num_workers=args.cpu)
    
    sample_cell = list(cell_data.keys())[0]
    args.num_feature = cell_data[sample_cell].x.shape[1]
    model = TGDRP(cluster_predefine, args).to(args.device)
    print('mean degree:{}'.format(len(edge_index[0]) / 706))
    summary_model(model)
  
    if args.mode == 'train':
        # if args.pretrain and args.weight_path != '':
        #     model.GNN_drug.load_state_dict(torch.load('./{}.pth'.format(args.weight_path), map_location=args.device)['model_state_dict'])
        
        criterion = nn.MSELoss()
        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        
        # log_folder = dir_hparam
        # log_folder = os.path.join(os.getcwd(), "logs", model._get_name())
        
        # if not os.path.exists(log_folder):
        #     os.makedirs(log_folder)
        # fitlog.set_log_dir(log_folder)
        # fitlog.add_hyper(args)
        # fitlog.add_hyper_in_file(__file__)

        stopper = EarlyStopping(mode='lower', patience=args.patience, filename=args.dir_param)
        for epoch in range(1, args.epochs + 1):
            print("=====Epoch {}".format(epoch))
            print("Training...")
            train_loss = train(model, train_loader, criterion, opt, args.device)
            # fitlog.add_loss(train_loss.item(), name='Train MSE', step=epoch)

            print('Evaluating...')
            # rmse, _, _, _, _ = validate(model, val_loader, args.device)
            # print("Validation rmse:{}".format(rmse))
            rmse, _, _, _, _ = validate(model, test_loader, args.device)
            print("Test rmse:{}".format(rmse))
            # fitlog.add_metric({'val': {'RMSE': rmse}}, step=epoch)

            early_stop = stopper.step(rmse, model)
            if early_stop:
                break

        print('EarlyStopping! Finish training!')
        print('Testing...')
        stopper.load_checkpoint(model)

        train_rmse, train_MAE, train_r2, train_r, _ = validate(model, train_loader, args.device)
        # val_rmse, val_MAE, val_r2, val_r, pred_valid = validate(model, val_loader, args.device)
        test_rmse, test_MAE, test_r2, test_r, pred_test = validate(model, test_loader, args.device)
        print('Train reslut: rmse:{} r2:{} r:{}'.format(train_rmse, train_r2, train_r))
        # print('Val reslut: rmse:{} r2:{} r:{}'.format(val_rmse, val_r2, val_r))
        print('Test reslut: rmse:{} r2:{} r:{}'.format(test_rmse, test_r2, test_r))
        
        # pred_to_csv(pred_valid, ic50_valid, dir_pred=args.dir_valid)
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)
        
        # Save hyper-parameters
        with open(args.dir_hparam, "wb") as f :
            pickle.dump(args, f)
        print("\n### Model hyper-parameters saved!!!")
        
        # fitlog.add_best_metric(
        #     {'epoch': epoch - args.patience,
        #      "train": {'RMSE': train_rmse, 'MAE': train_MAE, 'pearson': train_r, "R2": train_r2},
        #      "valid": {'RMSE': stopper.best_score, 'MAE': val_MAE, 'pearson': val_r, 'R2': val_r2},
        #      "test": {'RMSE': test_rmse, 'MAE': test_MAE, 'pearson': test_r, 'R2': test_r2}})

    else :
        # weight = "TGDRP_pre" if args.pretrain else "TGDRP"
        checkpoint = torch.load(args.weight_tgdrp, map_location=args.device)['model_state_dict']
        
        model.load_state_dict(checkpoint)
        print("Success!!")
        test_rmse, test_MAE, test_r2, test_r, pred_test = validate(model, test_loader, args.device)
        print('Test RMSE: {}, MAE: {}, R2: {}, R: {}'.format(round(test_rmse.item(), 4), round(test_MAE, 4), round(test_r2, 4), round(test_r, 4)))
        
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)
        
if __name__ == "__main__":
    main()


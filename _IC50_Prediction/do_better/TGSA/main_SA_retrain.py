import os
import torch
import torch.nn as nn
import numpy as np
import pandas as pd

from models.SA import SA
import argparse
import fitlog
import pickle

import pandas as pd
from utils import *
from utils_def import *
from utils_tgsa import *


def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42, help='random seed (default: 42)')
    parser.add_argument('--knn', type=int, default=5, help='k-nearest-neighbour')
    parser.add_argument('--layer_drug', type=int, default=3, help='layer for drug')
    parser.add_argument('--dim_drug', type=int, default=128, help='hidden dim for drug')
    parser.add_argument('--layer', type=int, default=3, help='number of GNN layer')
    parser.add_argument('--hidden_dim', type=int, default=8, help='hidden dim for cell')
    parser.add_argument('--batch_size', type=int, default=512, help='batch size (default: 512)')
    parser.add_argument('--lr', type=float, default=0.0001, help='learning rate (default: 0.0001)')
    parser.add_argument('--weight_decay', type=float, default=0, help='weight decay')
    parser.add_argument('--dropout_ratio', type=float, default=0.2, help='dropout ratio')
    parser.add_argument('--epochs', type=int, default=300, help='maximum number of epochs (default: 100)')
    parser.add_argument('--patience', type=int, default=10, help='patience for earlystopping (default: 7)')
    parser.add_argument('--edge', type=float, default='0.95', help='edge for gene graph')
    parser.add_argument('--mode', type=str, default='test', help='train or test')
    parser.add_argument('--pretrain', type=int, default=0, help='pretrain(0 or 1)')
    parser.add_argument('--weight_path', type=str, default='', help='filepath for pretrained weights')
    
    parser.add_argument("-nth", type=int, default=0)
    parser.add_argument("-cpu", type=int, default=4)
    parser.add_argument("-choice", type=int, default=0)
    
    parser.add_argument("-ic50", type=str, default="_data/IC50_GDSC2.txt")
    parser.add_argument('-cell', type=str, default="_data/cell_feature_all.npy", help='Cell')
    parser.add_argument('-drug', type=str, default="_data/drug_feature_graph.npy", help='Drug')
    
    parser.add_argument("-col_cell", type=str, default="Cell_BROAD")
    parser.add_argument("-col_drug", type=str, default="Drug")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50")
    
    parser.add_argument("-dir_test", type=str, default=None)
    parser.add_argument("-dir_param", type=str, default=None)
    parser.add_argument("-dir_tgsa", type=str, default=None)
    
    parser.add_argument("-dir_sim", type=str, default="_data/drug_cell_edges_5_knn_gdsc")
    parser.add_argument("-dir_cidx", type=str, default="_data/cell_id2idx_dict_gdsc")
    # parser.add_argument("-dir_didx", type=str, default="_data/drug_name2idx_dict_gdsc")
    
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
        # Test mode
        args.weight_tgdrp = args.dir_param
        args.weight_tgsa = args.dir_tgsa
    else :
        # Train mode, Pre-train mode
        args = create_dir_out_(args, suf_param="pth", tgsa=True, retrain=True)
    
    ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")
    cell_data = np.load(args.cell, allow_pickle=True).item()
    drug_data = np.load(args.drug, allow_pickle=True).item()
    
    edge_index = np.load("./_data/edge_index_PPI_{}.npy".format(args.edge))
    cluster_predefine = np.load("./_data/cluster_predefine_PPI_{}.npy".format(args.edge), allow_pickle=True).item()
    cluster_predefine = {i: j.to(args.device) for i, j in cluster_predefine.items()}
    
    if args.ic50=="IC50_GDSC.txt" :
        suffix = "_gdsc"
        # Cell 700, Drug 432, IC50 269519
    elif args.ic50=="IC50_GDSC1.txt" :
        suffix = "_gdsc1"
        # Cell 696, Drug 321, IC50 191356
    elif args.ic50=="IC50_GDSC2.txt" :
        suffix = "_gdsc2" 
        # Cell 697, Drug 236, IC50 144744
    else :
        suffix = "_gdsc"
    
    dict_dir = './_data/'
    dir_sim = "./_data/drug_cell_edges_5_knn" + suffix
    dir_cidx_dict = dict_dir+"cell_id2idx_dict"+suffix
    # dir_didx_dict = dict_dir+"drug_name2idx_dict"+suffix
    
    if args.mode=="test" and args.pretrain==0 :
        # Test model without pretrained weight... [GDSC1+2]
        dir_sim = args.dir_sim
        dir_cidx_dict = args.dir_cidx
        # dir_didx_dict = args.dir_didx
    
    print("Similarity KNN : {}".format(dir_sim))
    print("Cell idx dict : {}".format(dir_cidx_dict))
    # print("Drug idx dict : {}".format(dir_didx_dict))
    
    with open(dir_cidx_dict, "rb") as f :
        cell_name2idx = pickle.load(f)
    # with open(dir_didx_dict, "rb") as f :
    #     drug_name2idx = pickle.load(f)
    
    drug_data = {str(k):v for k, v in drug_data.items()}
    cell_data = {str(k):v for k, v in cell_data.items()}
    
    # TGSA github code does not utilize drug similarity graph...
    # No need to prepare drug similarity graph and the corresponding drug index...
    drug_name2idx = {k:v for k, v in zip(drug_data.keys(), range(len(drug_data)))}
    drug_data = {k:drug_data[k] for k in drug_name2idx.keys()}
    cell_data = {k:cell_data[k] for k in cell_name2idx.keys()}
    ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
    
    if args.mode=="train" :
        ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, retrain=True)
        # ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
        train_data = MyDataset_SA(drug_name2idx, cell_name2idx, ic50_train, args)
        # valid_data = MyDataset_SA(drug_name2idx, cell_name2idx, ic50_valid, args)
        test_data = MyDataset_SA(drug_name2idx, cell_name2idx, ic50_test, args)
        
        train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True, num_workers=args.cpu)
        # val_loader = DataLoader(valid_data, batch_size=args.batch_size, shuffle=False, num_workers=args.cpu)
        test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, num_workers=args.cpu)
    
    else :
        ic50_test = ic50_data
        test_data = MyDataset_SA(drug_name2idx, cell_name2idx, ic50_test, args)
        test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, num_workers=args.cpu)
    
    # train_loader, val_loader, test_loader = load_data_SA(args)
    drug_nodes_data, cell_nodes_data, drug_edges, cell_edges, parameter = load_graph_data_SA(args, drug_data, cell_data, edge_index, cluster_predefine, dir_sim)
    model = SA(drug_nodes_data, cell_nodes_data, drug_edges, cell_edges, args).to(args.device)
    
    print("TGDRP Parameters : {}".format(args.weight_tgdrp))
    print("TGSA Parameters : {}".format(args.weight_tgsa))
    
    if args.mode == "train":
        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        criterion = nn.MSELoss()
        # weight = "TGDRP_pre" if args.pretrain else "TGDRP"
        model.regression = parameter['regression']
        model.drug_emb = parameter['drug_emb']
        summary_model_sa(model)
        
        # log_folder = os.path.join(os.getcwd(), "logs", model._get_name())
        # if not os.path.exists(log_folder):
        #     os.makedirs(log_folder)
        # fitlog.set_log_dir(log_folder)
        # fitlog.add_hyper(args)
        # fitlog.add_hyper_in_file(__file__)
        
        set_random_seed(args.seed)
        stopper = EarlyStopping(mode='lower', patience=args.patience, filename=args.dir_param)
        print("train_num: {},  test_num: {}".format(len(train_loader.dataset), len(test_loader.dataset)))
        # print("train_num: {},  val_num: {}, test_num: {}".format(len(train_loader.dataset), 
        #                                                          len(val_loader.dataset),
        #                                                          len(test_loader.dataset)))
        
        for epoch in range(1, args.epochs + 1):
            print("=====Epoch {}".format(epoch))
            print("Training...")
            train_SA(train_loader, model, criterion, opt, args)
            print('Evaluating...')
            rmse, _, _, _, _ = validate_SA(test_loader, model, args)
            print("Test rmse:{}".format(rmse))
            # rmse, _, _, _, _ = validate_SA(val_loader, model, args)
            # print("Validation rmse:{}".format(rmse))
            early_stop = stopper.step(rmse, model)
            if early_stop: break

        print('EarlyStopping! Finish training!')
        print('Testing...')
        stopper.load_checkpoint(model)

        train_rmse, train_MAE, train_r2, train_r, pred_train = validate_SA(train_loader, model, args)
        # val_rmse, val_MAE, val_r2, val_r, pred_valid = validate_SA(val_loader, model, args)
        test_rmse, test_MAE, test_r2, test_r, pred_test = validate_SA(test_loader, model, args)
        print('Train reslut: rmse:{} r2:{} r:{}'.format(train_rmse, train_r2, train_r))
        # print('Val reslut: rmse:{} r2:{} r:{}'.format(val_rmse, val_r2, val_r))
        print('Test reslut: rmse:{} r2:{} r:{}'.format(test_rmse, test_r2, test_r))

        # fitlog.add_best_metric(
        #     {'epoch': epoch - args.patience,
        #     "train": {'RMSE': train_rmse, 'MAE': train_MAE, 'pearson': train_r, "R2": train_r2},
        #     "valid": {'RMSE': stopper.best_score, 'MAE': val_MAE, 'pearson': val_r, 'R2': val_r2},
        #     "test": {'RMSE': test_rmse, 'MAE': test_MAE, 'pearson': test_r, 'R2': test_r2}})
        # fitlog.finish()
        
        # pred_to_csv(pred_valid, ic50_valid, dir_pred=args.dir_valid)
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)
        
        # Save hyper-parameters
        with open(args.dir_hparam, "wb") as f :
            pickle.dump(args, f)
        print("\n### Model hyper-parameters saved!!!")

    elif args.mode == "test":
        # weight = "SA_pre" if args.pretrain else "SA"
        # model.load_state_dict(torch.load('./weights/{}.pth'.format(weight), map_location=args.device)['model_state_dict'])
        
        weight_sa = torch.load(args.weight_tgsa, map_location=args.device)['model_state_dict']
        model.load_state_dict(weight_sa)
        
        test_rmse, test_MAE, test_r2, test_r, pred_test = validate_SA(test_loader, model, args)
        print('Test RMSE: {}, MAE: {}, R2: {}, R: {}'.format(round(test_rmse.item(), 4), round(test_MAE, 4),
                                                             round(test_r2, 4), round(test_r, 4)))
        
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)

if __name__ == '__main__':
    main()

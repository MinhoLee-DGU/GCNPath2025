import os
import re
import sys
import pickle
import argparse

import torch
import torch.nn as nn
import numpy as np
import pandas as pd

from utils import *
from utils_def import *
from utils_drpreter import *

from Model.DRPreter import DRPreter
from Model.Similarity import Similarity
from torch_scatter import scatter_add


def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42, help='seed')
    parser.add_argument('--batch_size', type=int, default=128, help='batch size (default: 128)')
    parser.add_argument('--lr', type=float, default=0.0001, help='learning rate (default: 0.0001)')
    parser.add_argument('--layer', type=int, default=3, help='Number of cell layers')
    parser.add_argument('--hidden_dim', type=int, default=8, help='hidden dim for cell')
    parser.add_argument('--layer_drug', type=int, default=3, help='Number of drug layers')
    parser.add_argument('--dim_drug', type=int, default=128, help='hidden dim for drug (default: 128)')
    parser.add_argument('--dim_drug_cell', type=int, default=256, help='hidden dim for drug and cell (default: 256)')
    parser.add_argument('--dropout_ratio', type=float, default=0.1, help='Dropout ratio (default: 0.1)')
    parser.add_argument('--epochs', type=int, default=300, help='Maximum number of epochs (default: 300)')
    parser.add_argument('--patience', type=int, default=10, help='patience for early stopping (default: 10)')
    parser.add_argument('--mode', type=str, default='train', help='train, test')
    parser.add_argument('--edge', type=str, default='STRING', help='STRING, BIOGRID') # BIOGRID: removed
    parser.add_argument('--string_edge', type=float, default=0.99, help='Threshold for edges of cell line graph')
    parser.add_argument('--dataset', type=str, default='disjoint', help='2369joint, 2369disjoint, COSMIC')
    parser.add_argument('--trans', type=bool, default=True, help='Use Transformer or not')
    parser.add_argument('--sim', type=int, default=0, help='Construct homogeneous similarity networks or not')
    
    parser.add_argument("-nth", type=int, default=0)
    parser.add_argument("-cpu", type=int, default=4)
    parser.add_argument("-choice", type=int, default=0)
    
    parser.add_argument("-ic50", type=str, default="IC50_GDSC2.txt")
    parser.add_argument("-cell", type=str, default="_data/cell_feature_std_disjoint.npy")
    parser.add_argument("-drug", type=str, default="_data/drug_feature_graph.npy")
    
    parser.add_argument("-col_cell", type=str, default="Cell")
    parser.add_argument("-col_drug", type=str, default="Drug")
    parser.add_argument("-col_ic50", type=str, default="LN_IC50")
    
    parser.add_argument("-pretrain", type=int, default=0)
    parser.add_argument("-dir_test", type=str, default=None)
    parser.add_argument("-dir_param", type=str, default=None)
    
    return parser.parse_args()


def main():
    args = arg_parse()
    rpath = './'
    result_path = rpath + 'Result/'
    set_random_seed(args.seed)
    
    nth = args.nth
    choice = args.choice
    dir_ic50 = "_data"
    input_ic50 = os.path.join(dir_ic50, args.ic50)
    args.device = "cuda" if torch.cuda.is_available() else "cpu"
    
    if args.dir_test is not None and args.dir_param is not None :
        pass
    else : 
        args = create_dir_out_(args, suf_param="pth", retrain=True)
    
    ic50_data = pd.read_csv(input_ic50, header=0, sep="\t")
    edge_index = np.load(rpath+f'Data/Cell/edge_index_PPI_990_2369{args.dataset}.npy')
    drug_data = np.load(args.drug, allow_pickle=True).item() # pyg format of drug graph
    cell_data = np.load(args.cell, allow_pickle=True).item() # pyg data format of cell graph

    drug_data = {str(k):v for k, v in drug_data.items()}
    cell_data = {str(k):v for k, v in cell_data.items()}
    ic50_data = filt_ic50(ic50_data, cell_data, drug_data, args)
    
    if args.sim==1 :
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
            suffix = None
        
        dict_dir = './_data/'
        dir_sim = "./_data/drug_cell_edges_5_knn" + suffix
        
        dir_cidx_dict = dict_dir+"cell_id2idx_dict"+suffix
        dir_didx_dict = dict_dir+"drug_name2idx_dict"+suffix
        print("Cell idx dict : {}".format(dir_cidx_dict))
        print("Drug idx dict : {}".format(dir_didx_dict))
        
        with open(dir_cidx_dict, "rb") as f :
            cell_name2idx = pickle.load(f)
        with open(dir_didx_dict, "rb") as f :
            drug_name2idx = pickle.load(f)
            
        drug_data = {k:drug_data[k] for k in drug_name2idx.keys()}
        cell_data = {k:cell_data[k] for k in cell_name2idx.keys()}
    
    example = list(cell_data.values())[0]
    args.num_feature = example.x.shape[1] # 1
    args.num_genes = example.x.shape[0] # 4646
    print(f'num_feature: {args.num_feature}, num_genes: {args.num_genes}')
    
    
    if 'disjoint' in args.dataset:
        gene_list = scatter_add(torch.ones_like(example.x.squeeze()), example.x_mask.to(torch.int64)).to(torch.int)
        args.max_gene = gene_list.max().item()
        args.cum_num_nodes = torch.cat([gene_list.new_zeros(1), gene_list.cumsum(dim=0)], dim=0)
        args.n_pathways = gene_list.size(0)
        print('num_genes:{}, num_edges:{}'.format(args.num_genes, len(edge_index[0])))
        print('gene distribution: {}'.format(gene_list))
        print('mean degree:{}'.format(len(edge_index[0]) / args.num_genes))
    else:
        print('num_genes:{}, num_edges:{}'.format(args.num_genes, len(edge_index[0])))
        print('mean degree:{}'.format(len(edge_index[0]) / args.num_genes))
    
    
    if args.mode=="train" and not args.sim :
        ic50_train, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth, retrain=True)
        train_loader, test_loader = load_data_retrain(ic50_train, ic50_test, drug_data, cell_data, edge_index, args)
    elif args.mode=="test" and not args.sim :
        test_loader = load_data_test(ic50_data, drug_data, cell_data, edge_index, args)
    elif args.mode=="train" and args.sim :
        # ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, args, choice=choice, nth=nth)
        # train_loader, val_loader, test_loader = load_sim_data(ic50_train, ic50_valid, ic50_test, drug_name2idx, cell_name2idx, args)
        # drug_nodes_data, cell_nodes_data, drug_edges, cell_edges = load_sim_graph(drug_data, cell_data, edge_index, args.dir_drpreter, dir_sim, args)
        raise Exception("DRPreter SA does not provide inductive inferences...")
    else :
        # test_loader = load_sim_data_test(ic50_test, drug_name2idx, cell_name2idx, args)
        raise Exception("DRPreter SA does not provide inductive inferences...")
    
    if not args.sim :
        model = DRPreter(args).to(args.device)
        summary_model(model)
    else :
        model = Similarity(drug_nodes_data, cell_nodes_data, drug_edges, cell_edges, args).to(args.device)
        summary_model_sa(model)
    
    # # ---- [1] Pathway + Transformer ----
    # if not args.sim :
    #      ic50_train, ic50_valid, ic50_test = split_ic50(ic50_data, cell_data, drug_data, args, choice=choice, nth=nth)
    #     train_loader, val_loader, test_loader = load_data(ic50_train, ic50_valid, ic50_test, drug_data, cell_data, edge_index, args)
    #     print('train: {}, val: {}, test: {}'.format(len(train_loader.dataset), len(val_loader.dataset), len(test_loader.dataset)))
    #     model = DRPreter(args).to(args.device)
    #     summary_model(model)
    #     # print(model)
    #     
    # # ---- [2] Add similarity information after obtaining embeddings ----
    # else:
    #     # dict_dir = "./Data/similarity/dict/"
    #     knn_dir = "./Data/Similarity/edge/drug_cell_edges_5_knn"
    #     weight_dir = "./Results/{}/{}/DRPreter/param_{}.pth".format(dir1, dir2, nth)
    #     train_loader, val_loader, test_loader = load_sim_data(ic50_train, ic50_valid, ic50_test, drug_data, cell_data, args)
    #     print('train: {}, val: {}, test: {}'.format(len(train_loader.dataset), len(val_loader.dataset), len(test_loader.dataset)))
    #     
    #     drug_nodes_data, cell_nodes_data, drug_edges, cell_edges = load_sim_graph(drug_data, cell_data, edge_index, weight_dir, knn_dir, args)
    #     model = Similarity(drug_nodes_data, cell_nodes_data, drug_edges, cell_edges, args).to(args.device)
    #     summary_model_sa(model)
    #     # print(model)

        
# -----------------------------------------------------------------
        
    if args.mode == 'train':
        result_col = ('mse\trmse\tmae\tpcc\tscc')
        # result_type = 'results_sim' if args.sim==True else 'results'
        # results_path = get_path(args, result_path, result_type=result_type)
        results_path = args.dir_out + "/results_retrain_{}.txt".format(nth)
        
        with open(results_path, 'w') as f:
            f.write(result_col + '\n')
        criterion = nn.MSELoss()
        opt = torch.optim.Adam(model.parameters(), lr=args.lr)

        # state_dict_name = f'{rpath}weights/weight_sim_seed{args.seed}.pth' if args.sim==True else f'{rpath}weights/weight_seed{args.seed}.pth'
        stopper = EarlyStopping(mode='lower', patience=args.patience, filename=args.dir_param)
                
        for epoch in range(1, args.epochs + 1):
            print(f"===== Epoch {epoch} =====")
            train_loss = train(model, train_loader, criterion, opt, args)

            # mse, rmse, mae, pcc, scc, _ = validate(model, val_loader, args)
            mse, rmse, mae, pcc, scc, _ = validate(model, test_loader, args)
            results = [epoch, mse, rmse, mae, pcc, scc]
            save_results(results, results_path)
            
            # print(f"Validation mse: {mse}")
            # test_MSE, test_RMSE, test_MAE, test_PCC, test_SCC, _ = validate(model, test_loader, args)
            # print(f"Test mse: {test_MSE}")
            print(f"Test mse: {mse}")
            early_stop = stopper.step(mse, model)
            if early_stop:
                break

        print('EarlyStopping! Finish training!')
        print('Best epoch: {}'.format(epoch-stopper.counter))

        stopper.load_checkpoint(model)

        train_MSE, train_RMSE, train_MAE, train_PCC, train_SCC, _ = validate(model, train_loader, args)
        # val_MSE, val_RMSE, val_MAE, val_PCC, val_SCC, pred_valid = validate(model, val_loader, args)
        test_MSE, test_RMSE, test_MAE, test_PCC, test_SCC, pred_test = validate(model, test_loader, args)

        print('-------- DRPreter -------')
        print(f'sim: {args.sim}')
        print(f'##### Seed: {args.seed} #####')
        print('\t\tMSE\tRMSE\tMAE\tPCC\tSCC')
        print('Train result: {}\t{}\t{}\t{}\t{}'.format(r4(train_MSE), r4(train_RMSE), r4(train_MAE), r4(train_PCC), r4(train_SCC)))
        # print('Val result: {}\t{}\t{}\t{}\t{}'.format(r4(val_MSE), r4(val_RMSE), r4(val_MAE), r4(val_PCC), r4(val_SCC)))
        print('Test result: {}\t{}\t{}\t{}\t{}'.format(r4(test_MSE), r4(test_RMSE), r4(test_MAE), r4(test_PCC), r4(test_SCC)))
        
        # pred_to_csv(pred_valid, ic50_valid, dir_pred=args.dir_valid)
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)
        # df.to_csv(get_path(args, result_path, result_type=result_type+'_df', extension='csv'), sep='\t', index=0)

        
    else :
        # state_dict_name = f'{rpath}weights/weight_sim_seed{args.seed}.pth' if args.sim==True else f'{rpath}weights/weight_seed{args.seed}.pth'
        if args.pretrain==1 : args.dir_param = f'./weights/weight_seed{args.seed}.pth'
        model.load_state_dict(torch.load(args.dir_param, map_location=args.device)['model_state_dict'])

#         '''Get embeddings of specific drug and cell line pair'''
#         drug_name, cell_name = 'Bortezomib', 'ACH-000137' # 8MGBA
#         drug_emb, cell_emb = embedding(model, drug_name, cell_name, drug_data, cell_data, edge_index, args)
#         print(drug_emb, cell_emb)
        
        
        ''' Test results only '''
        test_MSE, test_RMSE, test_MAE, test_PCC, test_SCC, pred_test = validate(model, test_loader, args)
        pred_to_csv(pred_test, ic50_test, dir_pred=args.dir_test)
        
        print('-------- DRPreter -------')
        print(f'sim: {args.sim}')
        print(f'##### Seed: {args.seed} #####')
        print('\t\tMSE\tRMSE\tMAE\tPCC\tSCC')
        print('Test result: {}\t{}\t{}\t{}\t{}'.format(r4(test_MSE), r4(test_RMSE), r4(test_MAE), r4(test_PCC), r4(test_SCC)))
        
        
        '''GradCAM'''
        # ----- (1) Calculate gradient-based importance score for one cell line-drug pair -----        
        # drug_name, cell_name = 'Dihydrorotenone', 'ACH-001374'
        # gradcam_path =  get_path(args, rpath + 'GradCAM/', result_type=f'{drug_name}_{cell_name}_gradcam', extension='csv')
        
        # gene_dict = np.load(rpath + 'Data/Cell/cell_idx2gene_dict.npy', allow_pickle=True)
        
        # # Save importance score
        # sorted_cell_node_importance, indices = gradcam(model, drug_name, cell_name, drug_data, cell_data, edge_index, args)
        # idx2gene = [gene_dict[idx] for idx in indices]
        
        # sorted_cell_node_importance = list(sorted_cell_node_importance.cpu().detach().numpy())
        # indice = list(indices)
        
        # df = pd.DataFrame((zip(sorted_cell_node_importance, indice, idx2gene)), columns=['cell_node_importance','indice','idx2gene'])
        # # df.to_csv(gradcam_path, index=False)
        # print(*list(df['idx2gene'])[:30])
        
        # ----- (2) Calculate scores from total test set in 'inference.csv' -----
        # data = pd.read_excel(f'inference_seed{args.seed}.xlsx', sheet_name='test')
        # name = data[['Drug name', 'DepMap_ID']]
        
        # gene_dict = np.load(rpath + 'Data/Cell/cell_idx2gene_dict.npy', allow_pickle=True)
        
        # total_gene_df = pd.Series(list(range(len(data))))
        # for i in tqdm(range(len(data))):
        #     drug_name, cell_name = name.iloc[i]
        #     _, indices = gradcam(model, drug_name, cell_name, drug_data, cell_data, edge_index, args)
        #     idx2gene = [gene_dict[idx] for idx in indices]
        #     gene_df = pd.DataFrame(idx2gene)
        #     total_gene_df.loc[i] = ', '.join(list(gene_df.drop_duplicates(keep='first')[0])[:5])
        
        # data['Top5 genes'] = total_gene_df
        # data.to_excel(f'inference_seed{args.seed}_gradcam.xlsx', sheet_name='test')
        
            
        '''Visualize pathway-drug self-attention score from Transformer'''
        # ----- (1) For one cell line - drug pair -----
        
        # drug_name, cell_name = 'Rapamycin', 'ACH-000019'

        # # print(cell_name)
        # attn_score = attention_score(model, drug_name, cell_name, drug_data, cell_data, edge_index, args)
        # print(f'attn_score: {attn_score}')
        # print(f'attn_score.shape: {attn_score.shape}') # attn_score.shape: torch.Size([1, 35, 35])
        # # print(torch.sum(attn_score, axis=1))
        # with open(rpath+'Data/Cell/34pathway_score990.pkl', 'rb') as file:
        #     pathway_names = pickle.load(file).keys()
        # tks = [p[5:] for p in list(pathway_names)]
        # tks.append(drug_name)
        # # print(tks)
        # draw_pair_heatmap(attn_score, drug_name, cell_name, tks, args)
        
        # ----- (2) Heatmap of all cell lines of one drug -----
        
        # drug_name = 'Rapamycin'
        # data = pd.read_csv(f'./Data/{drug_name}.csv')
        # cell_list = list(data['DepMap_ID'])
        
        # result_dict = {}
        # total_result = np.full(35, 0)
        # for cell_name in tqdm(cell_list):
        #     attn_score = attention_score(model, drug_name, cell_name, drug_data, cell_data, edge_index, args)
        #     print(attn_score.shape)
        #     attn_score = torch.squeeze(attn_score).cpu().detach().numpy()
        #     print(np.sum(attn_score, axis=1))
        #     result_dict[cell_name] = attn_score[-1, :] # (35, 1)
        #     total_result = np.vstack([total_result, attn_score[-1, :]])
        
        # with open(rpath+'Data/Cell/34pathway_score990.pkl', 'rb') as file:
        #     pathway_names = pickle.load(file).keys()
        # xtks = [p[5:] for p in list(pathway_names)]
        # xtks.append(drug_name)
        # total_result = total_result[1:,:-1]
        # draw_drug_heatmap(total_result, drug_name, xtks, cell_list, args)
        

        '''Interpolation of unknown values'''
#         inference(model, drug_data, cell_data, edge_index, f'inference_seed{args.seed}.xlsx', args)
        
        
if __name__ == "__main__":
    main()

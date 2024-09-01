#!/usr/bin/env python

import os
import torch
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import RobustScaler, StandardScaler


def z_cap(data, max_abs=10, verbose=False) :
    
    if verbose :
        cond_min = sum(sum(data>max_abs))
        cond_max = sum(sum(data<-max_abs))
        
        print("### z-score capping")
        print("# z-score > {} : {}".format(max_abs, cond_min))
        print("# z-score < {} : {}".format(-max_abs, cond_max))
  
    data[data > max_abs] = max_abs
    data[data < -max_abs] = -max_abs
    return(data)


def scale_cell_data(cell_data, scaler=None, robust=True, max_abs=10) :
  
    if scaler is not None :
        cell_data_ = scaler.transform(cell_data)
    else :
        scaler = RobustScaler(unit_variance=True) if robust else StandardScaler() 
        cell_data_ = scaler.fit_transform(cell_data)
    
    cell_data = pd.DataFrame(cell_data_, index=cell_data.index, columns=cell_data.columns)
    cell_data = z_cap(cell_data, max_abs=max_abs)
    return cell_data, scaler


def process_cell(cell_data, path_dict, net_data, edge_type_dict=None,
                 col_src="Pathway1", col_dest="Pathway2", col_type="Edge_Type") : 
    
    cell_data_dict = {}
    from torch_geometric.data import Data
    
    def union_list(list_a, list_b) :
        list_ab = list(set(list_a) & set(list_b)) if list_a!=[] else list(set(list_b))
        return list_ab
    
    edge_type = None
    net_data[col_src] = [path_dict[key] for key in net_data[col_src]]
    net_data[col_dest] = [path_dict[key] for key in net_data[col_dest]]
    cell_data.columns = [path_dict[key] for key in cell_data.columns]
    
    if hasattr(net_data, col_type) :
        if edge_type_dict is None :
            edge_type_list = net_data[col_type].unique()
            edge_type_dict = dict(zip(edge_type_list, range(len(edge_type_list))))
        edge_type = np.array([edge_type_dict[k] for k in net_data[col_type]])
        edge_type = torch.from_numpy(edge_type)
    
    for cell in cell_data.index :
        cell_data_tp = cell_data.loc[[cell]].T
        node_feats = torch.from_numpy(cell_data_tp.values).float()
        edge_index = torch.from_numpy(net_data[[col_src, col_dest]].T.values)
        
        cell_data_tp = Data(x=node_feats, edge_index=edge_index, edge_type=edge_type)
        cell_data_dict[cell] = cell_data_tp
    
    return cell_data_dict, edge_type_dict


def process_cell_linear(cell_data) :
    
    cell_data_dict = {}
    for cell in cell_data.index :
        cell_data_tp = cell_data.loc[[cell]].values[0]
        cell_data_dict[cell] = torch.from_numpy(cell_data_tp).float()

    return cell_data_dict


def parse_cell() :
    import argparse
    parser = argparse.ArgumentParser()
    input_net_ = "data/net_data_biocarta/KNN5_STR9_Reg_Corr.csv"
    input_omics_ = "processed/cell_data_biocarta/SANGER_RNA_GSVA.csv"
    output_dir_ = "processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle"

    parser.add_argument("-undirect", type=int, default=0)
    parser.add_argument("-train", type=str, default=None)
    parser.add_argument("-net", type=str, default=input_net_)
    parser.add_argument("-omics", type=str, default=input_omics_)
    parser.add_argument("-out", type=str, default=output_dir_)
    return parser.parse_args()


if __name__ == '__main__' :
    args = parse_cell()
    input_train = args.train
    input_net = args.net
    input_omics = args.omics
    output_dir = args.out
    
    scaler = None
    net_data = None
    path_dict = None
    edge_type_dict = None
    input_net = input_net if input_net!="None" else None
    cell_data = pd.read_csv(input_omics, header=0, index_col=0, sep=",")
    
    
    # Call network data
    if input_net is not None :
        net_data = pd.read_csv(input_net, header=0, sep=",")
        print("# Graph Edges [Direct] : {}".format(net_data.shape[0]))
    
    if input_net is not None and args.undirect!=0 :
        columns = {"Pathway1":"Pathway2", "Pathway2":"Pathway1"}
        net_data_rev = net_data.copy().rename(columns=columns)
        net_data = pd.concat([net_data, net_data_rev])
        
        if not hasattr(net_data, "Edge_Type") :
            col_sub = ["Pathway1", "Pathway2"]
        else :
            col_sub = ["Pathway1", "Pathway2", "Edge_Type"]
        
        net_data = net_data.drop_duplicates(subset=col_sub, keep="first")
        print("# Graph Edges [Undirect] : {}".format(net_data.shape[0]))
    
    
    # Column alignment & Feature scaler
    # If train data exists, apply those 2 proceeses to external data
    if input_train is not None :
        with open(input_train, "rb") as f :
            _, path_dict, scaler, edge_type_dict = pickle.load(f)
        print("# Train data exists...")
        print("# Column alignment and Scaler from train data are applied...")
    
    if path_dict is not None :
        cell_data = cell_data[list(path_dict.keys())]
    else :
        path_dict = dict(zip(cell_data.columns, range(len(cell_data.columns))))
    
    cell_data, scaler = scale_cell_data(cell_data, scaler, robust=True)
    
    
    # Torch dictionary
    if net_data is not None : 
        cell_data, edge_type_dict = process_cell(cell_data, path_dict, net_data, edge_type_dict)
    else :
        cell_data = process_cell_linear(cell_data)
    
    
    # Save files
    with open(output_dir, "wb") as f :
        pickle.dump([cell_data, path_dict, scaler, edge_type_dict], f)
    
    print("### Total cells {} & pathways {}!".format(len(cell_data), len(path_dict)))
    print("### Cell processing succeeded!\n")
    

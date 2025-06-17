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


def perturb_graph(net_data, col_src="Pathway1", col_dest="Pathway2", perturb_seed=2021) :
    import networkx as nx
    G = nx.Graph()
    G.add_edges_from(net_data[[col_src, col_dest]].values)

    perturbed_G = nx.configuration_model([d for _, d in G.degree()], seed=perturb_seed)
    perturbed_G = nx.Graph(perturbed_G)
    
    perturbed_data = np.array(list(perturbed_G.edges))
    perturbed_data = pd.DataFrame(perturbed_data, columns=[col_src, col_dest])
    
    verbose = True
    if verbose :
        degree = [d for _, d in G.degree()]
        degree_pert = [d for _, d in perturbed_G.degree()]
        
        info_deg1 = degree==degree_pert
        info_deg2 = sum([i==j for i, j in zip(degree, degree_pert)])
        print(f"# Perturbed graph has the same degree? [{info_deg1}, {info_deg2}]")
        
    return perturbed_data


def process_cell(cell_data, path_dict, net_data, edge_type_dict=None,
                 col_src="Pathway1", col_dest="Pathway2", 
                 col_type="Edge_Type", col_weight=None, perturb_seed=None) : 
    
    cell_data_dict = {}
    from torch_geometric.data import Data
    # col_weight = ["Weight_STR9", "Weight_Reg", "Weight_Corr", "Weight_OVL"]
    
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
    
    if perturb_seed >= 0 :
        if edge_type is not None :
            perturbed_edges = []
            for etype, idx in edge_type_dict.items():
                sub_net_data = net_data[net_data[col_type] == etype].copy()
                perturbed_sub_data = perturb_graph(sub_net_data, col_src, col_dest, perturb_seed)
                perturbed_sub_data[col_type] = etype
                perturbed_edges.append(perturbed_sub_data)
            
            net_data = pd.concat(perturbed_edges, ignore_index=True)
            edge_type = np.array([edge_type_dict[k] for k in net_data[col_type]])
            edge_type = torch.from_numpy(edge_type)
        else:
            net_data = perturb_graph(net_data, col_src, col_dest, perturb_seed)
    
    for cell in cell_data.index :
        cell_data_tp = cell_data.loc[[cell]].T
        node_feats = torch.from_numpy(cell_data_tp.values).float()
        edge_index = torch.from_numpy(net_data[[col_src, col_dest]].T.values)
        
        cell_data_tp = Data(x=node_feats, edge_index=edge_index, edge_type=edge_type)
        if col_weight is not None and perturb_seed < 0 :
            cell_data_tp["edge_attr"] = torch.from_numpy(net_data[col_weight].values).float()
            if net_data[col_weight].ndim==1 : 
                cell_data_tp["edge_attr"] = cell_data_tp["edge_attr"].unsqueeze(-1)
        
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
    parser.add_argument("-perturb_seed", type=int, default=-1)
    parser.add_argument("-train", type=str, default=None)
    parser.add_argument("-net", type=str, default=input_net_)
    parser.add_argument("-omics", type=str, default=input_omics_)
    parser.add_argument("-out", type=str, default=output_dir_)
    parser.add_argument("-col1", type=str, default="Pathway1")
    parser.add_argument("-col2", type=str, default="Pathway2")
    parser.add_argument("-col_etype", type=str, default="Edge_Type")
    parser.add_argument("-col_weight", type=str, default=None)
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
    cell_data.columns = cell_data.columns.astype(str)

    
    # Call network data
    if input_net is not None :
        net_data = pd.read_csv(input_net, header=0, sep=",")
        net_data[args.col1] = net_data[args.col1].astype(str)
        net_data[args.col2] = net_data[args.col2].astype(str)
        print("# Graph Edges [Direct] : {}".format(net_data.shape[0]))
    
    if input_net is not None and args.undirect!=0 :
        # columns = {"Pathway1":"Pathway2", "Pathway2":"Pathway1"}
        columns = {args.col1:args.col2, args.col2:args.col1}
        net_data_rev = net_data.copy().rename(columns=columns)
        net_data = pd.concat([net_data, net_data_rev])
        
        if not hasattr(net_data, args.col_etype) :
            # col_sub = ["Pathway1", "Pathway2"]
            col_sub = [args.col1, args.col2]
        else :
            # col_sub = ["Pathway1", "Pathway2", "Edge_Type"]
            col_sub = [args.col1, args.col2, args.col_etype]
        
        if hasattr(net_data, args.col_weight) :
            # col_sub = ["Pathway1", "Pathway2", "Weight", ("Edge_Type")]
            col_sub = [args.col1, args.col2, args.col_etype, args.col_weight]
        
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
        # TODO Fill NAs for external data if exists...
        cell_data = cell_data[list(path_dict.keys())]
    else :
        if net_data is not None:
            path_int = set(net_data[args.col1]) | set(net_data[args.col2])
            path_int = path_int & set(cell_data.columns)
            cond = cell_data.columns.isin(path_int)
            cell_data = cell_data.loc[:, cond]
        path_dict = dict(zip(cell_data.columns, range(len(cell_data.columns))))
    
    cell_data, scaler = scale_cell_data(cell_data, scaler, robust=True)
    
    
    # Torch dictionary
    if net_data is not None : 
        cell_data, edge_type_dict = process_cell(cell_data, path_dict, net_data, edge_type_dict, 
                                                 col_src=args.col1, col_dest=args.col2, 
                                                 col_type=args.col_etype, col_weight=args.col_weight,
                                                 perturb_seed=args.perturb_seed)
    else :
        cell_data = process_cell_linear(cell_data)
    
    
    # Save files
    with open(output_dir, "wb") as f :
        pickle.dump([cell_data, path_dict, scaler, edge_type_dict], f)
    
    print("### Total cells {} & pathways {}!".format(len(cell_data), len(path_dict)))
    print("### Cell processing succeeded!\n")

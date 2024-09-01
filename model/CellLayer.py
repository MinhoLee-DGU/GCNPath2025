#!/usr/bin/env python

import torch
from torch.nn import ReLU, ELU, LeakyReLU
from torch.nn import Linear, Module, ModuleList

from utils.utils_gnn import *


class CellLayer(Module) :
    def __init__(self, gnn="GCN", mode=None, dim=[64, 64, 64], 
                 concat_gat=True, h_gat=4, n_rel=1,
                 n_bases=None, n_blocks=None, is_sorted=True,
                 act=ReLU(), norm=True, drop=0.2, norm_graph=False, **kwargs) :
                     
        super().__init__()
        self.gnn = gnn
        
        if gnn=="linear" :
            self.cell_layer = MLP(dim=dim, act=act, norm=norm, drop=drop)
        else :
            self.cell_layer = GNNLayer(mode=mode, gnn=gnn, dim=dim, act=act, norm=norm, drop=drop, 
                                       norm_graph=norm_graph, heads=h_gat, concat=concat_gat, edge_dim=None, 
                                       n_rel=n_rel, n_bases=n_bases, n_blocks=n_blocks, is_sorted=is_sorted, **kwargs)
        
        self.dim_fin = self.cell_layer.dim_fin
                                       
    def forward(self, data, edge_index=None, edge_type=None, **kwargs) :
            
        if self.gnn.lower() == "linear" :
            # Linear
            data = self.cell_layer(data)
            # (batch, dim[-1])
            return data
            
        else :
            # GNN
            data = self.cell_layer(data, edge_index, edge_attr=None, edge_type=edge_type)
            # (batch*n_path, dim[-1])
            return data
            

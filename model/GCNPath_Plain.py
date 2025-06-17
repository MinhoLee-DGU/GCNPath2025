#!/usr/bin/env python

import torch
from torch.nn import ModuleList, Module
from torch.nn import Linear, BatchNorm1d, Dropout, LayerNorm
from torch.nn import Tanh, ReLU, ELU, LeakyReLU, Mish, Softmax

import torch_geometric
from torch_geometric.nn import global_max_pool, global_mean_pool

from utils.utils_gnn import *
from model.CellLayer import *
from model.DrugLayer import *


class GCNPath(Module) :
    
    def __init__(self, gnn_cell="GCN", gnn_drug="GCN", mode_cell=None, mode_drug=None,
                 dim_cell=[16, 16], dim_drug=[64, 64], dim_pred=[128, 128], 
                 n_drug_edge=None, concat_gat=True, h_gat_cell=4, h_gat_drug=4, n_rel=1, 
                 attn_mode=0, dim_attn=16, dim_embed_cell=256, dim_embed_drug=256, dim_compact=32,
                 n_path=32, coef_ffnn=2, n_attn_layer=1, h_attn=4, act=ReLU(), norm=True, drop=0.2) :
        
        super().__init__()
        self.n_path = n_path
        self.n_subs = len(dim_drug[1:])
        self.dim_cell = dim_cell
        self.dim_drug = dim_drug
        
        self.dim_attn = dim_attn
        self.dim_embed_cell = dim_embed_cell
        self.dim_embed_drug = dim_embed_drug
        
        self.attn_mode = attn_mode
        self.gnn_cell = gnn_cell.lower() != "linear"
        self.gnn_drug = gnn_drug.lower() != "linear"
        
        # Cell & Drug [GNN or Linear]
        # Do not apply dropout in GNN...
        drop_cell = False if self.gnn_cell else drop
        drop_drug = False if self.gnn_drug else drop
        norm_cell = True if self.gnn_cell and norm not in ["batch", "True", True] else False
        norm_drug = True if self.gnn_drug and norm not in ["batch", "True", True] else False

        self.cell_layer = CellLayer(gnn=gnn_cell, mode=mode_cell, dim=dim_cell, 
                                    concat_gat=concat_gat, h_gat=h_gat_cell, 
                                    n_rel=n_rel, n_bases=None, n_blocks=None,
                                    act=act, norm=norm, drop=drop_cell, norm_graph=norm_cell)
        
        self.drug_layer = DrugLayer(gnn=gnn_drug, mode=mode_drug, dim=dim_drug, 
                                    concat_gat=concat_gat, h_gat=h_gat_drug, n_drug_edge=n_drug_edge, 
                                    attn_mode=attn_mode, act=act, norm=norm, drop=drop_drug, norm_graph=norm_drug)
                                        
        print("\n### Model Information")
        print("# Cell Layer : {} {}".format(gnn_cell.upper(), dim_cell))
        print("# Drug Layer : {} {}".format(gnn_drug.upper(), dim_drug))
        
        self.dim_cell_fin = self.cell_layer.dim_fin
        self.dim_drug_fin = self.drug_layer.dim_fin
        
        # Cell Compact Layer
        if self.gnn_cell :
            self.dim_compact = min(self.dim_cell_fin, dim_compact)
            dim_compact_cell = [self.dim_cell_fin, self.dim_compact]
            self.cell_compact_layer = MLP(dim=dim_compact_cell, act=act, norm="layer", drop=drop)
        
        # Cell & Drug Embedding Layer
        if self.gnn_cell :
            dim_embed_cell = [self.n_path * self.dim_compact, self.dim_embed_cell]
        else :
            dim_embed_cell = [self.dim_cell_fin, self.dim_embed_cell]
        
        dim_embed_drug = [self.dim_drug_fin, self.dim_embed_drug]
        self.cell_embed_layer = MLP(dim=dim_embed_cell, act=act, norm=True, drop=drop)
        self.drug_embed_layer = MLP(dim=dim_embed_drug, act=act, norm=True, drop=drop)
        
        # Final Layer for IC50
        dim_pred = [self.dim_embed_cell+self.dim_embed_drug, *dim_pred, 1]
        self.final_layer = MLP(dim=dim_pred, act=act, norm=True, drop=drop, last_nn=True)
    
    
    def summary(self, verbose=True) :
        if verbose : 
            print("\n### Model Summary\n{}\n".format(self))
    
        param_total = sum(p.numel() for p in self.parameters())
        param_cell = sum(p.numel() for p in self.cell_layer.parameters())
        param_drug = sum(p.numel() for p in self.drug_layer.parameters())
        
        print("# Model Parameters : {}".format(param_total))
        print("# Model Parameters : {} [Cell]".format(param_cell))
        print("# Model Parameters : {} [Drug]".format(param_drug))
        
        if hasattr(self, "attention") :
            param_attn_cell = sum(p.numel() for p in self.cell_pre_attn_layer.parameters())
            param_attn_drug = sum(p.numel() for p in self.drug_pre_attn_layer.parameters())
            param_attn = sum(p.numel() for p in self.attention.parameters())
            print("# Model Parameters : {} [Attn]".format(param_attn+param_attn_cell+param_attn_drug))
        
        if hasattr(self, "cell_compact_layer") :
            param_compact = sum(p.numel() for p in self.cell_compact_layer.parameters())
            print("# Model Parameters : {} [Cell_Compact]".format(param_compact))
                    
        if hasattr(self, "concat_embed_layer") :
            param_embed_cat = sum(p.numel() for p in self.concat_embed_layer.parameters())
            print("# Model Parameters : {} [Concat_Embed]".format(param_embed_cat))
        
        if all(hasattr(self, _) for _ in ["cell_embed_layer", "drug_embed_layer"]) :
            param_embed_cell = sum(p.numel() for p in self.cell_embed_layer.parameters())
            param_embed_drug = sum(p.numel() for p in self.drug_embed_layer.parameters())
            print("# Model Parameters : {} [Cell_Embed]".format(param_embed_cell))
            print("# Model Parameters : {} [Drug_Embed]".format(param_embed_drug))
        
        param_final = sum(p.numel() for p in self.final_layer.parameters())
        print("# Model Parameters : {} [Final]".format(param_final))
    
    
    def forward(self, cell_data, drug_data, return_attn=False, grad_cam=False) :
        
        cell = cell_data
        batch_cell = None
        edge_index_cell = None
        edge_type_cell = None
        edge_attr_cell = None
        
        drug = drug_data
        batch_drug = None
        edge_index_drug = None
        edge_attr_drug = None
        
        if self.gnn_cell :
            cell = cell_data.x
            batch_cell = cell_data.batch
            edge_index_cell = cell_data.edge_index
            edge_attr_cell = cell_data.edge_attr
            
            if hasattr(cell_data, "edge_type") :
                edge_type_cell = cell_data.edge_type
        
        if self.gnn_drug :
            drug = drug_data.x
            batch_drug = drug_data.batch
            edge_index_drug = drug_data.edge_index
            edge_attr_drug = drug_data.edge_attr
        
        # Cell Layer
        cell = self.cell_layer(cell, edge_index_cell, edge_type=edge_type_cell, 
                               edge_attr=edge_attr_cell, grad_cam=grad_cam)
        # Linear    [batch, dim_cell_fin]
        # GNN       [batch*n_path, dim_cell_fin]
        
        # Grad-CAM
        if grad_cam :
            cell, cell_explain = cell
            # [batch*n_path, dim_cell_fin, [batch*n_path, dim_cell[1]]
        
        # Drug Layer
        drug = self.drug_layer(drug, edge_index_drug, edge_attr=edge_attr_drug, batch=batch_drug)
        # Linear    [batch, dim_drug_fin]
        # GNN       [batch*n_atom, dim_drug_fin]
        
        # Plain Mode without Attention
        if self.gnn_cell :
            cell = self.cell_compact_layer(cell)
            # [batch*n_path, compact_cell]
            cell = cell.view(-1, self.n_path*self.dim_compact)
            # [batch, n_path*compact_cell]
        
        if self.gnn_drug :
            drug = global_max_pool(drug, batch_drug)
            # [batch, dim_drug_fin]
        
        cell = self.cell_embed_layer(cell)
        # [batch, dim_embed_cell]
        drug = self.drug_embed_layer(drug)
        # [batch, dim_embed_drug]
        concat = torch.cat([cell, drug], axis=1)
        # [batch, dim_embed_cell + dim_embed_drug]
        
        # Final Layer
        ic50 = self.final_layer(concat)
        # [batch, 1]
        
        if self.attn_mode!=0 and return_attn :
            return ic50, attn
        elif grad_cam :
            return ic50, cell_explain
        else :
            return ic50


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
from model.Transformer import *


class GCNPath_SSI_v2(Module) :
    
    def __init__(self, gnn_cell="GCN", gnn_drug="GCN", mode_cell=None, mode_drug=None,
                 dim_cell=[16, 16], dim_drug=[64, 64], dim_pred=[128, 128], bond_info=True,
                 n_drug_edge=None, concat_gat=True, h_gat_cell=4, h_gat_drug=4, n_rel=1, 
                 attn_mode=0, dim_attn=16, dim_compact=32, dim_embed_cell=256, dim_embed_drug=256,
                 n_path=32, coef_ffnn=4, n_attn_layer=1, h_attn=4, act=ReLU(), norm=True, drop=0.2) :
        
        super().__init__()
        self.n_path = n_path
        self.n_subs = len(dim_drug[1:])
        
        self.dim_attn = dim_attn
        self.dim_embed_cell = dim_embed_cell
        self.dim_embed_drug = dim_embed_drug
        
        self.attn_mode = attn_mode
        self.gnn_cell = gnn_cell.lower() != "linear"
        self.gnn_drug = gnn_drug.lower() != "linear"
        assert gnn_cell.lower()!="linear", "Attention forces cell encoder to be GNN..."
        
        # Cell & Drug [GNN or Linear]
        # Do not apply dropout in GNN...
        drop_cell = False if self.gnn_cell else drop
        # drop_drug = False if self.gnn_drug else drop
        
        drop_drug = drop
        norm_cell = True if self.gnn_cell and norm not in ["batch", "True", True] else False
        norm_drug = True if self.gnn_drug and norm not in ["batch", "True", True] else False
        
        # Cell & Drug [GNN or Linear]
        self.cell_layer = CellLayer(gnn=gnn_cell, mode=mode_cell, dim=dim_cell, 
                                    concat_gat=concat_gat, h_gat=h_gat_cell, 
                                    n_rel=n_rel, n_bases=None, n_blocks=None,
                                    act=act, norm=norm, drop=drop_cell, norm_graph=norm_cell)
                                    
        self.drug_layer = DrugLayer(gnn=gnn_drug, mode=mode_drug, dim=dim_drug, 
                                    concat_gat=concat_gat, h_gat=h_gat_drug, n_drug_edge=n_drug_edge, 
                                    attn_mode=attn_mode, act=act, norm=norm, drop=drop_drug, norm_graph=norm_drug, dmpnn=False)
        
        print("\n### Model Information")
        print("# Cell Layer : {} {}".format(gnn_cell.upper(), dim_cell))
        print("# Drug Layer : {} {}".format(gnn_drug.upper(), dim_drug))
        
        self.dim_cell_fin = self.cell_layer.dim_fin
        self.dim_drug_fin = self.drug_layer.dim_fin
        
        # Pre-Attention Layer
        dim_pre_attn_cell = [self.dim_cell_fin, dim_attn]
        dim_pre_attn_drug = [self.dim_drug_fin, dim_attn]
        self.cell_pre_attn_layer = MLP(dim=dim_pre_attn_cell, act=act, norm="layer", drop=drop)
        self.drug_pre_attn_layer = MLP(dim=dim_pre_attn_drug, act=act, norm="layer", drop=drop)
        
        # Attention Layer
        dim_ffnn = coef_ffnn * dim_attn
        attention_layer = EncoderLayer(n_heads=h_attn, dim_in=dim_attn, dim_ffnn=dim_ffnn, drop=drop, act_gelu=True)
        self.attention = Encoder(attention_layer, n_layer=n_attn_layer)
        
        # Cell Compact Layer
        self.dim_compact = min(dim_attn, dim_compact)
        dim_compact_cell = [dim_attn, self.dim_compact]
        self.cell_compact_layer = MLP(dim=dim_compact_cell, act=act, norm="layer", drop=drop)
        
        # Cell & Drug Embedding Layer
        dim_embed_cell = [self.n_path * self.dim_compact, self.dim_embed_cell]
        dim_embed_drug = [self.n_subs * (dim_attn + self.dim_drug_fin), self.dim_embed_drug]
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
        
        if all(hasattr(self, _) for _ in ["cell_embed_layer", "drug_embed_layer"]) :
            param_embed_cell = sum(p.numel() for p in self.cell_embed_layer.parameters())
            param_embed_drug = sum(p.numel() for p in self.drug_embed_layer.parameters())
            print("# Model Parameters : {} [Cell_Embed]".format(param_embed_cell))
            print("# Model Parameters : {} [Drug_Embed]".format(param_embed_drug))
        
        param_final = sum(p.numel() for p in self.final_layer.parameters())
        print("# Model Parameters : {} [Final]".format(param_final))
    
    
    def forward(self, cell_data, drug_data, drug_data_bond=None, return_attn=False) :
        
        cell = cell_data
        batch_cell = None
        edge_index_cell = None
        edge_type_cell = None
                
        drug = drug_data
        batch_drug = None
        edge_index_drug = None
        edge_attr_drug = None
        
        if self.gnn_cell :
            cell = cell_data.x
            batch_cell = cell_data.batch
            edge_index_cell = cell_data.edge_index
            
        if self.gnn_cell and hasattr(cell_data, "edge_type") :
            edge_type_cell = cell_data.edge_type
        
        if self.gnn_drug :
            drug = drug_data.x
            batch_drug = drug_data.batch
            edge_index_drug = drug_data.edge_index
            edge_attr_drug = drug_data.edge_attr
        
        
        # Cell Layer
        cell = self.cell_layer(cell, edge_index_cell, edge_type=edge_type_cell)
        # [batch*n_path, dim_cell_fin]
        
        # Drug Layer
        drug = self.drug_layer(drug, edge_index_drug, batch=batch_drug, edge_attr=edge_attr_drug)
        # [batch, n_subs, dim_drug_fin]
        
        # Dual-Attention of Cell-Subs & Drug-Path [SSI Module]
        cell = cell.view(-1, self.n_path, self.dim_cell_fin)
        # [batch, n_path, dim_cell_fin]
        cell = self.cell_pre_attn_layer(cell)
        # [batch, n_path, dim_attn]
        drug_attn = self.drug_pre_attn_layer(drug)
        # [batch, n_subs, dim_attn]
        concat = torch.cat([cell, drug_attn], axis=1)
        # [batch, n_path+n_subs, dim_attn]
        concat, attn = self.attention(concat)
        # [batch, n_path+n_subs, dim_attn], [batch, n_path+n_subs, n_path+n_subs]
        
        # Flatten instead of Global Pooling [Cell, Drug]
        cell = concat[:, :self.n_path, :]
        # [batch, n_path, dim_attn]
        drug_attn = concat[:, self.n_path:, :]
        # [batch, n_subs, dim_attn]
        cell = self.cell_compact_layer(cell)
        # [batch, n_path, dim_compact]
        cell = cell.view(-1, self.n_path*self.dim_compact)
        # [batch, n_path*dim_compact]
        drug_attn = drug_attn.view(-1, self.n_subs*self.dim_attn)
        # [batch, n_subs*dim_attn]

        cell = self.cell_embed_layer(cell)
        # [batch, dim_embed_cell]
        drug = drug.view(-1, self.n_subs*self.dim_drug_fin)
        # [batch, n_subs*dim_drug_fin]
        drug = torch.cat([drug, drug_attn], axis=1)
        # [batch, n_subs*(dim_drug_fin+dim_attn)]
        drug = self.drug_embed_layer(drug)
        # [batch, dim_embed_drug]
        concat = torch.cat([cell, drug], axis=1)
        # [batch, dim_embed_cell + dim_embed_drug]
        
        # Final Layer
        ic50 = self.final_layer(concat)
        # [batch, 1]
        
        if self.attn_mode!=0 and return_attn :
            return ic50, attn
        else :
            return ic50

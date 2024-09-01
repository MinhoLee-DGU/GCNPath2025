#!/usr/bin/env python

import torch
from torch_scatter import scatter
from torch.nn import Tanh, ReLU, ELU, LeakyReLU, PReLU
from torch.nn import Linear, BatchNorm1d, Parameter, Sequential, Module, ModuleList

from torch_geometric.utils import degree
from torch_geometric.nn.inits import glorot
from torch_geometric.nn.norm import LayerNorm
from torch_geometric.nn import GATv2Conv, GINConv, GraphConv, SAGPooling, global_mean_pool, global_add_pool

from utils.utils_gnn import *


class DrugLayer(Module) :
    
    def __init__(self, gnn="GAT", mode=None, dim=[128, 128, 128, 128], 
                 n_drug_edge=None, concat_gat=True, h_gat=4, attn_mode=3, 
                 act=ReLU(), norm=True, drop=0.2, dmpnn=False, iter_dmpnn=10, **kwargs) :
                     
        super().__init__()
        self.gnn = gnn
        self.attn_mode = attn_mode
        
        if self.attn_mode in [3, 4] :
            # Substructure-Aware Mode [SSI-DDI]
            n_layer = len(dim)-1
            use_gin = {"gat":0, "gin":1, "gine":2}
            use_gin = use_gin[gnn.lower()] if gnn.lower() in use_gin.keys() else 0
            self.drug_layer = SSIModule(n_layer=n_layer, dim_in=dim[0], dim_out=dim[1], 
                                        heads=h_gat, concat=concat_gat, edge_dim=n_drug_edge, 
                                        act=act, use_gin=use_gin, dmpnn=dmpnn, iter_dmpnn=iter_dmpnn, **kwargs)
        else :
            if gnn!="linear" :
                self.drug_layer = GNNLayer(mode=mode, gnn=gnn, dim=dim, act=act, norm=norm, 
                                           heads=h_gat, concat=concat_gat, edge_dim=n_drug_edge, **kwargs)
            else :
                self.drug_layer = MLP(dim=dim, act=act, norm=norm, drop=drop, **kwargs)
                
        self.dim_fin = self.drug_layer.dim_fin
        
    def forward(self, data, edge_index=None, batch=None, edge_attr=None, 
                edge_index_bond=None, edge_index_batch=None) :
        
        if self.gnn.lower() == "linear" :
            # Linear
            data = self.drug_layer(data)
            # [batch, dim[-1]]
            return data
        
        elif self.attn_mode not in [3, 4] :
            # GNN
            data = self.drug_layer(data, edge_index, edge_attr=edge_attr)
            # [batch*n_atoms, dim[-1]]
            return data
        
        else :
            # SSI-DDI
            drug, global_emb = self.drug_layer(data, edge_index, batch=batch, edge_attr=edge_attr, 
                                               edge_index_bond=edge_index_bond, edge_index_batch=edge_index_batch)
            # [batch*n_atoms, dim[1]], [batch, n_gat, dim[1]]
            return global_emb


class SSIModule(Module) :
    
    def __init__(self, n_layer=3, dim_in=128, dim_out=128, heads=4, concat=True, 
                 edge_dim=16, act=ELU(1.0), use_gin=0, dmpnn=False, iter_dmpnn=10) :
        
        super().__init__()
        self.dmpnn = dmpnn
        self.dim_fin = dim_out
        self.dim_out = dim_out
        self.layers = ModuleList()
        
        if self.dmpnn :
            self.pre_dmpnn = Sequential(
                Linear(dim_in, dim_out),
                PReLU(),
                Linear(dim_out, dim_out),
                BatchNorm1d(dim_out),
                PReLU(),
                Linear(dim_out, dim_out),
                BatchNorm1d(dim_out), 
            )
            dim_in = dim_out
        
        for _ in range(n_layer) :
            layer = SSILayer(dim_in, dim_out, heads=heads, concat=concat, act=act, 
                             edge_dim=edge_dim, use_gin=use_gin, dmpnn=dmpnn, iter_dmpnn=iter_dmpnn)
            self.layers.append(layer)
            dim_in = dim_out
            
    def forward(self, data, edge_index, batch, edge_attr=None, 
                edge_index_bond=None, edge_index_batch=None) :
        
        drug_embed_list = []
        if self.dmpnn : data = self.pre_dmpnn(data)
        
        for layer in self.layers :
            if not self.dmpnn :
                data, global_emb = layer(data, edge_index, batch=batch, edge_attr=edge_attr)
            else :
                data, global_emb = layer(data, edge_index, batch=batch, edge_attr=edge_attr, 
                                         edge_index_bond=edge_index_bond, edge_index_batch=edge_index_batch)
            # [batch*n_atom, dim_out], [batch, dim_out]
            drug_embed_list.append(global_emb)
            # [batch, dim_out] x n_layer
        
        drug_embed_list = torch.stack(drug_embed_list, dim=1)
        # [batch, n_layer, dim_out]
        return data, drug_embed_list
        

class SSILayer(Module) :
    
    def __init__(self, dim_in=128, dim_out=128, heads=4, concat=True, 
                 act=ELU(1.0), edge_dim=10, use_gin=0, dmpnn=False, iter_dmpnn=10) :
        
        super().__init__()
        self.dmpnn = dmpnn
        self.use_gin = use_gin
        
        if use_gin==0 :
            # [use_gin=0] GATv2Conv
            dim_mid = int(dim_out / heads) if concat else dim_out
            self.conv = GATv2Conv(dim_in, dim_mid, heads=heads, concat=concat, edge_dim=edge_dim)
            if concat and (dim_out % heads)!=0 :
                raise NotImplementedError('[SSI_DDI] GAT dim_out {} & heads {}'.format(dim_out, heads))
        elif use_gin==1 :
            # [use_gin=1] GINConv
            gin_linear = Sequential(Linear(dim_in, dim_out), ReLU(), Linear(dim_out, dim_out))
            self.conv = GINConv(gin_linear)
        else :
            # [use_gin>1] GINEConv
            gin_linear = Sequential(Linear(dim_in, dim_out), ReLU(), Linear(dim_out, dim_out))
            self.conv = GINEConv(gin_linear, edge_dim=edge_dim)
        if self.dmpnn :
            self.conv = DMPNN(dim_in, dim_out, edge_dim=edge_dim, n_iter=iter_dmpnn)
        
        self.act = act
        self.readout = SAGPooling(dim_out, min_score=-1)
        self.layer_norm = torch_geometric.nn.norm.LayerNorm(dim_out)
        self.global_pool = global_add_pool
    
    def forward(self, data, edge_index, batch, edge_attr=None,
                edge_index_bond=None, edge_index_batch=None) :
        
        if self.dmpnn :
            data = self.conv(data, edge_index, edge_attr, edge_index_bond, edge_index_batch)
            # [batch*n_atom, dim_out]
        elif self.use_gin==1 :
            data = self.conv(data, edge_index)
            # [batch*n_atom, dim_out]
        else :
            data = self.conv(data, edge_index, edge_attr=edge_attr)
            # [batch*n_atom, dim_out]
        
        attn_x, _, _, attn_batch, _, _ = self.readout(data, edge_index, batch=batch)
        # attn_x, attn_edge_index, attn_edge_attr, attn_batch, attn_perm, attn_scores
        # [batch*n_subs, dim_out]
        global_emb = self.global_pool(attn_x, attn_batch)
        # (batch, dim_out)
        
        data = self.layer_norm(data, batch=batch)
        # [batch*n_atom, dim_out]
        data = self.act(data)
        # [batch*n_atom, dim_out]
        
        return data, global_emb


class DMPNN(Module) :
    # DMPNN from SA-DDI
    
    def __init__(self, dim_in, dim_out, edge_dim, n_iter) :
        super().__init__()
        self.n_iter = n_iter
        
        self.act = ReLU()
        self.act_gout = Tanh()
        self.lin_u = Linear(dim_in, dim_out, bias=False)
        self.lin_v = Linear(dim_in, dim_out, bias=False)
        self.lin_edge = Linear(edge_dim, dim_out, bias=False)
    
        self.att = GlobalAttentionPool(dim_out)
        self.a = Parameter(torch.zeros(1, dim_out, n_iter))
        self.lin_gout = Linear(dim_out, dim_out)
        self.a_bias = Parameter(torch.zeros(1, 1, n_iter))

        glorot(self.a)
        # self.lin_block = LinearBlock(dim_out)

    def forward(self, x, edge_index, edge_attr, edge_index_bond, edge_index_batch) :
        
        # edge_index_bond : line_graph_edge_index
        
        edge_u = self.lin_u(x)
        # [batch*n_atom, dim_out]
        edge_v = self.lin_v(x)
        # [batch*n_atom, dim_out]
        edge_uv = self.lin_edge(edge_attr)
        # [batch*n_edge, dim_out]
        edge_attr = (edge_u[edge_index[0]] + edge_v[edge_index[1]] + edge_uv) / 3
        # [batch*n_edge, dim_out]
        out = edge_attr
        # [batch*n_edge, dim_out]
        
        out_list = []
        gout_list = []
        for n in range(self.n_iter):
            out = scatter(out[edge_index_bond[0]], edge_index_bond[1], dim_size=edge_attr.size(0), dim=0, reduce='add')
            # [batch*n_edge, dim_out]
            out = edge_attr + out
            # [batch*n_edge, dim_out]
            gout = self.att(out, edge_index_bond, edge_index_batch)
            # [batch, dim_out]
            out_list.append(out)
            # [batch*n_edge, dim_out] x n_iter
            gout_list.append(self.act_gout(self.lin_gout(gout)))
            # [batch, dim_out] x n_iter
        
        ### Substructure attention [Eq. 3]
        gout_all = torch.stack(gout_list, dim=-1)
        # [batch, dim_out, n_iter]
        out_all = torch.stack(out_list, dim=-1)
        # [batch*n_edge, dim_out, n_iter]
        
        ### Substructure attention [Eq. 4]
        scores = (gout_all * self.a).sum(1, keepdim=True) + self.a_bias
        # [batch, 1, n_iter]
        scores = torch.softmax(scores, dim=-1)
        # [batch, 1, n_iter]
        
        ### Weighted sum of bond-level hidden features across all steps [Eq. 5]
        scores = scores.repeat_interleave(degree(edge_index_batch, dtype=edge_index_batch.dtype), dim=0)
        # [batch*n_edge, 1, n_iter]
        out = (out_all * scores).sum(-1)
        # [batch*n_edge, dim_out]
        
        ### Return to node-level hidden features [Eq. 6, 7]
        x = x + scatter(out, edge_index[1], dim_size=x.size(0), dim=0, reduce='add')
        # [batch*n_atom, dim_out]
        # x = self.lin_block(x)
        # [batch*n_atom, dim_out]
        
        return x
        

class GlobalAttentionPool(Module) :
    
    def __init__(self, hidden_dim) :
        super().__init__()
        self.conv = GraphConv(hidden_dim, 1)

    def forward(self, x, edge_index, batch) :
        x_conv = self.conv(x, edge_index)
        scores = torch_geometric.utils.softmax(x_conv, batch, dim=0)
        gx = global_add_pool(x * scores, batch)

        return gx

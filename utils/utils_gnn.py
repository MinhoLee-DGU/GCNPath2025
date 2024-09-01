#!/usr/bin/env python

import numpy

import torch
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint
from torch.nn import ModuleList, Module, Sequential
from torch.nn import Tanh, ReLU, ELU, LeakyReLU, Softmax
from torch.nn import Linear, Dropout, BatchNorm1d, LayerNorm, InstanceNorm1d

import torch_geometric
from torch_geometric.nn import global_max_pool, global_mean_pool
from torch_geometric.nn.models import DeepGCNLayer, JumpingKnowledge
from torch_geometric.nn.conv import GCNConv, GATv2Conv, GINConv, GINEConv
from torch_geometric.nn.conv import RGCNConv, FastRGCNConv, RGATConv


def norm_layer(norm_type, dim) :
    if norm_type in ["batch", "True", True] :
        layer = BatchNorm1d(dim)
    elif norm_type == "layer" :
        layer = LayerNorm(dim)
    elif norm_type == "instance" :
        layer = InstanceNorm1d(dim)
    elif norm_type in ["none", "None", "False", None, False] :
        layer = None
    else:
        raise NotImplementedError('normalization layer [%s] is not found' % norm_type)
    return layer


def norm_layer_graph(norm_type, dim, **kwargs) :
    # forward(x, batch, batch_size)
    if norm_type == "layer" :
        layer = torch_geometric.nn.norm.LayerNorm(dim, **kwargs)
    elif norm_type == "instance" :
        layer = torch_geometric.nn.norm.InstanceNorm(dim, **kwargs)
    elif norm_type == "graph" :
        layer = torch_geometric.nn.norm.GraphNorm(dim, **kwargs)
    elif norm_type == "graph_size" :
        layer = torch_geometric.nn.norm.GraphSizeNorm()
    elif norm_type in ["none", "None", "False", None, False] :
        layer = None
    else:
        raise NotImplementedError('normalization layer [%s] is not found' % norm_type)
    return layer


def global_pool(x, batch_data) :
    # x : (batch x n_nodes, feat)
    max_pool = global_max_pool(x, batch_data)
    # (batch, feat)
    mean_pool = global_mean_pool(x, batch_data)
    # (batch, feat)
    pool = torch.cat([max_pool, mean_pool], axis=1)
    # (batch, 2*feat)
    return(pool)


class GNN(Module) :

    def __init__(self, dim_in=16, dim_out=16, gnn='GCN', 
                 edge_dim=None, heads=4, concat=True, 
                 n_rel=1, n_bases=None, n_blocks=None, is_sorted=True, **kwargs) :
        super(GNN, self).__init__()
        self.gnn_ = gnn.lower()

        if self.gnn_ in ["gin", "gine"] :
            gin_linear = Sequential(Linear(dim_in, dim_out), ReLU(), Linear(dim_out, dim_out))
                
        if self.gnn_ in ["gat", "rgat"] :
            if concat and (dim_out % heads)!=0 :
                raise NotImplementedError("GAT dim_out {} & heads {}".format(dim_out, heads))
            elif concat :
                dim_out_gat = int(dim_out/heads)
            else : 
                dim_out_gat = dim_out
            
        if self.gnn_ in ["rgcn", "rgat"] and n_blocks is not None :
            if (dim_in % n_blocks)!=0 or (dim_out % n_blocks)!=0 : n_blocks = None

        if self.gnn_=="gcn" :
            self.gnn = GCNConv(dim_in, dim_out, **kwargs)
        elif self.gnn_=="gat" :
            self.gnn = GATv2Conv(dim_in, dim_out_gat, edge_dim=edge_dim, heads=heads, concat=concat, **kwargs)
        elif self.gnn_=="gin" :
            self.gnn = GINConv(gin_linear, **kwargs)
        elif self.gnn_=="gine" :
            self.gnn = GINEConv(gin_linear, edge_dim=edge_dim, **kwargs)
        elif self.gnn_=="rgcn" :
            self.gnn = RGCNConv(dim_in, dim_out, num_relations=n_rel, 
                                num_bases=n_bases, num_blocks=n_blocks, is_sorted=is_sorted, **kwargs)
        elif self.gnn_=="rgat" :
            self.gnn = RGATConv(dim_in, dim_out_gat, num_relations=n_rel, 
                                num_bases=n_bases, num_blocks=n_blocks, is_sorted=is_sorted, 
                                edge_dim=edge_dim, heads=heads, concat=concat, **kwargs)
        else :
            raise NotImplementedError('GNN {} is not found...'.format(gnn))

    def forward(self, x, edge_index, edge_attr=None, edge_type=None) :
        if self.gnn_ in ["gcn", "gin"] :
            x = self.gnn(x, edge_index)
        elif self.gnn_ in ["gat", "gine"] :
            x = self.gnn(x, edge_index, edge_attr=edge_attr)
        elif self.gnn_ == "rgcn" :
            x = self.gnn(x, edge_index, edge_type=edge_type)
        elif self.gnn_ == "rgat" :
            x = self.gnn(x, edge_index, edge_attr=edge_attr, edge_type=edge_type)
        else :
            raise NotImplementedError('GNN {} is not found...'.format(self.gnn_))
        return x


class DeepGCNLayerDef(Module) :
    
    def __init__(self, conv, norm=None, act=None, block='res+', dropout=0.0, ckpt_grad=False) :
        super().__init__()
        self.act = act
        self.conv = conv
        self.norm = norm
        self.dropout = dropout
        self.ckpt_grad = ckpt_grad
        
        norm_graph = (
            torch_geometric.nn.norm.LayerNorm, 
            torch_geometric.nn.norm.InstanceNorm, 
            torch_geometric.nn.norm.GraphNorm, 
            torch_geometric.nn.norm.GraphSizeNorm
        )
        
        self.block = block.lower()
        assert self.block in ['res+', 'res', 'dense', 'plain']
        self.norm_graph = isinstance(norm, norm_graph)
    
    def reset_parameters(self):
        self.conv.reset_parameters()
        if self.norm is not None :
            self.norm.reset_parameters()

    def forward(self, *args, **kwargs) :
        args = list(args)
        x = args.pop(0)

        if self.block=="res+" :
            h = x
            if self.norm is not None :
                if not self.norm_graph :
                    h = self.norm(h)
                else :
                    h = self.norm(h, batch=kwargs.get("batch"))
            if self.act is not None :
                h = self.act(h)
            h = F.dropout(h, p=self.dropout, training=self.training)
            if self.conv is not None and self.ckpt_grad and h.requires_grad :
                h = checkpoint(self.conv, h, *args, use_reentrant=True, **kwargs)
            else :
                h = self.conv(h, *args, **kwargs)

            return x + h

        else:
            if self.conv is not None and self.ckpt_grad and x.requires_grad :
                h = checkpoint(self.conv, x, *args, use_reentrant=True, **kwargs)
            else :
                h = self.conv(x, *args, **kwargs)
            if self.norm is not None :
                if not self.norm_graph :
                    h = self.norm(h)
                else :
                    h = self.norm(h, batch=kwargs.get("batch"))
            if self.act is not None :
                h = self.act(h)

            if self.block=="res" :
                h = x + h
            elif self.block=="dense" :
                h = torch.cat([x, h], dim=-1)
            elif self.block=="plain" :
                pass

            return F.dropout(h, p=self.dropout, training=self.training)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(block={self.block})'


class GNNLayer(Module) :
  
    def __init__(self, dim=[16, 16, 16, 16], gnn="gcn", 
                 mode=None, act=ReLU(), norm="batch", drop=0.2, 
                 edge_dim=None, heads=4, concat=True, last_nn=False, 
                 n_rel=1, n_bases=None, n_blocks=None, is_sorted=True, norm_graph=False, **kwargs) :
        
        super(GNNLayer, self).__init__()
        self.gnn = gnn
        module_layer = []
        self.mode = mode if mode is not None else "plain"
        mode_ = self.mode if self.mode not in ["jk_max", "jk_cat"] else "plain"
        norm_layer_def = norm_layer if not norm_graph else norm_layer_graph
        
        gnn_1st = GNN(dim_in=dim[0], dim_out=dim[1], gnn=self.gnn,
                      edge_dim=edge_dim, heads=heads, concat=concat, n_rel=n_rel, 
                      n_bases=n_bases, n_blocks=n_blocks, is_sorted=is_sorted, **kwargs)
        
        norm_1st = norm_layer_def(norm, dim=dim[1])
        module_layer.append(DeepGCNLayerDef(conv=gnn_1st, norm=norm_1st, act=act, dropout=drop, block="plain", **kwargs))
        
        if self.mode in ["plain", "jk_max", "jk_cat"] :
            dim_in = dim[1:-1]
            dim_out = dim[2:]
        elif self.mode in ["res", "res+"] :
            dim_in = [dim[1] for i in range(len(dim)-2)]
            dim_out = [dim[1] for i in range(len(dim)-2)]
        else :
            dim_in = [dim[1]*(i+1) for i in range(len(dim)-2)]
            dim_out = [dim[1] for i in range(len(dim)-2)]
        
        for i in range(len(dim_in)) :
            gnn = GNN(dim_in=dim_in[i], dim_out=dim_out[i], gnn=self.gnn,
                      edge_dim=edge_dim, heads=heads, concat=concat, n_rel=n_rel, 
                      n_bases=n_bases, n_blocks=n_blocks, is_sorted=is_sorted, **kwargs)
            
            dim_norm = dim_in[i] if self.mode=="res+" else dim_out[i]
            norm_ = norm_layer_def(norm, dim=dim_norm)
            module_layer.append(DeepGCNLayerDef(conv=gnn, norm=norm_, act=act, dropout=drop, block=mode_, **kwargs))
        
        if last_nn :
            del module_layer[-1]
            gnn_last = GNN(dim_in=dim_in[-1], dim_out=dim_out[-1], gnn=self.gnn,
                           edge_dim=edge_dim, heads=heads, concat=concat, n_rel=n_rel, 
                           n_bases=n_bases, n_blocks=n_blocks, is_sorted=is_sorted, **kwargs)
            module_layer.append(DeepGCNLayerDef(conv=gnn, norm=None, act=None, dropout=None, block=mode_, **kwargs))
        
        self.dim_in = dim_in
        self.dim_out = dim_out
        self.module_layer = Sequential(*module_layer)
        
        if self.mode=="jk_max" :
            self.jknet = JumpingKnowledge("max")
        elif self.mode=="jk_cat" :
            self.jknet = JumpingKnowledge("cat")
        else : 
            pass
        
        if self.mode=="jk_cat" :
            dim_fin = sum(dim[1:])
        elif self.mode=="dense" :
            dim_fin = dim[1] * len(dim[1:])
        elif self.mode in ["res", "res+"] :
            dim_fin = dim[1]
        else :
            dim_fin = dim[-1]
        self.dim_fin = dim_fin
        
    def forward(self, x, edge_index=None, edge_attr=None, edge_type=None, **kwargs) :
        x_list = list()
        
        for i, layer in enumerate(self.module_layer) :
            x = layer(x, edge_index, edge_attr=edge_attr, edge_type=edge_type)
            
            if self.mode in ["jk_max", "jk_cat"] :
                x_list.append(x)
        
        if self.mode in ["jk_max", "jk_cat"] :
            x = self.jknet(x_list)
        
        return x


class MLP(Module) :
    
    def __init__(self, dim=[128, 128, 128], act=ReLU(), 
                 drop=0.0, norm=None, last_nn=False, norm_3d=False) :
        super().__init__()
        
        self.dim = dim
        self.dim_fin = dim[-1]
        self.norm_3d = norm_3d
        self.layers = ModuleList()
        # norm_3d is needed when treating 3D Tensor 
        # with BatchNorm1d normalizing features of last dimesion
        
        for i in range(1, len(dim)) :
            self.layers.append(Linear(dim[i-1], dim[i]))
            
            if i==len(dim)-1 and last_nn :
                pass
            else :
                if norm not in [None, False, "None", "False"] :
                    self.layers.append(norm_layer(norm, dim[i]))
                if act not in [None, False, "None", "False"] :
                    self.layers.append(act)
                if drop not in [None, False, "None", "False"] and drop > 0 :
                    self.layers.append(Dropout(p=drop))
            
    def forward(self, x) :
        if self.norm_3d :
            num_element = x.shape[1]
            x = x.view(-1, self.dim[0])
            # [batch*num_element, dim[0]]
        
        for layer in self.layers :
            x = layer(x)
            # [batch, dim[-1]] or [batch*num_element, dim[-1]]
        
        if self.norm_3d :
            x = x.view(-1, num_element, self.dim[-1])
            # [batch, num_element, dim[-1]]
        
        return x


# class MLP(Sequential) :
# 
#     def __init__(self, dim=[128, 128, 128], act=ReLU(), drop=0.0, norm=None, last_nn=False) :
#         mlp = []
#         for i in range(1, len(dim)) :
#             mlp.append(Linear(dim[i-1], dim[i]))
# 
#             if i==len(dim)-1 and last_nn :
#                 pass
#             else :
#                 if norm not in [None, False, "None", "False"] :
#                     mlp.append(norm_layer(norm, dim[i]))
#                 if act not in [None, False, "None", "False"] :
#                     mlp.append(act)
#                 if drop not in [None, False, "None", "False"] and drop > 0 :
#                     mlp.append(Dropout(p=drop))
# 
#         self.mlp = mlp
#         super(MLP, self).__init__(*self.mlp)

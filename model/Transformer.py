#!/usr/bin/env python

import copy
import torch
from torch.nn import ReLU, Tanh, GELU, Softmax, Parameter
from torch.nn import Linear, Dropout, LayerNorm, Module, ModuleList, MultiheadAttention

from utils.utils_gnn import *


class DualAttentionLayer(Module) :
    # https://doi.org/10.48550/arXiv.1606.00061
    # [2016] Hierarchical Question-Image Co-Attention for Visual Question Answering
    
    def __init__(self, dim_cell=16, dim_drug=64, dim_attn=16, 
                 scale=True, add_norm=True, drop=0.1, coef_ffnn=0) :
        super().__init__()
        self.act = Tanh()
        self.add_norm = add_norm
        self.add_ffnn = coef_ffnn!=0
        self.act_attn = Softmax(dim=1)
        self.scale = dim_attn**0.5 if scale else 1
        
        self.drop = Dropout(p=drop)
        self.proj_feat = Linear(dim_cell, dim_drug)
        self.proj_path = Linear(dim_cell, dim_attn)
        self.proj_subs = Linear(dim_drug, dim_attn)
        self.attn_path = Linear(dim_attn, 1)
        self.attn_subs = Linear(dim_attn, 1)
        
        if self.add_norm :
            self.norm_cell = LayerNorm(dim_cell)
            self.norm_drug = LayerNorm(dim_drug)
            
        if self.add_ffnn :
            dim_cell_ffnn = int(coef_ffnn*dim_cell)
            dim_drug_ffnn = int(coef_ffnn*dim_drug)
            self.norm_cell_ffnn = LayerNorm(dim_cell)
            self.norm_drug_ffnn = LayerNorm(dim_drug)
            self.ffnn_cell = FeedForward(dim_in=dim_cell, dim_ffnn=dim_cell_ffnn, drop=drop, act_gelu=False)
            self.ffnn_drug = FeedForward(dim_in=dim_drug, dim_ffnn=dim_drug_ffnn, drop=drop, act_gelu=False)
        
    def forward(self, cell, drug) :
        # cell : [batch, n_path, dim_cell]
        # drug : [batch, n_subs, dim_drug]
        
        # Feature Affinity Matrix 
        # P = Tanh(C*Wb*D^T)
        proj_feat = self.proj_feat(cell)
        # [batch, n_path, dim_drug] < [batch, n_path, dim_cell] x [dim_cell, dim_drug]
        proj_feat = torch.matmul(proj_feat, drug.transpose(2, 1))
        # [batch, n_path, n_subs] < [batch, n_path, dim_drug] x [batch, dim_drug, n_subs]
        proj_feat = self.act(proj_feat)
        # [batch, n_path, n_subs]
        
        # Project Cell to Drug Latent Space 
        # Ca = Tanh(C*Wc + P*D*Wd)
        proj_cell1 = self.proj_path(cell)
        # [batch, n_path, dim_attn] < [batch, n_path, dim_cell] x [dim_cell, dim_attn]
        proj_cell2 = torch.matmul(proj_feat, drug)
        # [batch, n_path, dim_drug] < [batch, n_path, n_subs] x [batch, n_subs, dim_drug]
        proj_cell2 = self.proj_subs(proj_cell2)
        # [batch, n_path, dim_attn] < [batch, n_path, dim_drug] x [dim_drug, dim_attn]
        proj_cell = self.act(proj_cell1 + proj_cell2)
        # [batch, n_path, dim_attn]
        
        # Project Drug to Cell Latent Space
        # Da = Tanh(D*Wd + P^T*C*Wc)
        proj_drug1 = self.proj_subs(drug)
        # [batch, n_subs, dim_attn] < [batch, n_subs, dim_drug] x [dim_drug, dim_attn]
        proj_drug2 = torch.matmul(proj_feat.transpose(2, 1), cell)
        # [batch, n_subs, dim_cell] < [batch, n_subs, n_path] x [batch, n_path, dim_cell]
        proj_drug2 = self.proj_path(proj_drug2)
        # [batch, n_subs, dim_attn] < [batch, n_subs, dim_cell] x [dim_cell, dim_attn]
        proj_drug = self.act(proj_drug1 + proj_drug2)
        # [batch, n_subs, dim_attn]
        
        # Attention Score of Pathway & Substructure
        # Ac = wc * Ca / dim_attn^1/2
        # Ad = wd * Da / dim_attn^1/2
        attn_path = self.attn_path(proj_cell)
        # [batch, n_path, 1] < [batch, n_path, dim_attn] x [dim_attn, 1]
        attn_subs = self.attn_subs(proj_drug)
        # [batch, n_subs, 1] < [batch, n_subs, dim_attn] x [dim_attn, 1]
        attn_path = self.act_attn(attn_path / self.scale)
        # [batch, n_path, 1]
        attn_subs = self.act_attn(attn_subs / self.scale)
        # [batch, n_subs, 1]
        
        # Apply Attention Score
        # C' = Ac (dot) C
        # D' = Ad (dot) D
        cell_after = torch.mul(attn_path, cell)
        # [batch, n_path, dim_cell] < [batch, n_path, 1] dot [batch, n_path, dim_cell]
        drug_after = torch.mul(attn_subs, drug)
        # [batch, n_subs, dim_drug] < [batch, n_subs, 1] dot [batch, n_subs, dim_drug]
        
        if self.add_norm :
            cell_after = self.norm_cell(cell + self.drop(cell_after))
            # [batch, n_path, dim_cell]
            drug_after = self.norm_drug(drug + self.drop(drug_after))
            # [batch, n_subs, dim_drug]
        
        if self.add_ffnn :
            cell_after_ = self.ffnn_cell(cell_after)
            # [batch, n_path, dim_cell]
            drug_after_ = self.ffnn_drug(drug_after)
            # [batch, n_subs, dim_drug]
            cell_after = self.norm_cell_ffnn(cell_after + self.drop(cell_after_))
            # [batch, n_path, dim_cell]
            drug_after = self.norm_drug_ffnn(drug_after + self.drop(drug_after_))
            # [batch, n_subs, dim_drug]
        
        attn_path = attn_path.squeeze(-1)
        # [batch, n_path]
        attn_subs = attn_subs.squeeze(-1)
        # [batch, n_subs]
        
        return cell_after, drug_after, attn_path, attn_subs


class DualAttention(Module) :

    def __init__(self, attention_layer, n_layer=2) :
        super().__init__()
        self.layers = ModuleList([copy.deepcopy(attention_layer) for _ in range(n_layer)])

    def forward(self, cell, drug) :
        for layer in self.layers :
            cell, drug, attn_path, attn_subs = layer(cell, drug)
        return cell, drug, attn_path, attn_subs


# class DualEncoderLayer(Module) :
# 
#     def __init__(self, n_heads=4, dim_cell=48, dim_drug=128, n_enc_layer=1,
#                  n_path=192, n_subs=3, coef_ffnn=2, drop=0.1, act_gelu=False) :
#         super().__init__()
#         self.act = Tanh()
#         self.n_path = n_path
#         self.n_subs = n_subs
#         # self.dim_cell = dim_cell
#         # self.dim_drug = dim_drug
# 
#         self.dim_cell_attn = 32
#         self.dim_drug_attn = dim_drug
#         dim_ffnn_cell = int(coef_ffnn * self.dim_cell_attn)
#         dim_ffnn_drug = int(coef_ffnn * self.dim_drug_attn)
# 
#         self.proj_cell = Linear(dim_cell, self.dim_cell_attn)
#         self.proj_drug = Linear(dim_drug, self.dim_drug_attn)
#         self.proj_path = Linear(self.dim_cell_attn*n_path, self.dim_drug_attn)
#         self.proj_subs = Linear(self.dim_drug_attn*n_subs, self.dim_cell_attn)
# 
#         enc_layer_cell = EncoderLayer(n_heads=n_heads, dim_in=self.dim_cell_attn, dim_ffnn=dim_ffnn_cell, drop=drop, act_gelu=act_gelu)
#         enc_layer_drug = EncoderLayer(n_heads=n_heads, dim_in=self.dim_drug_attn, dim_ffnn=dim_ffnn_drug, drop=drop, act_gelu=act_gelu)
#         self.enc_cell = Encoder(enc_layer_cell, n_layer=n_enc_layer)
#         self.enc_drug = Encoder(enc_layer_drug, n_layer=n_enc_layer)
# 
#     def forward(self, cell, drug) :
#         # cell : [batch, n_path, dim_cell]
#         # drug : [batch, n_subs, dim_drug]
# 
#         cell = self.proj_cell(cell)
#         # [batch, n_path, dim_cell_attn]
#         drug = self.proj_drug(drug)
#         # [batch, n_subs, dim_drug_attn]
#         cell_global = cell.view(-1, self.n_path*self.dim_cell_attn)
#         # [batch, n_path*dim_cell_attn]
#         drug_global = drug.view(-1, self.n_subs*self.dim_drug_attn)
#         # [batch, n_subs*dim_drug_attn]
# 
#         cell_global = self.proj_path(cell_global)
#         # [batch, dim_drug_attn]
#         drug_global = self.proj_subs(drug_global)
#         # [batch, dim_cell_attn]
#         cell_global = self.act(cell_global.unsqueeze(-2))
#         # [batch, 1, dim_drug_attn]
#         drug_global = self.act(drug_global.unsqueeze(-2))
#         # [batch, 1, dim_cell_attn]
# 
#         cell_attn = torch.cat([cell, drug_global], axis=-2)
#         # [batch, n_path+1, dim_cell_attn]
#         drug_attn = torch.cat([drug, cell_global], axis=-2)
#         # [batch, n_subs+1, dim_drug_attn]
#         cell_attn, attn_path = self.enc_cell(cell_attn)
#         # [batch, n_path+1, dim_cell_attn] & [batch, n_path+1, n_path+1]
#         drug_attn, attn_subs = self.enc_drug(drug_attn)
#         # [batch, n_subs+1, dim_drug_attn] & [batch, n_subs+1, n_subs+1]
# 
#         # cell_attn = cell_attn[:, :self.n_path, :]
#         # # [batch, n_path, dim_cell_attn]
#         # drug_attn = drug_attn[:, :self.n_subs, :]
#         # # [batch, n_subs, dim_drug_attn]
# 
#         return cell_attn, drug_attn, attn_path, attn_subs
# 
# 
# class DualEncoder(Module) :
# 
#     def __init__(self, attention_layer, n_layer=2) :
#         super().__init__()
#         self.dim_cell_attn = attention_layer.dim_cell_attn
#         self.dim_drug_attn = attention_layer.dim_drug_attn
#         self.layers = ModuleList([copy.deepcopy(attention_layer) for _ in range(n_layer)])
# 
#     def forward(self, cell, drug) :
#         for layer in self.layers :
#             cell, drug, attn_path, attn_subs = layer(cell, drug)
#         return cell, drug, attn_path, attn_subs


class FeedForward(Module) :

    def __init__ (self, dim_in=128, dim_ffnn=256, drop=0.1, act_gelu=False) :
        super().__init__()
        self.dropout = Dropout(p=drop)
        self.fc1 = Linear(dim_in, dim_ffnn)
        self.fc2 = Linear(dim_ffnn, dim_in)
        self.act = ReLU() if not act_gelu else GELU()

    def forward(self, x) :
        out = self.fc1(x)
        out = self.act(out)
        out = self.dropout(out)
        out = self.fc2(out)
        return out


class EncoderLayer(Module) :

    def __init__(self, n_heads=4, dim_in=128, dim_ffnn=256, drop=0.1, act_gelu=False) :
        super().__init__()
        self.dropout = Dropout(p=drop)
        self.layerNorm1 = LayerNorm(dim_in)
        self.layerNorm2 = LayerNorm(dim_in)
        self.ffnn = FeedForward(dim_in, dim_ffnn, drop, act_gelu)
        self.self_attn = MultiheadAttention(embed_dim=dim_in, num_heads=n_heads, dropout=drop, batch_first=True)

    def forward(self, source) :
        out, attn = self.self_attn(query=source, key=source, value=source, attn_mask=None)
        out = self.layerNorm1(source + self.dropout(out))
        out_ = self.ffnn(out)
        out = self.layerNorm2(out + self.dropout(out_))
        return out, attn


class Encoder(Module) :

    def __init__(self, encoder_layer, n_layer=2) :
        super().__init__()
        self.layers = ModuleList([copy.deepcopy(encoder_layer) for _ in range(n_layer)])

    def forward(self, source) :
        for layer in self.layers :
            source, attn = layer(source)
        return source, attn


# class DecoderLayer(Module) :
# 
#     def __init__(self, n_heads=4, dim_in=128, dim_ffnn=256, drop=0.1, act_gelu=False) :
#         super().__init__()
#         self.dropout = Dropout(p=drop)
#         self.layerNorm1 = LayerNorm(dim_in)
#         self.layerNorm2 = LayerNorm(dim_in)
#         self.layerNorm3 = LayerNorm(dim_in)
#         self.ffnn = FeedForward(dim_in, dim_ffnn, drop, act_gelu)
#         self.self_attn = MultiheadAttention(embed_dim=dim_in, num_heads=n_heads, dropout=drop, batch_first=True)
#         self.cross_attn = MultiheadAttention(embed_dim=dim_in, num_heads=n_heads, dropout=drop, batch_first=True)
#         
#     def forward(self, source, target) :
#         out, _ = self.self_attn(query=target, key=target, value=target, attn_mask=None)
#         out = self.layerNorm1(target + self.dropout(out))
#         out_, attn = self.cross_attn(query=target, key=source, value=source, attn_mask=None)
#         out = self.layerNorm2(target + self.dropout(out_))
#         out_ = self.ffnn(out)
#         out = self.layerNorm3(out + self.dropout(out_))
#         return out, attn
# 
# 
# class Decoder(Module) :
#   
#     def __init__(self, decoder_layer, n_layers=2) :
#         super().__init__()
#         self.layers = ModuleList([copy.deepcopy(decoder_layer) for _ in range(n_layers)])
#     
#     def forward(self, source, target) :
#         for layer in self.layers :
#             target, attn = layer(source, target)
#         return target, attn
# 
# 
# class Transformer(Module) :
#     def __init__(self, Encoder, Decoder, n_enc=2, n_dec=2) :
#         super().__init__()
#         self.encoder = Encoder
#         self.decoder = Decoder
#     
#     def forward(self, source, target) :
#         source, self_attn = self.encoder(source)
#         target, cross_attn = self.encoder(source, target)
#         return source, target, self_attn, cross_attn
# 
# 
# if __name__ == '__main__' :
#     x = torch.randn(16, 4, 128)
#     encoder_layer = EncoderLayer(act_gelu=True)
#     encoder = Encoder(encoder_layer)
# 
#     a, b = encoder_layer(x)
#     print("# Encoder_Layer [In] : {}".format(x.shape))
#     print("# Encoder_Layer [Out] : {}".format(a.shape))
#     print("# Encoder_Layer [Attn] : {}".format(b.shape))
# 
#     c, d = encoder(x)
#     print("# Encoder [In] : {}".format(x.shape))
#     print("# Encoder [Out] : {}".format(c.shape))
#     print("# Encoder [Attn] : {}".format(d.shape))

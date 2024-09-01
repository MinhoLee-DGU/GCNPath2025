import math
import datetime

import torch
from torch import nn
import torch.nn.functional as F

class CoA_SSI_DDI(nn.Module):
    def __init__(self, dim_attn=64):
        super().__init__()
        # dim_attn_ = dim_attn
        dim_attn_ = dim_attn//2
        self.w_q = nn.Parameter(torch.zeros(dim_attn, dim_attn_))
        self.w_k = nn.Parameter(torch.zeros(dim_attn, dim_attn_))
        self.bias = nn.Parameter(torch.zeros(dim_attn_))
        self.a = nn.Parameter(torch.zeros(dim_attn_))

        nn.init.xavier_uniform_(self.w_q)
        nn.init.xavier_uniform_(self.w_k)
        nn.init.xavier_uniform_(self.bias.view(*self.bias.shape, -1))
        nn.init.xavier_uniform_(self.a.view(*self.a.shape, -1))
    
    def forward(self, cell, drug):
        cell_attn = cell @ self.w_k
        # [batch, n_path, dim_attn/2] < [batch, n_path, dim_attn] x [dim_attn, dim_attn/2]
        drug_attn = drug @ self.w_q
        # [batch, n_subs, dim_attn/2] < [batch, n_subs, dim_attn] x [dim_attn, dim_attn/2]
        e_activations = cell_attn.unsqueeze(-2) + drug_attn.unsqueeze(-3) + self.bias
        # [batch, n_path, n_subs, dim_attn/2]
        attentions = torch.tanh(e_activations) @ self.a
        # [batch, n_path, n_subs]
        
        return attentions


class Final_SSI_DDI(nn.Module) :
    def __init__(self) :
        super().__init__()
    
    def forward(self, cell, drug, attentions) :
        # breakpoint()
        cell = F.normalize(cell, dim=-1)
        # [batch, n_path, dim_attn]
        drug = F.normalize(drug, dim=-1)
        # [batch, n_subs, dim_attn]
        concat = torch.matmul(cell, drug.transpose(-2, -1))
        # [batch, n_path, n_subs]
        concat = torch.mul(attentions, concat)
        # [batch, n_path, n_subs]
        # concat = concat.sum(dim=(-2, -1))
        # # [batch]
        
        return concat

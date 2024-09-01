#!/usr/bin/env python

import os
import math
import numpy as np
import pandas as pd

import keras
from keras.utils import Sequence


def morgan_fp(mol, nBits=256, radius=2) :
    from rdkit.Chem import AllChem
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    return np.array(fp)


def process_drug(smiles, drugs=None, nBits=512, radius=2, return_df=True) :
  
    drug_dict = {}
    from rdkit import Chem
    if drugs is None : drugs = ['drug_{}'.format(i) for i in range(len(smiles))]
    
    for drug, smi in zip(drugs, smiles) :
        mol = Chem.MolFromSmiles(smi)
        if mol is None : 
            print("Drug {} is not processed...".format(drug))
            pass
        else :
            drug_data_tp = morgan_fp(mol, nBits=nBits, radius=radius)
            drug_dict[drug] = drug_data_tp
    
    if return_df :
        drug_name = list(drug_dict.keys())
        drug_dict = np.array(list(drug_dict.values()))
        drug_dict = pd.DataFrame(drug_dict, index=drug_name)
    
    return drug_dict


class DataLoader(Sequence) :
    
    def __init__(self, ic50_data, data, args, batch_size=1024, shuffle=True) :
        
        self.data = data        
        self.ic50_data = ic50_data
        self.batch_size = batch_size
        
        self.shuffle = shuffle
        self.col_cell = args.col_cell
        self.col_drug = args.col_drug
        self.col_ic50 = args.col_ic50
        self.on_epoch_end()

    def __len__(self) :
        return int(np.ceil(len(self.ic50_data) / self.batch_size))

    def __getitem__(self, idx) :
        idx_start = idx*self.batch_size
        idx_end = (idx+1)*self.batch_size
        # if idx_end>len(self.indices) : idx_end = len(self.indices)
      
        indices = self.indices[idx_start:idx_end]
        cells = self.ic50_data[self.col_cell].iloc[indices]
        drugs = self.ic50_data[self.col_drug].iloc[indices]
        
        Batch_C = [self.data[i].loc[cells].values for i in range(186)]
        Batch_D = self.data[186].loc[drugs].values
        batch_y = self.ic50_data[self.col_ic50].iloc[indices]
        batch_y = np.array(batch_y.astype("float32"))
        return [*Batch_C, Batch_D], batch_y
        
    def on_epoch_end(self) :
        self.indices = np.arange(len(self.ic50_data))
        if self.shuffle :
            np.random.shuffle(self.indices)


def filt_ic50_(ic50_data, data, args) :
    
    ic50_data.dropna(subset=[args.col_cell, args.col_drug], inplace=True)
    if args.col_cell in ["COSMIC_ID"] :
        ic50_data = ic50_data.astype({args.col_cell:"int"})
    ic50_data = ic50_data.astype({args.col_cell:"str"})
    ic50_data = ic50_data.astype({args.col_drug:"str"})
    
    # cell_list = list(cell_data.keys())
    # drug_list = list(drug_data.keys())
    cell_list = list(data[0].index)
    drug_list = list(data[-1].index)

    cell_list = list(set(ic50_data[args.col_cell]) & set(cell_list))
    drug_list = list(set(ic50_data[args.col_drug]) & set(drug_list))
    ic50_data = ic50_data[ic50_data[args.col_cell].isin(cell_list)]
    ic50_data = ic50_data[ic50_data[args.col_drug].isin(drug_list)]

    print("# [Cell] {}".format(len(cell_list)))
    print("# [Drug] {}".format(len(drug_list)))
    print("# [IC50] {}".format(ic50_data.shape[0]))

    return ic50_data


def index_to_str(data) :
    for i in range(len(data)) :
        data[i].index = data[i].index.astype(str)
    return(data)


if __name__ == '__main__' :
    # dir_drug = "_data/SMILES_GDSC.csv"
    # dir_res = "_data/ProcessedFile/drug.csv"
    # 
    # drug_data = pd.read_csv(dir_drug)
    # drug_data = process_drug(drug_data["SMILES_CAN"], drug_data["Drug_CID"])
    # drug_data.index.name = "Drug name"
    # drug_data.to_csv(dir_res)
     
    dir_drug = "_data/SMILES_ChEMBL.csv"
    dir_res = "_data/ProcessedFile/drug_chembl.csv"
    
    drug_data = pd.read_csv(dir_drug)
    drug_data = process_drug(drug_data["Canonical_SMILES"], drug_data["Molecule_ChEMBL_ID"])
    drug_data.index.name = "Drug name"
    drug_data.to_csv(dir_res)

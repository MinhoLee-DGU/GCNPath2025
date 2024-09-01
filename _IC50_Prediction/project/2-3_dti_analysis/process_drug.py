import os
import re
import pickle
import numpy as np
import pandas as pd

import torch


def morgan_fp(mol, nBits=256, radius=2) :
    from rdkit.Chem import AllChem
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    return np.array(fp)


def pubchem_fp(mol) :
    from deepchem import feat
    featurizer = feat.PubChemFingerprint()
    fp = featurizer.featurize(mol)
    return np.array(fp[0])


def process_drug(smiles, drugs=None, nBits=256, radius=2) :
  
    drug_dict = {}
    from rdkit import Chem
    from torch_geometric.data import Data
    if drugs is None : drugs = ['drug_{}'.format(i) for i in range(len(smiles))]
    
    for drug, smi in zip(drugs, smiles) :
        mol = Chem.MolFromSmiles(smi)
        if mol is None : 
            print("Drug {} is not processed...".format(drug))
            pass
        else :
            mol_feats = morgan_fp(mol, nBits=nBits, radius=radius)
            # mol_feats = pubchem_fp(mol)
            
            drug_data_tp = torch.Tensor(mol_feats)
            drug_dict[drug] = drug_data_tp
    
    return drug_dict


def parse_drug():
    import argparse
    parser = argparse.ArgumentParser()
    smi_ = "../../processed_data/drug_data/GDSC/SMILES_GDSC.csv"
    out_dir_="../../processed_data/drug_data/GDSC/GDSC_Drug_Morgan.pickle"

    parser.add_argument("-smi", type=str, default=smi_)
    parser.add_argument("-out_dir", type=str, default=out_dir_)
    # parser.add_argument("-drug_feat", type=int, default=0)
    parser.add_argument("-col_names", type=str, default="Drug_CID")
    parser.add_argument("-col_smiles", type=str, default="SMILES_CAN")
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_drug()
    input_smi = args.smi
    output_dir = args.out_dir
    # drug_feat = args.drug_feat
    col_names = args.col_names
    col_smiles = args.col_smiles
    
    sep = input_smi.split('/')[(-1)].split('.')[(-1)]
    sep = ',' if sep == 'csv' else '\t'
    drug_data = pd.read_csv(input_smi, sep=sep)
    
    if col_names is not None:
        drugs = drug_data[col_names].astype(str)
    else:
        drugs = None
    
    smiles = drug_data[col_smiles]
    drug_data = process_drug(smiles, drugs)
    # drug_data = process_drug(smiles, drugs, drug_feat=drug_feat)
    
    with open(output_dir, 'wb') as (f):
        pickle.dump(drug_data, f)
    
    # if drug_feat==0 :
    fp_value = [_.numpy().astype(int) for _ in drug_data.values()]
    fp_value = np.stack(fp_value)
    
    index = list(drug_data.keys())
    columns = list("FP{}".format(i) for i in range(fp_value.shape[1]))
    drug_data = pd.DataFrame(fp_value, index=index, columns=columns)
    
    output_dir = re.sub(".pickle", ".csv", output_dir)
    drug_data.to_csv(output_dir)
    
    print("### Total drugs {}".format(len(drug_data)))   
    print("### Drug processing succeeded!")


import os
import re
import pickle
import numpy as np
import pandas as pd

import torch
from dgllife.utils import *
from functools import partial


def morgan_fp(mol, nBits=256, radius=2) :
    from rdkit.Chem import AllChem
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    return np.array(fp)


def pubchem_fp(mol) :
    from deepchem import feat
    featurizer = feat.PubChemFingerprint()
    fp = featurizer.featurize(mol)
    return np.array(fp[0])


def get_edge_index(mol) :
    bond_f, bond_r = [], []
    for bond in mol.GetBonds():
        atom_from, atom_to = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_f = bond_f + [atom_from, atom_to]
        bond_r = bond_r + [atom_to, atom_from]
    else:
        edge_index = np.asarray([bond_f, bond_r], dtype=int)
        return edge_index


def process_drug(smiles, drugs=None, drug_feat=3, use_chirality=True, 
                 nBits=256, radius=2, coord_3d=False, coord_seed=2021) :
  
    drug_dict = {}
    from rdkit import Chem
    from torch_geometric.data import Data
    drugs = drugs if drugs is not None else smiles
    
    for drug, smi in zip(drugs, smiles) :
        mol = Chem.MolFromSmiles(smi)
        if mol is None : 
            print("Drug {} is not processed...".format(drug))
            pass
        else :
            try :
                if drug_feat==0 :
                    mol_feats = morgan_fp(mol, nBits=nBits, radius=radius)
                    # mol_feats = pubchem_fp(mol)
                    
                elif drug_feat==1 :
                    from deepchem import feat
                    featurizer = feat.ConvMolFeaturizer(use_chirality=use_chirality)
                    mol_feats = featurizer.featurize(mol)
                    node_feats = mol_feats[0].atom_features
                    edge_index = get_edge_index(mol)
                    
                elif drug_feat==2 :
                    from deepchem import feat
                    featurizer = feat.MolGraphConvFeaturizer(use_edges=True, use_chirality=use_chirality)
                    mol_feats = featurizer.featurize(mol)
                    node_feats = mol_feats[0].node_features
                    edge_index = mol_feats[0].edge_index
                    edge_feats = mol_feats[0].edge_features
                    
                elif drug_feat == 3 :
                    node_featurizer = DrugAtomFeat()
                    edge_featurizer = AttentiveFPBondFeaturizer()
                    
                    if coord_3d :
                        from rdkit.Chem import AllChem
                        from dgllife.utils import get_mol_3d_coordinates
                        status = AllChem.EmbedMolecule(mol)
                        # status = AllChem.EmbedMolecule(mol, AllChem.ETKDG(), randomSeed=coord_seed)
                        mol.GetConformer()
                        position = get_mol_3d_coordinates(mol)
                        
                        # conf = mol.GetConformer()
                        # position = conf.GetAtomPosition(atom_idx)

                    mol_feats = mol_to_bigraph(mol, node_featurizer=node_featurizer, edge_featurizer=edge_featurizer)

                    node_feats = mol_feats.ndata['h'].numpy()
                    edge_feats = mol_feats.edata['e'].numpy()
                    edge_index = mol_feats.edges()
                    edge_index = torch.stack([edge_index[0], edge_index[1]], axis=0).numpy()
                
                elif drug_feat==4 :
                    from utils.utils_sa_ddi import generate_drug_data
                    mol_feats = generate_drug_data(mol)
                
                elif drug_feat==5 :
                    from utils.utils_atom_bond import atom_bond_feat
                    drug_data_tp = atom_bond_feat(mol)
                
                else :
                    raise Exception('drug_feat : 0 [Morgan], 1 [ConvMol], 2 [MolGraph], 3 [Custom]')
                  
                if drug_feat not in [0, 4, 5] :
                    node_feats = torch.from_numpy(node_feats).float()
                    edge_index = torch.from_numpy(edge_index).to(torch.int64)
                    
                if drug_feat==0 :
                    drug_data_tp = torch.Tensor(mol_feats)
                elif drug_feat==1 :
                    drug_data_tp = Data(x=node_feats, edge_index=edge_index)
                elif drug_feat==4 :
                    drug_data_tp = mol_feats
                elif drug_feat==5 :
                    pass
                else :
                    edge_feats = torch.from_numpy(edge_feats).float()
                    drug_data_tp = Data(x=node_feats, edge_index=edge_index, edge_attr=edge_feats)
                    if coord_3d : drug_data_tp.pos = position

                drug_dict[drug] = drug_data_tp
            
            except :
                print("# Drug {} is not processed for unknown reason...".format(drug))
                print("# Drug {} : {}".format(drug, smi))
    
    return drug_dict


class DrugAtomFeat(BaseAtomFeaturizer) :

    def __init__(self, atom_data_field='h', SYMBOLS=None) :
        if SYMBOLS is None:
            SYMBOLS = [
             'C', 'N', 'O', 'H', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 
             'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'K', 'Sn', 
             'Ag', 'Co', 'Se', 'Zn', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Mn', 
             'Cr', 'Pt', 'Hg', '*', 'UNK']
             
        atom_feat = ConcatFeaturizer([
          partial(atom_type_one_hot, allowable_set=SYMBOLS,
          encode_unknown=True),
        atom_mass,
        atom_is_in_ring,
        atom_is_aromatic,
        atom_degree_one_hot,
        atom_total_num_H_one_hot,
        atom_hybridization_one_hot,
        atom_chirality_type_one_hot,
        atom_formal_charge_one_hot,
        atom_explicit_valence_one_hot,
        atom_implicit_valence_one_hot,
        atom_num_radical_electrons_one_hot])
        super().__init__(featurizer_funcs = {atom_data_field : atom_feat})


def parse_drug():
    import argparse
    parser = argparse.ArgumentParser()
    smi_ = "data/drug_data/SMILES_GDSC.csv"
    out_dir_="processed/drug_data/GDSC_Drug_Custom.pickle"

    parser.add_argument("-smi", type=str, default=smi_)
    parser.add_argument("-out_dir", type=str, default=out_dir_)
    parser.add_argument("-drug_feat", type=int, default=3)
    parser.add_argument("-coord_3d", type=int, default=0)
    parser.add_argument("-col_names", type=str, default="Drug_CID")
    parser.add_argument("-col_smiles", type=str, default="SMILES_CAN")
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_drug()
    input_smi = args.smi
    output_dir = args.out_dir
    drug_feat = args.drug_feat
    coord_3d = args.coord_3d
    col_names = args.col_names
    col_smiles = args.col_smiles
    
    sep = input_smi.split('/')[(-1)].split('.')[(-1)]
    sep = ',' if sep == 'csv' else '\t'
    coord_3d = coord_3d!=0
    drug_data = pd.read_csv(input_smi, sep=sep)
    
    if col_names is not None:
        drugs = drug_data[col_names].astype(str)
    else:
        drugs = None
    
    smiles = drug_data[col_smiles]
    drug_data = process_drug(smiles, drugs, drug_feat=drug_feat, coord_3d=coord_3d)
    
    with open(output_dir, 'wb') as (f):
        pickle.dump(drug_data, f)
    
    if drug_feat==0 :
        fp_value = [_.numpy().astype(int) for _ in drug_data.values()]
        fp_value = np.stack(fp_value)
        
        index = list(drug_data.keys())
        columns = list("FP{}".format(i) for i in range(fp_value.shape[1]))
        drug_data = pd.DataFrame(fp_value, index=index, columns=columns)
        
        output_dir = re.sub(".pickle", ".csv", output_dir)
        drug_data.to_csv(output_dir)
    
    print("### Total drugs {}".format(len(drug_data)))   
    print("### Drug processing succeeded!")


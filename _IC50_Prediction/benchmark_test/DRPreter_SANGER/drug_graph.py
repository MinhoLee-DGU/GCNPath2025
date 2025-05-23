from rdkit import Chem
import numpy as np
import pandas as pd
import torch
import torch_geometric
from torch_geometric.data import Data
from dgllife.utils import *


def atom_to_feature_vector(atom):
    """
    Converts rdkit atom object to feature list of indices
    :param mol: rdkit atom object
    :return: list
    8 features are canonical, 2 features are from OGB
    """
    featurizer_funcs = ConcatFeaturizer([atom_type_one_hot,
                                         atom_degree_one_hot,
                                         atom_implicit_valence_one_hot,
                                         atom_formal_charge,
                                         atom_num_radical_electrons,
                                         atom_hybridization_one_hot,
                                         atom_is_aromatic,
                                         atom_total_num_H_one_hot,
                                         atom_is_in_ring,
                                         atom_chirality_type_one_hot,
                                         ])
    atom_feature = featurizer_funcs(atom)
    return atom_feature


def bond_to_feature_vector(bond):
    """
    Converts rdkit bond object to feature list of indices
    :param mol: rdkit bond object
    :return: list
    """
    featurizer_funcs = ConcatFeaturizer([bond_type_one_hot,
                                         # bond_is_conjugated,
                                         # bond_is_in_ring,
                                         # bond_stereo_one_hot,
                                         ])
    bond_feature = featurizer_funcs(bond)

    return bond_feature


def smiles2graph(mol):
    """
    Converts SMILES string or rdkit's mol object to graph Data object without remove salt
    :input: SMILES string (str)
    :return: graph object
    """

    if isinstance(mol, Chem.rdchem.Mol):
        pass
    else:
        mol = Chem.MolFromSmiles(mol)

    # atoms
    atom_features_list = []
    for atom in mol.GetAtoms():
        atom_features_list.append(atom_to_feature_vector(atom))
    x = np.array(atom_features_list, dtype=np.int64)

    # bonds
#     num_bond_features = 3  # bond type, bond stereo, is_conjugated
    num_bond_features = 1  # bond type
    if len(mol.GetBonds()) > 0:  # mol has bonds
        edges_list = []
        edge_features_list = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            edge_feature = bond_to_feature_vector(bond)

            # add edges in both directions
            edges_list.append((i, j))
            edge_features_list.append(edge_feature)
            edges_list.append((j, i))
            edge_features_list.append(edge_feature)

        # data.edge_index: Graph connectivity in COO format with shape [2, num_edges]
        edge_index = np.array(edges_list, dtype=np.int64).T

        # data.edge_attr: Edge feature matrix with shape [num_edges, num_edge_features]
        edge_attr = np.array(edge_features_list, dtype=np.int64)

    else:  # mol has no bonds
        edge_index = np.empty((2, 0), dtype=np.int64)
        edge_attr = np.empty((0, num_bond_features), dtype=np.int64)

    graph = Data(x=torch.tensor(x, dtype=torch.float),
                 edge_index=torch.tensor(edge_index, dtype=torch.long),
                 edge_attr=torch.tensor(edge_attr), dtype=torch.float)

    return graph


def save_drug_graph(file_in, file_out, col_drug="Drug_CID", col_smi="SMILES_CAN"):
    smiles = pd.read_csv(file_in)
    # smiles = pd.read_csv('Data/Drug/drug_smiles.csv')
    
    drug_dict = {}
    for i in range(len(smiles)):
        drug_dict[smiles[col_drug][i]] = smiles2graph(smiles[col_smi][i])
        # drug_dict[smiles.iloc[i, 0]] = smiles2graph(smiles.iloc[i, 2])
    
    np.save(file_out, drug_dict)
    # np.save('Data/Drug/drug_feature_graph.npy', drug_dict)

    return drug_dict


if __name__ == '__main__':
    # graph = smiles2graph('O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5')
    # print(graph.x.shape)
    # print(graph.edge_attr.shape)
    
    dir_drug = "_data/SMILES_GDSC.csv"
    dir_drug_ = "_data/drug_feature_graph.npy"
    save_drug_graph(dir_drug, dir_drug_, col_drug="Drug_CID", col_smi="SMILES_CAN")

    dir_drug = "_data/SMILES_ChEMBL.csv"
    dir_drug_ = "_data/drug_chembl.npy"
    save_drug_graph(dir_drug, dir_drug_, col_drug="Molecule_ChEMBL_ID", col_smi="Canonical_SMILES")


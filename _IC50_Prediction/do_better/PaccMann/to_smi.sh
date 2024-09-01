#!/usr/bin/bash

# col_names=ChEMBL_ID
# col_smiles=Smiles
# input_dir=_data/SMILES_ChEMBL_COAD.csv
# output_dir=_data/SMILES_ChEMBL_COAD.smi
# python to_smi.py -input_dir $input_dir -output_dir $output_dir -col_names $col_names -col_smiles $col_smiles

col_names=Molecule_ChEMBL_ID
col_smiles=Canonical_SMILES
input_dir=_data/SMILES_ChEMBL.csv
output_dir=_data/SMILES_ChEMBL.smi
python to_smi.py -input_dir $input_dir -output_dir $output_dir -col_names $col_names -col_smiles $col_smiles

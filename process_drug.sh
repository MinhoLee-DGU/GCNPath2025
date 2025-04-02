#!/usr/bin/bash

feat=3
export PYTHONNOUSERSITE=1

col_names=Drug_CID
col_smi=SMILES_CAN
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
smi=data/drug_data/SMILES_GDSC.csv
out=processed/drug_data/GDSC_Drug_Graph.pickle
python process_drug.py -smi $smi -out_dir $out $option

col_names=Drug_CID
col_smi=SMILES_CAN
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
smi=data/drug_data/SMILES_CCLE.csv
out=processed/drug_data/CCLE_Drug_Graph.pickle
python process_drug.py -smi $smi -out_dir $out $option

col_names=Drug_CID
col_smi=SMILES
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
smi=data/drug_data/TCGA_Drug_Info.csv
out=processed/drug_data/TCGA_Drug_Graph.pickle
python process_drug.py -smi $smi -out_dir $out $option

col_names=Molecule_ChEMBL_ID
col_smi=Canonical_SMILES
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
smi=data/drug_data/SMILES_ChEMBL.csv
out=processed/drug_data/ChEMBL_Drug_Graph.pickle
python process_drug.py -smi $smi -out_dir $out $option


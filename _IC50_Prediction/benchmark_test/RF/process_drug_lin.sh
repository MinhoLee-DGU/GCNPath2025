#!/usr/bin/bash
export PYTHONNOUSERSITE=1

feat=0
# col_names=Drug_CID
# col_smi=SMILES_CAN
# option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
# 
# smi=_data/SMILES_GDSC.csv
# out=_data/GDSC_Drug_Morgan.pickle
# python process_drug.py -smi $smi -out_dir $out $option
# 
# col_names=Molecule_ChEMBL_ID
# col_smi=Canonical_SMILES
# option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"
# 
# smi=_data/SMILES_ChEMBL.csv
# out=_data/ChEMBL_Drug_Morgan.pickle
# python process_drug.py -smi $smi -out_dir $out $option

col_names=Drug_CID
col_smi=SMILES_Can
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"

smi=_data/TCGA_Drug_Info.csv
out=_data/TCGA_Drug_Morgan.pickle
python process_drug.py -smi $smi -out_dir $out $option
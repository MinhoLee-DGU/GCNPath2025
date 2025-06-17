#!/usr/bin/bash
export PYTHONNOUSERSITE=1

feat=0
col_names=Drug_CID
col_smi=SMILES_CAN
option="-drug_feat $feat -col_names $col_names -col_smi $col_smi"

smi=data/drug_data/SMILES_GDSC.csv
out=processed/drug_data/GDSC_Drug_Morgan.pickle
python process_drug.py -smi $smi -out_dir $out $option


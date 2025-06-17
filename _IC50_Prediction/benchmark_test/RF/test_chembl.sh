#!/usr/bin/bash

node=0
test_type=Normal

# ic50=_data/IC50_ChEMBL.txt
ic50=_data/IC50_ChEMBL_Time.txt
drug=_data/ChEMBL_Drug_Morgan.pickle
dir_cell=_data

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

use_slurm=1
seed_list=$(seq 2021 2030)
dir_model=Results/IC50_GDSC/$test_type

for seed in ${seed_list[@]} 
do
    param=${dir_model}/param_retrain_seed${seed}.joblib
    out_time=${dir_model}/log_time_chembl_seed${seed}.csv
    out_time_=None
    
    col_cell=SANGER_MODEL_ID
    jname=ChEMBL_S${seed}
    out_file=${dir_model}/pred_chembl_seed${seed}.csv
    cell=${dir_cell}/SANGER_RNA.pickle
    bash test_write.sh $ic50 $cell $drug $param $out_file \
        $out_time $col_cell $col_drug $col_ic50 $node $jname $use_slurm
    
    col_cell=SANGER_MODEL_ID
    jname=ChEMBL_CCLE_S${seed}
    out_file=${dir_model}/pred_chembl_ccle_seed${seed}.csv
    cell=${dir_cell}/CCLE_RNA.pickle
    bash test_write.sh $ic50 $cell $drug $param $out_file \
        $out_time_ $col_cell $col_drug $col_ic50 $node $jname $use_slurm
    
    col_cell=SANGER_MODEL_ID
    jname=ChEMBL_GDSC_S${seed}
    out_file=${dir_model}/pred_chembl_gdsc_seed${seed}.csv
    cell=${dir_cell}/GDSC_RNA.pickle
    bash test_write.sh $ic50 $cell $drug $param $out_file \
        $out_time_ $col_cell $col_drug $col_ic50 $node $jname $use_slurm
done

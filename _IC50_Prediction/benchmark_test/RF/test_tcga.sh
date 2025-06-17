#!/usr/bin/bash

node=0
test_type=Normal

cell=_data/TCGA_RNA.pickle
cell_cb=_data/TCGA_RNA_CB.pickle
ic50=_data/TCGA_Response.csv
drug=_data/TCGA_Drug_Morgan.pickle

col_cell=Sample
col_drug=Drug_CID
col_ic50=0

use_slurm=1
seed_list=$(seq 2021 2030)
dir_model=Results/IC50_GDSC/$test_type
out_time=None

for seed in ${seed_list[@]}
do
    param=${dir_model}/param_retrain_seed${seed}.joblib
    
    jname=TCGA_S${seed}
    out_file=${dir_model}/pred_tcga_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $param $out_file \
        $out_time $col_cell $col_drug $col_ic50 $node $jname $use_slurm
    
    jname=TCGA_S${seed}_CB
    out_file=${dir_model}/pred_tcga_seed${seed}_combat.csv
    bash test_write.sh $ic50 $cell_cb $drug $param $out_file \
        $out_time $col_cell $col_drug $col_ic50 $node $jname $use_slurm
done

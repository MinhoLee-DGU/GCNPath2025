#!/usr/bin/bash

gpu=1
node=0
test_type=Normal

ic50=data/ic50_data/IC50_CCLE.txt
cell=processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle
drug=processed/drug_data/CCLE_Drug_Graph.pickle
dir_model=results/IC50_GDSC/$test_type/RGCN

col_cell=Cell_BROAD_ID
col_drug=Drug_CID
col_ic50=LN_IC50

use_slurm=1
seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]}
do
    jname=GCN_${test_type}_CCLE_Seed${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_file=${dir_model}/pred_ccle_seed${seed}.csv
    out_cam=${dir_model}/gcam_ccle_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done


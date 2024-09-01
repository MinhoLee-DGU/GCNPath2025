#!/usr/bin/bash

gpu=1
node=0
test_type=Normal

ic50=data/ic50_data/IC50_GDSC_Rest.txt
cell=processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
drug=processed/drug_data/GDSC_Drug_Custom.pickle
dir_model=results/IC50_GDSC/$test_type/RGCN

col_cell=Cell
col_drug=Drug
col_ic50=0

use_slurm=1
seed_list=$(seq 2021 2030)
for seed in ${seed_list[@]}
do
    jname=GCN_${test_type}_Rest_Seed${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_file=${dir_model}/pred_rest_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done

# fold_list=$(seq 0 9)
# for nth in ${fold_list[@]}
# do
#     jname=GCN_${test_type}_Rest_${nth}
#     param=${dir_model}/param_retrain_${nth}.pt
#     hparam=${dir_model}/hyper_param_${nth}.pickle
#     out_file=${dir_model}/pred_rest_${nth}.csv
#     bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
#         $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
# done

#!/usr/bin/bash

gpu=0
node=0
test_type=Normal
model_type=RGCN

ic50=data/ic50_data/TCGA_Response.csv
cell=processed/cell_data_biocarta/TCGA_RNA_KNN5_STR9_Reg_Corr.pickle
drug=processed/drug_data/TCGA_Drug_Custom.pickle
dir_model=results/IC50_GDSC/$test_type/$model_type

col_cell=Sample
col_drug=Drug_CID
col_ic50=0
use_slurm=1

# fold_list=$(seq 0 9)
# for nth in ${fold_list[@]}
# do
#     jname=GCN_${model_type}_TCGA_${nth}
#     param=${dir_model}/param_retrain_${nth}.pt
#     hparam=${dir_model}/hyper_param_retrain_${nth}.pickle
#     out_file=${dir_model}/pred_tcga_${nth}.csv
#     
#     bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
#         $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
# done

seed_list=$(seq 2021 2030)
for seed in ${seed_list[@]}
do
    jname=GCN_${model_type}_TCGA_Seed${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_file=${dir_model}/pred_tcga_seed${seed}.csv
    
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done


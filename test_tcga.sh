#!/usr/bin/bash

gpu=0
node=0
test_type=Normal
model_type=RGCN

ic50=data/ic50_data/TCGA_Response.csv
cell=processed/cell_data_biocarta/TCGA_RNA_KNN5_STR9_Reg_Corr.pickle
cell_cb=processed/cell_data_biocarta/TCGA_RNA_KNN5_STR9_Reg_Corr_CB.pickle
drug=processed/drug_data/TCGA_Drug_Graph.pickle
dir_model=results/IC50_GDSC/$test_type/$model_type

out_cam=None
col_cell=Sample
col_drug=Drug_CID
col_ic50=0

use_slurm=1
seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]}
do
    # TCGA
    jname=TCGA_S${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_file=${dir_model}/pred_tcga_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm

    # TCGA [ComBat O]
    jname=TCGA_S${seed}_CB
    out_file=${dir_model}/pred_tcga_seed${seed}_combat.csv
    bash test_write.sh $ic50 $cell_cb $drug $param $hparam $out_file \
        $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done


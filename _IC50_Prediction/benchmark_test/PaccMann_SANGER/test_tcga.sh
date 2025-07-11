#!/usr/bin/bash

gpu=0

ic50=_data/TCGA_Response.csv
drug=_data/SMILES_TCGA.smi
dir_in=Results/IC50_GDSC/Normal

col_cell=Sample
col_drug=Drug_CID
col_ic50=0

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    # TCGA
    f_name=TCGA_S${seed}
    cell=_data/EXP_TCGA.csv
    dir_test=$dir_in/pred_tcga_seed${seed}.csv
    dir_param=$dir_in/pred_backup/param_retrain_seed${seed}.pt.tar
    dir_hparam=$dir_in/pred_backup/hyper_param_retrain_seed${seed}.json
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu

    # TCGA [ComBat]
    f_name=TCGA_S${seed}_CB
    cell=_data/EXP_TCGA_ComBat.csv
    dir_test=$dir_in/pred_tcga_seed${seed}_combat.csv
    dir_param=$dir_in/pred_backup/param_retrain_seed${seed}.pt.tar
    dir_hparam=$dir_in/pred_backup/hyper_param_retrain_seed${seed}.json
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

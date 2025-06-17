#!/usr/bin/bash

gpu=0

cell=_data/cell_feature_std_disjoint_tcga.npy
cell_cb=_data/cell_feature_std_disjoint_tcga_combat.npy
drug=_data/drug_tcga.npy
ic50=_data/TCGA_Response.csv
dir_in=Results/IC50_GDSC/Normal/DRPreter

col_cell=Sample
col_drug=Drug_CID
col_ic50=0

seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do
    # TCGA
    # f_name=TCGA_S${seed}
    # dir_test=${dir_in}/pred_tcga_seed${seed}.csv
    # dir_param=${dir_in}/pred_backup/param_retrain_seed${seed}.pth
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu

    # TCGA [ComBat O]
    f_name=TCGA_S${seed}_CB
    dir_test=${dir_in}/pred_tcga_seed${seed}_combat.csv
    dir_param=${dir_in}/pred_backup/param_retrain_seed${seed}.pth
    bash test_write.sh $ic50 $cell_cb $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

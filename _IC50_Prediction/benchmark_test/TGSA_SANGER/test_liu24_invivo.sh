#!/usr/bin/bash

gpu=0

cell=_data_liu24_invivo/cell_feature_all.npy
cell_cb=_data_liu24_invivo/cell_feature_all_combat.npy
drug=_data_liu24_invivo/drug_liu_2024.npy
ic50=Liu_2024_in_vivo.txt
dir_in=Results/IC50_GDSC/Normal/TGDRP

col_cell=Sample
col_drug=Drug_CID
col_ic50=0

seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do
    # ComBat X
    f_name=Liu_S${seed}
    dir_test=${dir_in}/pred_liu24_seed${seed}.csv
    dir_param=${dir_in}/pred_backup/param_retrain_seed${seed}.pth
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu

    # ComBat O
    f_name=Liu_S${seed}_CB
    dir_test=${dir_in}/pred_liu24_seed${seed}_combat.csv
    dir_param=${dir_in}/pred_backup/param_retrain_seed${seed}.pth
    bash test_write.sh $ic50 $cell_cb $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

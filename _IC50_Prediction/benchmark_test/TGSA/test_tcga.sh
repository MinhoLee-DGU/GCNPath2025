#!/usr/bin/bash

gpu=0

cell=_data_tcga/cell_feature_all_cncat.npy
cell_cb=_data_tcga/cell_feature_all_cncat_combat.npy
drug=_data/drug_tcga.npy
ic50=TCGA_Response.csv
dir_in=Results/IC50_GDSC/Normal/TGDRP

col_cell=Sample
col_drug=Drug_CID
col_ic50=0

seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do
    # TCGA
    f_name=TCGA_S${seed}
    dir_test=${dir_in}/pred_tcga_seed${seed}.csv
    dir_param=${dir_in}/param_retrain_seed${seed}.pth
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu

    # TCGA [ComBat O]
    f_name=TCGA_S${seed}_CB
    dir_test=${dir_in}/pred_tcga_seed${seed}_combat.csv
    dir_param=${dir_in}/param_retrain_seed${seed}.pth
    bash test_write.sh $ic50 $cell_cb $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

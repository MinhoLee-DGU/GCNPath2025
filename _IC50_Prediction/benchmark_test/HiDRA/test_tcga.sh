#!/usr/bin/bash

drug=drug_tcga.csv
ic50=TCGA_Response.csv

col_cell=Sample
col_drug=Drug_CID
col_ic50=0
dir_in=Results/IC50_GDSC/Normal

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    dir_param=$dir_in/param_retrain_seed${seed}.h5
    
    # TCGA
    # f_name=TCGA_S${seed}
    # cell=_data/ProcessedFile_TCGA
    # dir_test=$dir_in/pred_tcga_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name

    # TCGA [ComBat O]
    f_name=TCGA_S${seed}_CB
    cell=_data/ProcessedFile_TCGA_ComBat
    dir_test=$dir_in/pred_tcga_seed${seed}_combat.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name
done

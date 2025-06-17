#!/usr/bin/bash

# ic50=IC50_ChEMBL.txt
ic50=IC50_ChEMBL_Time.txt
drug=drug_chembl.csv

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID
dir_in=Results/IC50_GDSC/Normal

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    dir_param=$dir_in/param_retrain_seed${seed}.h5
    
    # GDSC
    col_cell=COSMIC_ID
    f_name=ChEMBL_S${seed}
    cell=_data/ProcessedFile
    dir_test=$dir_in/pred_chembl_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name
    
    # CCLE
    # col_cell=BROAD_ID
    # f_name=ChEMBL_CCLE_S${seed}
    # cell=_data/ProcessedFile_CCLE
    # dir_test=$dir_in/pred_chembl_ccle_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name
    
    # SANGER
    # col_cell=SANGER_MODEL_ID
    # f_name=ChEMBL_SANGER_S${seed}
    # cell=_data/ProcessedFile_SANGER
    # dir_test=$dir_in/pred_chembl_sanger_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name
done

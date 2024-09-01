#!/usr/bin/bash

ic50=IC50_ChEMBL.txt
drug=drug_chembl.csv
dir_in=Results/IC50_GDSC/Normal

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

mode=$(seq 0 9)
for nth in ${mode[@]}
do
    dir_param=$dir_in/param_$nth.h5
    
    # GDSC
    col_cell=COSMIC_ID
    f_name=HiDRA_ChEMBL_$nth
    cell=_data/ProcessedFile
    dir_test=$dir_in/pred_chembl_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_test $col_cell $col_drug $col_ic50 $f_name
    
    # CCLE
    col_cell=BROAD_ID
    f_name=HiDRA_ChEMBL_CCLE_$nth
    cell=_data/ProcessedFile_CCLE
    dir_test=$dir_in/pred_chembl_ccle_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_test $col_cell $col_drug $col_ic50 $f_name
    
    # SANGER
    col_cell=SANGER_MODEL_ID
    f_name=HiDRA_ChEMBL_SANGER_$nth
    cell=_data/ProcessedFile_SANGER
    dir_test=$dir_in/pred_chembl_sanger_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_test $col_cell $col_drug $col_ic50 $f_name
done

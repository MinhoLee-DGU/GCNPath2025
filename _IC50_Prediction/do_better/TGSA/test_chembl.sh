#!/usr/bin/bash

ic50=IC50_ChEMBL.txt
drug=_data/drug_chembl.npy
dir_in=Results/IC50_GDSC/Normal/TGDRP

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

gpu=1
mode=$(seq 0 9)
for nth in ${mode[@]}
do
    dir_param=$dir_in/param_retrain_$nth.pth
    
    # CCLE
    col_cell=BROAD_ID
    f_name=TGDRP_ChEMBL_$nth
    cell=_data/cell_feature_all.npy
    dir_test=$dir_in/pred_chembl_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # SANGER
    col_cell=BROAD_ID
    f_name=TGDRP_ChEMBL_SANGER_$nth
    cell=_data/cell_feature_all.npy
    dir_test=$dir_in/pred_chembl_sanger_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # GDSC
    col_cell=BROAD_ID
    f_name=TGDRP_ChEMBL_GDSC_$nth
    cell=_data_gdsc/cell_feature_all.npy
    dir_test=$dir_in/pred_chembl_gdsc_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

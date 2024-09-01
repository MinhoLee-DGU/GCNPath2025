#!/usr/bin/bash

ic50=IC50_ChEMBL.txt
drug=_data/drug_chembl.npy
dir_in=Results/IC50_GDSC/Normal/DRPreter

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

gpu=1
mode=$(seq 0 9)

for nth in ${mode[@]}
do  
    # SANGER
    col_cell=SANGER_MODEL_ID
    f_name=DRP_ChEMBL_${nth}
    cell=_data/cell_feature_std_disjoint.npy
    
    dir_test=${dir_in}/pred_chembl_${nth}.csv
    dir_param=${dir_in}/param_retrain_${nth}.pth
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # CCLE
    col_cell=BROAD_ID
    f_name=DRP_ChEMBL_CCLE_${nth}
    cell=_data/cell_feature_std_disjoint_ccle.npy
    
    dir_test=${dir_in}/pred_chembl_ccle_${nth}.csv
    dir_param=${dir_in}/param_retrain_${nth}.pth
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # GDSC
    col_cell=COSMIC_ID
    f_name=DRP_ChEMBL_GDSC_${nth}
    cell=_data/cell_feature_std_disjoint_gdsc.npy
    
    dir_test=${dir_in}/pred_chembl_gdsc_${nth}.csv
    dir_param=${dir_in}/param_retrain_${nth}.pth
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

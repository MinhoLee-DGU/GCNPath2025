#!/usr/bin/bash

ic50=_data/IC50_ChEMBL.txt
drug=_data/SMILES_ChEMBL.smi

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID
dir_in=Results/IC50_GDSC/Normal

gpu=1
fold_list=$(seq 0 9)
for nth in ${fold_list[@]}
do  
    dir_param=$dir_in/param_retrain_${nth}.pt.tar
    dir_hparam=$dir_in/hyper_param_retrain_${nth}.json
    
    # SANGER
    col_cell=SANGER_MODEL_ID
    f_name=PM_ChEMBL_${nth}
    cell=_data/EXP.csv
    
    dir_test=$dir_in/pred_chembl_${nth}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
        
    # CCLE
    col_cell=BROAD_ID
    f_name=PM_ChEMBL_CCLE_${nth}
    cell=_data/EXP_CCLE.csv
    
    dir_test=$dir_in/pred_chembl_ccle_${nth}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu

    # GDSC
    col_cell=COSMIC_ID
    f_name=PM_ChEMBL_GDSC_${nth}
    cell=_data/EXP_GDSC.csv
    
    dir_test=$dir_in/pred_chembl_gdsc_${nth}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

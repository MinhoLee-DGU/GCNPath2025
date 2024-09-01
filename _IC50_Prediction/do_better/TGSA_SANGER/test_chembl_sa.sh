#!/usr/bin/bash

ic50=IC50_ChEMBL.txt
drug=_data/drug_chembl.npy

dir_in=Results/IC50_GDSC/Normal/TGDRP
dir_in2=Results/IC50_GDSC/Normal/TGSA

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

mode=$(seq 0 9)
for nth in ${mode[@]}
do
    dir_param=$dir_in/param_retrain_$nth.pth
    dir_tgsa=$dir_in2/param_retrain_$nth.pth

    # SANGER
    dir_sa=_data
    col_cell=SANGER_MODEL_ID
    f_name=TGSA_ChEMBL_$nth
    
    cell=${dir_sa}/cell_feature_all.npy
    dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    dir_test=$dir_in2/pred_chembl_$nth.csv
    bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
        $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
    
    # CCLE
    dir_sa=_data_ccle
    col_cell=SANGER_MODEL_ID
    f_name=TGSA_ChEMBL_CCLE_$nth
    
    cell=${dir_sa}/cell_feature_all.npy
    dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    dir_test=$dir_in2/pred_chembl_ccle_$nth.csv
    bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
        $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
        
    # GDSC
    dir_sa=_data_gdsc
    col_cell=SANGER_MODEL_ID
    f_name=TGSA_ChEMBL_GDSC_$nth
    
    cell=${dir_sa}/cell_feature_all.npy
    dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    dir_test=$dir_in2/pred_chembl_gdsc_$nth.csv
    bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
        $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
done

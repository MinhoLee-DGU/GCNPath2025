#!/usr/bin/bash

# ic50=IC50_ChEMBL.txt
ic50=IC50_ChEMBL_Time.txt
drug=_data/drug_chembl.npy

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

dir_in=Results/IC50_GDSC/Normal/TGDRP
dir_in2=Results/IC50_GDSC/Normal/TGSA

seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do  
    dir_param=$dir_in/param_retrain_seed${seed}.pth
    dir_tgsa=$dir_in2/param_retrain_seed${seed}.pth

    # SANGER
    col_cell=SANGER_MODEL_ID
    f_name=TGSA_ChEMBL_S${seed}
    dir_test=$dir_in2/pred_chembl_seed${seed}.csv
    
    dir_sa=_data
    cell=${dir_sa}/cell_feature_all.npy
    dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
        $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
    
    # CCLE
    # col_cell=SANGER_MODEL_ID
    # f_name=TGSA_ChEMBL_CCLE_S${seed}
    # dir_test=$dir_in2/pred_chembl_ccle_seed${seed}.csv
    
    # dir_sa=_data_ccle
    # cell=${dir_sa}/cell_feature_all.npy
    # dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    # dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    # bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
    #     $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
    
    # GDSC
    # col_cell=SANGER_MODEL_ID
    # f_name=TGSA_ChEMBL_GDSC_S${seed}
    # dir_test=$dir_in2/pred_chembl_gdsc_seed${seed}.csv
    
    # dir_sa=_data_gdsc
    # cell=${dir_sa}/cell_feature_all.npy
    # dir_sim=${dir_sa}/drug_cell_edges_5_knn_gdsc
    # dir_cidx=${dir_sa}/cell_id2idx_dict_gdsc
    
    # bash test_write_sa.sh $ic50 $cell $drug $dir_param $dir_tgsa \
    #     $dir_test $col_cell $col_drug $col_ic50 $dir_sim $dir_cidx $f_name
done

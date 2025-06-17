#!/usr/bin/bash

# ic50=IC50_ChEMBL.txt
ic50=IC50_ChEMBL_Time.txt
drug=_data/drug_chembl.npy
dir_in=Results/IC50_GDSC/Normal/TGDRP

gpu=1
col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do
    dir_param=$dir_in/param_retrain_seed${seed}.pth
    
    # SANGER
    col_cell=SANGER_MODEL_ID
    f_name=TGDRP_ChEMBL_S${seed}
    cell=_data/cell_feature_all.npy
    dir_test=$dir_in/pred_chembl_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # CCLE
    # col_cell=SANGER_MODEL_ID
    # f_name=TGDRP_ChEMBL_CCLE_S${seed}
    # cell=_data_ccle/cell_feature_all.npy
    # dir_test=$dir_in/pred_chembl_ccle_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # GDSC
    # col_cell=SANGER_MODEL_ID
    # f_name=TGDRP_ChEMBL_GDSC_S${seed}
    # cell=_data_gdsc/cell_feature_all.npy
    # dir_test=$dir_in/pred_chembl_gdsc_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

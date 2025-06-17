#!/usr/bin/bash

gpu=1
# ic50=_data/IC50_ChEMBL.txt
ic50=_data/IC50_ChEMBL_Time.txt
drug=_data/SMILES_ChEMBL.smi

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    dir_in=Results/IC50_GDSC/Normal
    dir_param=$dir_in/param_retrain_seed${seed}.pt.tar
    dir_hparam=$dir_in/hyper_param_retrain_seed${seed}.json
    
    # SANGER
    col_cell=SANGER_MODEL_ID
    cell=_data/EXP.csv
    f_name=ChEMBL_S${seed}
    dir_test=$dir_in/pred_chembl_seed${seed}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
        $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # GDSC
    # col_cell=COSMIC_ID
    # cell=_data/EXP_GDSC.csv
    # f_name=ChEMBL_GDSC_S${seed}
    # dir_test=$dir_in/pred_chembl_gdsc_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
    
    # CCLE
    # col_cell=BROAD_ID
    # cell=_data/EXP_CCLE.csv
    # f_name=ChEMBL_CCLE_S${seed}
    # dir_test=$dir_in/pred_chembl_ccle_seed${seed}.csv
    # bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam \
    #     $dir_test $col_cell $col_drug $col_ic50 $f_name $gpu
done

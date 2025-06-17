#!/usr/bin/bash

dix=$1
# ic50=_data/IC50_ChEMBL.txt
ic50=_data/IC50_ChEMBL_Time.txt
cell=cell_mut_matrix.npy
drug=drug_onehot_smiles_chembl_$dix.npy

col_ic50=LN_IC50
col_cell=COSMIC_ID
col_drug=Molecule_ChEMBL_ID
dir_in=Results/IC50_GDSC/Normal

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    f_name=ChEMBL_S${seed}_${dix}
    dir_param=$dir_in/param_retrain_seed${seed}.ckpt
    dir_test=$dir_in/pred_chembl_seed${seed}_${dix}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param \
        $dir_test $col_cell $col_drug $col_ic50 $f_name
done


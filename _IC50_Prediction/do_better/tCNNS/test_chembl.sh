#!/usr/bin/bash

dix=$1
ic50=_data/IC50_ChEMBL.txt
cell=cell_mut_matrix.npy
drug=drug_onehot_smiles_chembl_$dix.npy

col_cell=COSMIC_ID
col_drug=Molecule_ChEMBL_ID
col_ic50=LN_IC50
dir_in=Results/IC50_GDSC/Normal

mode=$(seq 0 9)
for nth in ${mode[@]}
do
    f_name=tCNNS_ChEMBL_${nth}_${dix}
    dir_param=$dir_in/param_retrain_$nth.ckpt
    dir_test=$dir_in/pred_chembl_${nth}_${dix}.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_test $col_cell $col_drug $col_ic50 $f_name
done

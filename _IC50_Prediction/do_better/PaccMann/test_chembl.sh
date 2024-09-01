#!/usr/bin/bash

ic50=_data/IC50_ChEMBL.txt
cell=data/gene_expression/gdsc-rnaseq_gene-expression.csv
drug=_data/SMILES_ChEMBL.smi

col_cell=Cell_Line_Name
col_drug=Molecule_ChEMBL_ID
col_ic50=LN_IC50
dir_in=Results/IC50_GDSC/Normal

mode=$(seq 0 9)
for nth in ${mode[@]}
do
    f_name=PM_ChEMBL_$nth
    dir_param=$dir_in/param_$nth.pt.tar
    dir_hparam=$dir_in/hyper_param_$nth.json
    dir_test=$dir_in/pred_chembl_$nth.csv
    bash test_write.sh $ic50 $cell $drug $dir_param $dir_hparam $dir_test $col_cell $col_drug $col_ic50 $f_name
done

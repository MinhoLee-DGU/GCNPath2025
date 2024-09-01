#!/usr/bin/bash

choice=0
test_name=chembl
ic50=IC50_ChEMBL.txt

cell=cell_feat.csv
drug=drug_dict_chembl

col_cell=COSMIC_ID
col_drug=Molecule_ChEMBL_ID
col_ic50=LN_IC50

gpu=0
mode=$(seq 0 9)
for nth in ${mode[@]}
do
    bash test_write.sh $ic50 $choice $nth $test_name \
        $cell $drug $col_cell $col_drug $col_ic50 $gpu 
done

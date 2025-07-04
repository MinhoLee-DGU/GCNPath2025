#!/usr/bin/bash

nth=0
choice=0
test_name=chembl
# ic50=IC50_ChEMBL.txt
ic50=IC50_ChEMBL_Time.txt

cell=cell_feat.csv
drug=drug_dict_chembl

col_cell=COSMIC_ID
col_drug=Molecule_ChEMBL_ID
col_ic50=LN_IC50

gpu=1
seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    bash test_write_.sh $ic50 $choice $nth $test_name \
        $cell $drug $col_cell $col_drug $col_ic50 $gpu $seed
done


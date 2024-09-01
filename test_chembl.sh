#!/usr/bin/bash

gpu=1
node=0
test_type=Normal
model_type=RGCN

ic50=data/ic50_data/IC50_ChEMBL.txt
drug=processed/drug_data/ChEMBL_Drug_Custom.pickle
dir_cell=processed/cell_data_biocarta

col_drug=Molecule_ChEMBL_ID
col_ic50=LN_IC50

use_slurm=1
fold_list=$(seq 0 9)
dir_model=results/IC50_GDSC/$test_type/$model_type

for nth in ${fold_list[@]}
do
    param=${dir_model}/param_retrain_${nth}.pt
    hparam=${dir_model}/hyper_param_retrain_${nth}.pickle
    
    col_cell=SANGER_MODEL_ID
    jname=GCN_${model_type}_ChEMBL_${nth}
    out_file=${dir_model}/pred_chembl_${nth}.csv
    cell=${dir_cell}/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
    
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm

    col_cell=BROAD_ID
    jname=GCN_${model_type}_ChEMBL_CCLE_${nth}
    out_file=${dir_model}/pred_chembl_ccle_${nth}.csv
    cell=${dir_cell}/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle
    
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm

    col_cell=SANGER_MODEL_ID
    jname=GCN_${model_type}_ChEMBL_GDSC_${nth}
    out_file=${dir_model}/pred_chembl_gdsc_${nth}.csv
    cell=${dir_cell}/GDSC_RNA_KNN5_STR9_Reg_Corr.pickle
    
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done

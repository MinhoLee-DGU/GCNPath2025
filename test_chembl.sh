#!/usr/bin/bash

gpu=1
node=glu
test_type=Normal
model_type=RGCN

# ic50=data/ic50_data/IC50_ChEMBL.txt
ic50=data/ic50_data/IC50_ChEMBL_Time.txt
drug=processed/drug_data/ChEMBL_Drug_Graph.pickle
dir_cell=processed/cell_data_biocarta

col_ic50=LN_IC50
col_drug=Molecule_ChEMBL_ID
out_cam=None

use_slurm=1
seed_list=$(seq 2021 2030)
dir_model=results/IC50_GDSC/$test_type/$model_type

for seed in ${seed_list[@]} 
do
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    out_cam="None -out_time ${dir_model}/log_time_chembl_seed${seed}.csv"
    
    col_cell=SANGER_MODEL_ID
    jname=ChEMBL_S${seed}
    out_file=${dir_model}/pred_chembl_seed${seed}.csv
    cell=${dir_cell}/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        "$out_cam" $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
    
    # col_cell=BROAD_ID
    # jname=ChEMBL_CCLE_S${seed}
    # out_file=${dir_model}/pred_chembl_ccle_seed${seed}.csv
    # cell=${dir_cell}/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle
    # bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
    #     $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
    
    # col_cell=SANGER_MODEL_ID
    # jname=ChEMBL_GDSC_S${seed}
    # out_file=${dir_model}/pred_chembl_gdsc_seed${seed}.csv
    # cell=${dir_cell}/GDSC_RNA_KNN5_STR9_Reg_Corr.pickle
    # bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
    #     $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done

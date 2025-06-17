#!/usr/bin/bash

gpu=0
node=0
test_type=Normal

cell=processed/cell_data_biocarta/Liu24_invivo_RNA_KNN5_STR9_Reg_Corr.pickle
# cell=processed/cell_data_biocarta/Liu24_invivo_RNA_KNN5_STR9_Reg_Corr_CB.pickle
# cell=processed/cell_data_biocarta/Liu24_invivo_Prot_KNN5_STR9_Reg_Corr.pickle
 
ic50=data/ic50_data/Liu_2024_in_vivo.txt
drug=processed/drug_data/Liu_2024_Drug_Graph.pickle
dir_model=results/IC50_GDSC/$test_type/RGCN

out_cam=None
col_cell=Sample
col_drug=Drug_CID
col_ic50=0

use_slurm=1
seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]}
do
    jname=GCN_${test_type}_Liu_Seed${seed}
    param=${dir_model}/param_retrain_seed${seed}.pt
    hparam=${dir_model}/hyper_param_retrain_seed${seed}.pickle
    
    out_file=${dir_model}/pred_liu24_invivo_seed${seed}.csv
    # out_file=${dir_model}/pred_liu24_invivo_seed${seed}_combat.csv
    # out_file=${dir_model}/pred_liu24_invivo_prot_seed${seed}.csv
    
    bash test_write.sh $ic50 $cell $drug $param $hparam $out_file \
        $out_cam $col_cell $col_drug $col_ic50 $gpu $node $jname $use_slurm
done


#!/usr/bin/bash

ic50_th=$1
choice=$2

gpu=1
node=0
retrain=0
seed=2021

gnn_cell=4
gnn_drug=2
mode_cell=3
mode_drug=3

dir_out=RGCN
dir_cell=processed/cell_data_biocarta
dir_drug=processed/drug_data

case $gnn_cell in
    0) cell=${dir_cell}/SANGER_RNA_Lin.pickle ;;
    # *) cell=${dir_cell}/SANGER_RNA_KNN5_Pert2025.pickle ;;
    *) cell=${dir_cell}/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle ;;
esac

case $gnn_drug in
    0) drug=${dir_drug}/GDSC_Drug_Morgan.pickle ;;
    # 0) drug=${dir_drug}/GDSC_Drug_SMILESVec.pickle ;;
    *) drug=${dir_drug}/GDSC_Drug_Custom.pickle ;;
    # *) drug=${dir_drug}/GDSC_Drug_SA.pickle ;;
esac

case $ic50_th in
    0) IC50=IC50_GDSC.txt ;;
    1) IC50=IC50_GDSC1.txt ;;
    2) IC50=IC50_GDSC2.txt ;;
esac

ic50_=(IC50_GDSC IC50_GDSC1 IC50_GDSC2)
choice_=(Normal Cell_Blind Drug_Blind Strict_Blind)
dir_out=results/${ic50_[$ic50_th]}/${choice_[$choice]}/${dir_out}

ic50=data/ic50_data/$IC50
data="-cell $cell -drug $drug -ic50 $ic50"
option_attn="-attn_mode 0 -dim_attn 32 -n_attn_layer 1 -coef_ffnn 4 -h_attn 4"
option="-gnn_cell $gnn_cell -gnn_drug $gnn_drug -mode_cell $mode_cell -mode_drug $mode_drug"
option="$option -n_hid_cell 3 -n_hid_drug 3 -n_hid_pred 2 -act 3"

case $choice in
    3) fold_list=$(seq 0 24) ;;
    *) fold_list=$(seq 0 9) ;;
esac

use_slurm=1
for nth in ${fold_list[@]}
do
    bash train_write.sh "$data" $dir_out $choice $nth \
        "$option" "$option_attn" $gpu $node \
        $ic50_th $seed $retrain $use_slurm
done

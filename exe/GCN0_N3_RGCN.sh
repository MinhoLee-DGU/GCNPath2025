#!/usr/bin/bash
#SBATCH -J GCN0_N3
#SBATCH -c 4
#SBATCH -p glu
#SBATCH --gres=gpu:1
#SBATCH -x lysine
#SBATCH --mem=4GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out
dir_conda=`conda info --base`
export PATH=$dir_conda/envs/GCNPath/bin
export PYTHONNOUSERSITE=1
python train.py -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle -drug processed/drug_data/GDSC_Drug_Custom.pickle -ic50 data/ic50_data/IC50_GDSC.txt -gnn_cell 4 -gnn_drug 2 -mode_cell 3 -mode_drug 3 -n_hid_cell 3 -n_hid_drug 3 -n_hid_pred 2 -act 3 -attn_mode 0 -dim_attn 32 -n_attn_layer 1 -coef_ffnn 4 -h_attn 4     -out_dir results/IC50_GDSC/Normal/RGCN -choice 0 -nth 3     -cpu 3 -seed_model 2021

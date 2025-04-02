#!/usr/bin/bash

export PYTHONNOUSERSITE=1

undirect=0
path_suf_ori=biocarta
path_suf=biocarta
pathway=data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt

pert_seed=2021
net=data/net_data_${path_suf_ori}/KNN5_STR9_Reg_Corr.csv
omics=processed/cell_data_${path_suf}/SANGER_RNA_GSVA.csv
out=processed/cell_data_${path_suf}/SANGER_RNA_KNN5_Pert${pert_seed}.pickle
python process_cell.py -net $net -omics $omics -out $out -undirect $undirect -perturb_seed $pert_seed


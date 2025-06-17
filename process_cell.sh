#!/usr/bin/bash

undirect=0
path_suf_ori=biocarta
path_suf=biocarta
pathway=data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt

export PYTHONNOUSERSITE=1


### 1-1. Compress RNA data from gene- to pathway-level using GSVA

### SANGER RNA-Seq TPM
# omics=data/cell_data/SANGER_RNA_TPM.csv
# omics=data/cell_data/SANGER_RNA_TPM_Filt.csv
# gsva=processed/cell_data_${path_suf}/SANGER_RNA_GSVA.csv
# Rscript process_cell_gsva.R $omics $pathway $gsva

# ### CCLE RNA-Seq TPM
# omics=data/cell_data/CCLE_RNA_TPM.csv
# gsva=processed/cell_data_${path_suf}/CCLE_RNA_GSVA.csv
# Rscript process_cell_gsva.R $omics $pathway $gsva

# ### GDSC Microarray
# omics=data/cell_data/GDSC_RNA_RMA.csv
# gsva=processed/cell_data_${path_suf}/GDSC_RNA_GSVA.csv
# Rscript process_cell_gsva.R $omics $pathway $gsva

### TCGA RNA-Seq TPM
# omics=data/cell_data/TCGA_RNA_TPM.csv
# gsva=processed/cell_data_${path_suf}/TCGA_RNA_GSVA.csv
# Rscript process_cell_gsva.R $omics $pathway $gsva


### 1-2. Transform GSVA pathway data into graph format
### [-train] Process external data with the scaler fitted by train data 

# net=data/net_data_${path_suf_ori}/KNN5_STR9_Reg_Corr.csv
net=data/net_data_${path_suf_ori}/KNN5_STR9_Reg_Corr.csv
train=processed/cell_data_${path_suf}/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle

### SANGER RNA-Seq TPM
# omics=processed/cell_data_${path_suf}/SANGER_RNA_GSVA.csv
# out=processed/cell_data_${path_suf}/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
# python process_cell.py -net $net -omics $omics -out $out -undirect $undirect 
# -col_weight Weight

### CCLE RNA-Seq TPM
# omics=processed/cell_data_${path_suf}/CCLE_RNA_GSVA.csv
# out=processed/cell_data_${path_suf}/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle
# python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect

### GDSC Microarray RMA
# omics=processed/cell_data_${path_suf}/GDSC_RNA_GSVA.csv
# out=processed/cell_data_${path_suf}/GDSC_RNA_KNN5_STR9_Reg_Corr.pickle
# python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect

### TCGA RNA-Seq TPM
# omics=processed/cell_data_${path_suf}/TCGA_RNA_GSVA.csv
# out=processed/cell_data_${path_suf}/TCGA_RNA_KNN5_STR9_Reg_Corr.pickle
# python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect

### TCGA RNA-Seq TPM [ComBat]
# omics=processed/cell_data_${path_suf}/TCGA_RNA_GSVA_ComBat.csv
# out=processed/cell_data_${path_suf}/TCGA_RNA_KNN5_STR9_Reg_Corr_CB.pickle
# python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect

# ### Liu et al. (2024) RNA-Seq TPM
omics=processed/cell_data_${path_suf}/Liu24_invivo_RNA_GSVA.csv
out=processed/cell_data_${path_suf}/Liu24_invivo_RNA_KNN5_STR9_Reg_Corr.pickle
python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect
# 
# ### Liu et al. (2024) RNA-Seq TPM [ComBat]
omics=processed/cell_data_${path_suf}/Liu24_invivo_RNA_GSVA_ComBat.csv
out=processed/cell_data_${path_suf}/Liu24_invivo_RNA_KNN5_STR9_Reg_Corr_CB.pickle
python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect
# 
# ### Liu et al. (2024) Proteome
omics=processed/cell_data_${path_suf}/Liu24_invivo_Prot_GSVA.csv
out=processed/cell_data_${path_suf}/Liu24_invivo_Prot_KNN5_STR9_Reg_Corr.pickle
python process_cell.py -train $train -net $net -omics $omics -out $out -undirect $undirect
 

#!/usr/bin/env Rscript

source("utils/utils_gsva.R")

args = commandArgs(trailingOnly=TRUE)
default_rna = "data/cell_data/SANGER_RNA_TPM.csv"
default_path = "data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt"
default_out = "processed/cell_data_biocarta/SANGER_RNA_GSVA.csv"
# If the format of the RNA is Genes x Cell, please set cell_row=F

dir_rna = args_def(args[1], default_rna)     # Directory of Input RNA CSV file
dir_path = args_def(args[2], default_path)   # Directory of Input Pathway GMT file
dir_out = args_def(args[3], default_out)     # Directory of Output GSVA CSV file
cell_row = args_def(args[4], T)              # TPM row represent cells/genes [T/F]? (Default : TRUE)
fill_na = args_def(args[5], 0)               # Value to fill missing values in RNA data (Default : 0)
cores = args_def(args[6], T)                 # CPU threads to run GSVA (Default : TRUE, a half of threads)


# Implement GSVA in RNA Data
Path_List = read_gmt(dir_path)

RNA = read.csv(dir_rna, row.names=1, check.names=F)
RNA[is.na(RNA)] = fill_na

RNA_GSVA = gsva_def(RNA, Path_List, filt_genes=F, method="gsva", cores=cores, cell_row=cell_row)

dir_out_ = strsplit(dir_out, "/")[[1]]
dir_out_ = paste0(dir_out_[1:(length(dir_out_)-1)], collapse="/")
dir.create(dir_out_, recursive=T)

write.csv(RNA_GSVA, file=dir_out, row.names=T)


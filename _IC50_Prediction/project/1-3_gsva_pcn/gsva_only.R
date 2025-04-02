#!/usr/bin/env Rscript

suppressMessages(library(cogena))

source("../functions.R")
loadings()

args_def = function(args, value_def=NULL) {
  args_x = args
  if (toupper(args_x) %in% c("T", "TRUE")) args_x = T
  if (toupper(args_x) %in% c("F", "FALSE")) args_x = F
  args_x = ifelse(!is.na(args_x), args_x, value_def)
  return(args_x)
}

dir = "../../processed_data/path_data"
rna_def = sprintf("%s/KEGG_Total_Ent.gmt", dir)

dir = "../../processed_data/cell_data/SANGER_Passports"
path_def = sprintf("%s/TPM_Ent.csv", dir)

dir = mkdir("../../processed_data/cell_data/KEGG")
out_def = sprintf("%s/SANGER_RNA_GSVA.csv", dir)

args = commandArgs(trailingOnly=TRUE)
# If the format of the RNA is Genes x Cell, please set cell_row=F

dir_rna = args_def(args[1], rna_def)
dir_path = args_def(args[2], path_def)
dir_out = args_def(args[3], out_def)
cell_row = args_def(args[4], T)
cores = args_def(args[5], T)



# Implement GSVA in SANGER RNA-Seq Data

gsva_def = function(Omics, Path_List, method="gsva", cell_row=T,
                    filt_genes=T, do_qnorm=F, cores=T, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cells x Genes if cell_row=T]
  # Input : Omics [Genes x Cells if cell_row=F]
  # Output : Omics [Cells x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  suppressMessages(library(stats))
  gene_list = Path_List %>% unlist %>% unique
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  
  n_before = length(gene_list)
  if (cell_row) {
    n_after = sum(colnames(Omics) %in% gene_list)
  } else n_after = sum(rownames(Omics) %in% gene_list)
  
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    if (cell_row) {
      n_before = ncol(Omics)
      Omics = Omics[, colnames(Omics) %in% gene_list]
    } else {
      n_before = nrow(Omics)
      Omics = Omics[rownames(Omics) %in% gene_list, ]
    }
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  # GSVA input Omics as the format [Genes x Cells]
  if (cell_row) Omics = Omics %>% t %>% as.matrix
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = Omics %>% normalize.quantiles(keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = Omics_Path %>% t %>% as.data.frame
  return(Omics_Path)
}

Path_List = gmt2list(dir_path)
RNA = read.csv(dir_rna, row.names=1, check.names=F)
RNA_GSVA = RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", cell_row=cell_row, cores=cores)
write.csv(RNA_GSVA, file=dir_out, row.names=T)



# # Implement GSVA in CCLE RNA-Seq Data
# dir = "../../processed_data/path_data"
# rna_def = sprintf("%s/KEGG_Total_Ent.gmt", dir)
# 
# dir = "../../processed_data/cell_data/CCLE_DepMap"
# path_def = sprintf("%s/TPM_BROAD_ID_Ent.csv", dir)
# 
# dir = mkdir("../../processed_data/cell_data/KEGG")
# out_def = sprintf("%s/CCLE_RNA_GSVA_BROAD.csv", dir)

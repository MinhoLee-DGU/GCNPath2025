#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))

source("../functions.R")
loadings()

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
  if (cell_row) {
    Omics = Omics %>% t %>% as.matrix
  } else Omics = Omics %>% as.matrix
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = Omics %>% normalize.quantiles(keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = Omics_Path %>% t %>% as.data.frame
  return(Omics_Path)
}

singscore_def = function(Omics, Path_List, path_names=NULL, filt_genes=T, 
                         direction=F, center=F, stableGenes=NULL, ...) {
  
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  
  suppressMessages(library(GSEABase))
  suppressMessages(library(singscore))
  gene_list = Path_List %>% unlist %>% unique
  
  n_before = length(gene_list)
  n_after = sum(colnames(Omics) %in% gene_list)
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# The omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    n_before = ncol(Omics)
    Omics = Omics[, colnames(Omics) %in% gene_list]
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  if (is.null(path_names)) path_names = names(Path_List)
  Path_List = Path_List %>% lapply(GeneSet)
  Omics = Omics %>% t %>% rankGenes(stableGenes=stableGenes)
  for (i in 1:length(path_names)) Path_List[[i]]@setName = path_names[[i]]
  
  Omics_Path = multiScore(Omics, upSetColc=Path_List, knownDirection=direction, centerScore=center, ...)
  Omics_Path = Omics_Path[[1]] %>% t %>% as.data.frame
  return(Omics_Path)
}

examine_cell_pair = function(Omics_Path1, Omics_Path2, method, description=NULL, cores=T) {
  
  nrmse = function(x, y) RMSE_Norm(x, y, na.rm=T)
  Corr_Info = corr_pair(Omics_Path1, Omics_Path2, by_row=T, into_wide=F)
  RMSE_Info = stat_pair(Omics_Path1, Omics_Path2, by_row=T, into_wide=F, method=nrmse, cores=cores)
  
  by = c("Cell1"="Cell1", "Cell2"="Cell2")
  colnames(Corr_Info) = c("Cell1", "Cell2", "Corr")
  colnames(RMSE_Info) = c("Cell1", "Cell2", "NRMSE")
  
  description = ifelse(!is.null(description), description, method)
  Cell_Info = merge(RMSE_Info, Corr_Info, by=by)
  Cell_Info = Cell_Info %>% mutate(Method=method, Description=description)
  
  return(Cell_Info)
}

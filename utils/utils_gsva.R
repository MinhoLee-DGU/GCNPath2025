#!/usr/bin/env Rscript

args_def = function(args, default=NULL) {
  args_x = args
  if (toupper(args_x) %in% c("T", "TRUE")) args_x = T
  if (toupper(args_x) %in% c("F", "FALSE")) args_x = F
  args_x = ifelse(!is.na(args_x), args_x, default)
  return(args_x)
}

read_gmt = function(file) {
  Path_List = scan(file, what="", sep="\n", quiet=T)
  Path_List = strsplit(Path_List, "\t")
  names(Path_List) = sapply(Path_List, `[[`, 1)
  Path_List = lapply(Path_List, `[`, c(-1, -2))
  return(Path_List)
}

gsva_def = function(Omics, Path_List, method="gsva", cell_row=T,
                    filt_genes=T, do_qnorm=F, cores=T, fill_na=0, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cells x Genes if cell_row=T]
  # Input : Omics [Genes x Cells if cell_row=F]
  # Output : Omics [Cells x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  gene_list = unique(unlist(Path_List))
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  
  n_before = length(gene_list)
  if (cell_row) {
    n_after = sum(colnames(Omics) %in% gene_list)
  } else n_after = sum(rownames(Omics) %in% gene_list)
  
  int_ratio = round(100 * n_after/n_before, 2)
  print(sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before))
  
  if (filt_genes) {
    if (cell_row) {
      n_before = ncol(Omics)
      Omics = Omics[, colnames(Omics) %in% gene_list]
    } else {
      n_before = nrow(Omics)
      Omics = Omics[rownames(Omics) %in% gene_list, ]
    }
    print(sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after))
  } else print("# Genes are not filtered from omics data...")
  
  # GSVA input Omics as the format [Genes x Cells]
  if (cell_row) Omics = as.matrix(t(Omics))
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = normalize.quantiles(Omics, keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = as.data.frame(t(Omics_Path))
  
  if (all(names(Path_List) %in% colnames(Omics_Path))) {
    print(sprintf("# All %s pathways are processed...", length(Path_List)))
  } else {
    col_na = setdiff(names(Path_List), colnames(Omics_Path))
    Omics_Path[, col_na] = fill_na
    Omics_Path = Omics_Path[, names(Path_List)]
    
    print(sprintf("# The following %s pathways are not processed...", length(col_na)))
    print(sprintf("# Those pathways in all samples were filled with %s...", fill_na))
    print("# Utilizing the data in prediction is not recommended...")
    print("# We will add imputation function in upcoming version...")
    print(sprintf(" > %s", col_na))
  }
  
  return(Omics_Path)
}

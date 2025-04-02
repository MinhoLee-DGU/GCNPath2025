#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(ggpubr))
suppressMessages(library(cogena))
suppressMessages(library(graphite))
suppressMessages(library(reshape2))

source("../functions.R")
loadings()

args = commandArgs(trailingOnly=TRUE)
choice = args[1] %>% as.numeric

file = c("../../processed_data/cell_data/SANGER_Passports/TPM_Ent.csv", 
         "../../processed_data/cell_data/CCLE_DepMap/TPM_Ent.csv", 
         "../../processed_data/cell_data/GDSC/RNA_Array_Ent.csv", 
         "../../processed_data/cell_data/TCGA/RNA_Ent.csv")[choice]

omics = c("SANGER", "CCLE", "GDSC", "TCGA")[choice]

RNA = fread_def(file, check_names=F, header=T)

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.entrez.gmt", dir)
Path_List = gmt2list(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

dir_out = mkdir("../../processed_data/cell_data/BIOCARTA/Gene_Filtering/%s", omics)


gsva_def = function(Omics, Path_List, method="gsva", 
                    filt_genes=T, do_qnorm=F, cores=10, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  suppressMessages(library(stats))
  gene_list = Path_List %>% unlist %>% unique
  
  n_before = length(gene_list)
  n_after = sum(colnames(Omics) %in% gene_list)
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    n_before = ncol(Omics)
    Omics = Omics[, colnames(Omics) %in% gene_list]
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  Omics = Omics %>% t %>% as.matrix
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = Omics %>% normalize.quantiles(keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = Omics_Path %>% t %>% as.data.frame
  return(Omics_Path)
}

singscore_def = function(Omics, Path_List, path_names=NULL, filt_genes=T, 
                         direction=T, center=T, stableGenes=NULL, ...) {
  
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



### GSVA vs ssGSEA & singscore

rna_genes = RNA %>% colnames
path_genes = Path_List %>% unlist %>% unique
sum(path_genes %in% rna_genes)   # 99.80% [1506/1509]

Path_List %>% sapply(length) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))/length(x)) %>% hist
# [5, 81], [5, 81], Mostly above 0.95


### 2-1. GSVA is sensitive to the number of total genes in omics data? [X]

examine_cell_pair = function(Omics_Path_Ref, Omics_Path, col_common=NULL, method=NULL, cores=T) {
  
  nrmse = function(x, y) RMSE_Norm(x, y, na.rm=T)
  if (is.null(col_common)) col_common = colnames(Omics_Path)
  Corr_Info = corr_pair(Omics_Path[, col_common], Omics_Path_Ref[, col_common], by_row=T, into_wide=F)
  RMSE_Info = stat_pair(Omics_Path[, col_common], Omics_Path_Ref[, col_common], by_row=T, into_wide=F, method=nrmse, cores=cores)
  
  by = c("Cell1"="Cell1", "Cell2"="Cell2")
  colnames(Corr_Info) = c("Cell1", "Cell2", "Corr")
  colnames(RMSE_Info) = c("Cell1", "Cell2", "NRMSE")
  
  Cell_Info = merge(RMSE_Info, Corr_Info, by=by) %>% mutate(Method=method)
  return(Cell_Info)
}


# Define stable genes as those with lowest sd
stable = getStableGenes(n_stable=3000, type="carcinoma")   # 1000
idx = match(stable, Anno_Genes$HGNC_SYMBOL)
stable = Anno_Genes$ENTREZ_ID[idx] %>% na.omit   # 1000

sum(stable %in% colnames(RNA))   # 1000
sum(stable %in% path_genes)      # 119

stable_100 = stable %>% 
  intersect(path_genes) %>% 
  intersect(colnames(RNA)) %>% head(100)   # 100

# match(stable_100, stable)   # 56, ..., 937
common_genes = intersect(path_genes, colnames(RNA))

RNA_ZNorm = RNA[, common_genes] %>% t %>% scale %>% t %>% as.data.frame
RNA_ssGSEA = RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
RNA_ssGSEA_NN = RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
RNA_GSVA = RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
RNA_GSVA_QN = RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", do_qnorm=T)
RNA_SING = RNA %>% singscore_def(Path_List, filt_genes=T)
RNA_SING_ST = RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)

RNA_ZNorm_NF = RNA %>% t %>% scale %>% t %>% as.data.frame
RNA_ssGSEA_NF = RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea")
RNA_ssGSEA_NN_NF = RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea", ssgsea.norm=F)
RNA_GSVA_NF = RNA %>% gsva_def(Path_List, filt_genes=F, method="gsva")
RNA_GSVA_QN_NF = RNA %>% gsva_def(Path_List, filt_genes=F, method="gsva", do_qnorm=T)
RNA_SING_NF = RNA %>% singscore_def(Path_List, filt_genes=F)
RNA_SING_ST_NF = RNA %>% singscore_def(Path_List, filt_genes=F, stableGenes=stable_100)

RNA_ZNorm_NF = RNA_ZNorm_NF[, common_genes]

file = sprintf("%s/Omics.RData", dir_out)
save(RNA_ZNorm, RNA_ssGSEA, RNA_ssGSEA_NN, RNA_GSVA, RNA_GSVA_QN, RNA_SING, RNA_SING_ST, 
     RNA_ZNorm_NF, RNA_ssGSEA_NF, RNA_ssGSEA_NN_NF, RNA_GSVA_NF, RNA_GSVA_QN_NF, RNA_SING_NF, RNA_SING_ST_NF, file=file)


methods = c("Z_Sample", "ssGSEA", "ssGSEA_Norm_X", 
            "GSVA", "GSVA_QN", "singscore", "stingscore")

# Set the filtering non-pathway genes as the reference data
Corr_Filt_ZNorm = examine_cell_pair(RNA_ZNorm, RNA_ZNorm_NF, methods[1])
Corr_Filt_ssGSEA = examine_cell_pair(RNA_ssGSEA, RNA_ssGSEA_NF, methods[2])
Corr_Filt_ssGSEA_NN = examine_cell_pair(RNA_ssGSEA_NN, RNA_ssGSEA_NN_NF, methods[3])
Corr_Filt_GSVA = examine_cell_pair(RNA_GSVA, RNA_GSVA_NF, methods[4])
Corr_Filt_GSVA_QN = examine_cell_pair(RNA_GSVA_QN, RNA_GSVA_QN_NF, methods[5])
Corr_Filt_SING = examine_cell_pair(RNA_SING, RNA_SING_NF, methods[6])
Corr_Filt_SING_ST = examine_cell_pair(RNA_SING_ST, RNA_SING_ST_NF, methods[7])

file = sprintf("%s/Corr_Filt.RData", dir_out)
save(Corr_Filt_ZNorm, Corr_Filt_ssGSEA, Corr_Filt_ssGSEA_NN, 
     Corr_Filt_GSVA, Corr_Filt_GSVA_QN, Corr_Filt_SING, Corr_Filt_SING_ST, file=file)



# Example [SIDM00001]
Ex_Filt_ZNorm = data.frame(Filt_O=as.numeric(RNA_ZNorm[1, ]), 
                           Filt_X=as.numeric(RNA_ZNorm_NF[1, ]))
Ex_Filt_ssGSEA = data.frame(Filt_O=as.numeric(RNA_ssGSEA[1, ]), 
                            Filt_X=as.numeric(RNA_ssGSEA_NF[1, ]))
Ex_Filt_ssGSEA_NN = data.frame(Filt_O=as.numeric(RNA_ssGSEA_NN[1, ]), 
                               Filt_X=as.numeric(RNA_ssGSEA_NN_NF[1, ]))
Ex_Filt_GSVA = data.frame(Filt_O=as.numeric(RNA_GSVA[1, ]), 
                          Filt_X=as.numeric(RNA_GSVA_NF[1, ]))
Ex_Filt_GSVA_QN = data.frame(Filt_O=as.numeric(RNA_GSVA_QN[1, ]), 
                             Filt_X=as.numeric(RNA_GSVA_QN_NF[1, ]))
Ex_Filt_SING = data.frame(Filt_O=as.numeric(RNA_SING[1, ]), 
                          Filt_X=as.numeric(RNA_SING_NF[1, ]))
Ex_Filt_SING_ST = data.frame(Filt_O=as.numeric(RNA_SING_ST[1, ]), 
                             Filt_X=as.numeric(RNA_SING_ST_NF[1, ]))


cell_name = rownames(RNA_SING_ST)[1]
methods = c("ssGSEA", "GSVA", "singscore")
main = sprintf("%s/Example Scatter [%s, SIDM00001, Genes filtering]", dir, methods)

plot_ex_filt = function(Ex_Filt, main=NULL, cell_name=NULL, save=T, ...) {
  xlab = sprintf("%s [Pathway Genes]", )
  ylab = sprintf("%s [Total Genes]", )
  Ex_Filt %>% plot_def(Filt_O, Filt_X, main=main, xlab=xlab, ylab=ylab, xy_line=T, save=save, ...)
}

Ex_Filt_ZNorm %>% plot_ex_filt(main=main[1], cell_name=cell_name)
Ex_Filt_ssGSEA %>% plot_ex_filt(main=main[2], cell_name=cell_name)
Ex_Filt_ssGSEA_NN %>% plot_ex_filt(main=main[3], cell_name=cell_name)
Ex_Filt_GSVA %>% plot_ex_filt(main=main[4], cell_name=cell_name)
Ex_Filt_GSVA_QN %>% plot_ex_filt(main=main[5], cell_name=cell_name)
Ex_Filt_SING %>% plot_ex_filt(main=main[6], cell_name=cell_name)
Ex_Filt_SING_ST %>% plot_ex_filt(main=main[7], cell_name=cell_name)

file = sprintf("%s/Example_Filt.RData", dir_out)
save(Ex_Filt_ZNorm, Ex_Filt_ssGSEA, Ex_Filt_ssGSEA_NN, 
     Ex_Filt_GSVA, Ex_Filt_GSVA_QN, Ex_Filt_SING, Ex_Filt_SING_ST, file=file)

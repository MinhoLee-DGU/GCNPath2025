#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(igraph))
suppressMessages(library(ggpubr))
suppressMessages(library(cogena))
suppressMessages(library(reshape2))

source("../functions.R")
loadings()


dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
SANGER_RNA = fread_def(file, check_names=F, header=T)

dir = "../../processed_data/cell_data/CCLE_DepMap"
file = sprintf("%s/TPM_Ent.csv", dir)
CCLE_RNA = fread_def(file, check_names=F, header=T)

dir = "../../processed_data/cell_data/GDSC"
file = sprintf("%s/RNA_Array_Ent.csv", dir)
GDSC_RNA = fread_def(file, check_names=F, header=T)

tcga = T
# Details are listed in project/4-2_tcga_response/code.R

if (tcga) {
  dir = "../../processed_data/cell_data/TCGA"
  file = sprintf("%s/TPM_Ent.csv", dir)
  TCGA_RNA = fread_def(file, check_names=F, header=T)
}

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.entrez.gmt", dir)
Path_List = gmt2list(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/Anno_Cells.csv", dir)
# Anno_Cells = read.csv(file)



### Define statistic functions
# Those functions are already defined in 0_functions.R

# RMSE_Norm [NRMSE-IQR]
# RNSE Normalized by the interquartile range [IQR] of observations

# stat_pair
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of two data.frames in pairwise manner

# stat_pair_apply
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of two data.frames in matched rows or columns

# stat_self
# Apply a statistic function (such as RMSE, R2, Euclidean distance)
# to all rows or columns of one data.frame in pairwise manner

# corr_self
# Apply a correlation function (method="pearson", use="pairwise.complete.obs")
# to all rows or columns of one data.frame in pairwise manner



##### 2. Calculate pathway activity scores [GSVA]

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

rna_genes = SANGER_RNA %>% colnames
path_genes = Path_List %>% unlist %>% unique
sum(path_genes %in% rna_genes)   # 99.80% [1506/1509]

Path_List %>% sapply(length) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))) %>% range
Path_List %>% sapply(function(x) length(intersect(x, rna_genes))/length(x)) %>% hist
# [5, 81], [5, 81], Mostly above 0.95


### 2-1. GSVA is sensitive to the number of total genes in omics data? [X]

examine_cell_pair = function(Omics_Path_Ref, Omics_Path, method=NULL, col_common=NULL, calc_rmse=F, cores=T) {
  
  if (is.null(col_common)) col_common = colnames(Omics_Path)
  Corr_Info = corr_pair(Omics_Path[, col_common], Omics_Path_Ref[, col_common], by_row=T, into_wide=F)
  colnames(Corr_Info) = c("Cell1", "Cell2", "Corr")
  Cell_Info = Corr_Info
  
  if (calc_rmse) {
    nrmse = function(x, y) RMSE_Norm(x, y, na.rm=T)
    RMSE_Info = stat_pair(Omics_Path[, col_common], Omics_Path_Ref[, col_common], by_row=T, into_wide=F, method=nrmse, cores=cores)
    colnames(RMSE_Info) = c("Cell1", "Cell2", "NRMSE")
    by = c("Cell1"="Cell1", "Cell2"="Cell2")
    Cell_Info = merge(RMSE_Info, Corr_Info, by=by)
  } 
  
  Cell_Info = Cell_Info %>% mutate(Method=method)
  return(Cell_Info)
}

examine_cell_each = function(Omics_Path_Ref, Omics_Path, method=NULL, col_common=NULL) {
  
  if (!is.null(col_common)) {
    col = col_common
  } else col = intersect(colnames(Omics_Path), colnames(Omics_Path_Ref))
  row = intersect(rownames(Omics_Path), rownames(Omics_Path_Ref))
  
  # Omics_Path = Omics_Path[row, col]
  # Omics_Path_Ref = Omics_Path_Ref[row, col]
  Omics_Path = Omics_Path[row, col, drop=F] %>% as.data.table
  Omics_Path_Ref = Omics_Path_Ref[row, col, drop=F] %>% as.data.table
  
  corr = sapply(1:length(row), function(i) {
    cor(as.numeric(Omics_Path_Ref[i, ]), 
        as.numeric(Omics_Path[i, ]))
  })
  
  Cell_Info = data.frame(Cell=row, Corr=corr, Method=method)
  return(Cell_Info)
}


# Define stable genes as those with lowest sd
stable = getStableGenes(n_stable=1000, type="carcinoma")   # 1000
idx = match(stable, Anno_Genes$HGNC_SYMBOL)
stable = Anno_Genes$ENTREZ_ID[idx] %>% na.omit   # 1000

sum(stable %in% colnames(SANGER_RNA))   # 1000
sum(stable %in% colnames(GDSC_RNA))     # 953
sum(stable %in% colnames(CCLE_RNA))     # 999
sum(stable %in% colnames(TCGA_RNA))     # 1000
sum(stable %in% path_genes)             # 119

stable_100 = stable %>% 
  intersect(path_genes) %>% 
  intersect(colnames(SANGER_RNA)) %>% head(100)   # 100

scale_robust = function(x, na.rm=T, constant=1.4826) {
  (x-median(x, na.rm=na.rm)) / mad(x, constant=constant, na.rm=na.rm)
}

# match(stable_100, stable)   # 56, ..., 937
common_genes = intersect(path_genes, colnames(SANGER_RNA))   # 1506

SANGER_ZNorm = SANGER_RNA[, common_genes] %>% t %>% scale %>% t %>% as.data.frame
SANGER_ZNormR = SANGER_RNA[, common_genes] %>% apply(1, scale_robust) %>% t %>% as.data.frame
SANGER_ssGSEA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
SANGER_ssGSEA_NN = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
SANGER_GSVA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
SANGER_GSVA_QN = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", do_qnorm=T)
SANGER_SING = SANGER_RNA %>% singscore_def(Path_List, filt_genes=T)
SANGER_SING_ST = SANGER_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)

SANGER_ZNorm_NF = SANGER_RNA %>% t %>% scale %>% t %>% as.data.frame
SANGER_ZNormR_NF = SANGER_RNA %>% apply(1, scale_robust) %>% t %>% as.data.frame
SANGER_ssGSEA_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea")
SANGER_ssGSEA_NN_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="ssgsea", ssgsea.norm=F)
SANGER_GSVA_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="gsva")
SANGER_GSVA_QN_NF = SANGER_RNA %>% gsva_def(Path_List, filt_genes=F, method="gsva", do_qnorm=T)
SANGER_SING_NF = SANGER_RNA %>% singscore_def(Path_List, filt_genes=F)
SANGER_SING_ST_NF = SANGER_RNA %>% singscore_def(Path_List, filt_genes=F, stableGenes=stable_100)

identical(as.numeric(SANGER_RNA[1, common_genes]) %>% scale_robust, as.numeric(SANGER_ZNormR[1, ]))   # T
identical(as.numeric(SANGER_RNA[1, ]) %>% scale_robust, as.numeric(SANGER_ZNormR_NF[1, ]))   # T

methods = c("Z_Sample", "ZRobust_Sample", "ssGSEA", "ssGSEA_Norm_X", 
            "GSVA", "GSVA_QN", "singscore", "stingscore")

# Set the filtering non-pathway genes as the reference data
Corr_Filt_ZNorm = examine_cell_each(SANGER_ZNorm, SANGER_ZNorm_NF[, common_genes], methods[1])
Corr_Filt_ZNormR = examine_cell_each(SANGER_ZNormR, SANGER_ZNormR_NF[, common_genes], methods[2])
Corr_Filt_ssGSEA = examine_cell_each(SANGER_ssGSEA, SANGER_ssGSEA_NF, methods[3])
Corr_Filt_ssGSEA_NN = examine_cell_each(SANGER_ssGSEA_NN, SANGER_ssGSEA_NN_NF, methods[4])
Corr_Filt_GSVA = examine_cell_each(SANGER_GSVA, SANGER_GSVA_NF, methods[5])
Corr_Filt_GSVA_QN = examine_cell_each(SANGER_GSVA_QN, SANGER_GSVA_QN_NF, methods[6])
Corr_Filt_SING = examine_cell_each(SANGER_SING, SANGER_SING_NF, methods[7])
Corr_Filt_SING_ST = examine_cell_each(SANGER_SING_ST, SANGER_SING_ST_NF, methods[8])

Corr_Filt = Reduce(rbind, list(Corr_Filt_ZNorm, Corr_Filt_ZNormR, 
                               Corr_Filt_ssGSEA, Corr_Filt_ssGSEA_NN, 
                               Corr_Filt_GSVA, Corr_Filt_GSVA_QN, 
                               Corr_Filt_SING, Corr_Filt_SING_ST))

Corr_Filt$Method=factor(Corr_Filt$Method, level=methods)


width = 20
height = 16
legend = "Cell-Pair"

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Gene Filtering")
file = sprintf("%s/Cell-Pair Comparison [Gene Filtering, PCC]", dir)

ylab = "Cell-Pair PCC"
pos = position_dodge(width=0.85)
# color = RColorBrewer::brewer.pal(8, "Reds")[8]

Corr_Filt %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Method, main=file, pos=pos, 
              ylab=ylab, width=width, height=height, legend=F,
              alpha=0.9, hjust=1, vjust=1, dpi=600, axis_tl=30, axis_tx=20, 
              point=F, violin=F, raster=T, save=T, save_svg=T)


# Example [SIDM00001]
Ex_Filt_ZNorm = data.frame(Filt_O=as.numeric(SANGER_ZNorm[1, ]), 
                           Filt_X=as.numeric(SANGER_ZNorm_NF[1, common_genes]))
Ex_Filt_ZNormR = data.frame(Filt_O=as.numeric(SANGER_ZNormR[1, ]), 
                            Filt_X=as.numeric(SANGER_ZNormR_NF[1, common_genes]))
Ex_Filt_ssGSEA = data.frame(Filt_O=as.numeric(SANGER_ssGSEA[1, ]), 
                            Filt_X=as.numeric(SANGER_ssGSEA_NF[1, ]))
Ex_Filt_ssGSEA_NN = data.frame(Filt_O=as.numeric(SANGER_ssGSEA_NN[1, ]), 
                               Filt_X=as.numeric(SANGER_ssGSEA_NN_NF[1, ]))
Ex_Filt_GSVA = data.frame(Filt_O=as.numeric(SANGER_GSVA[1, ]), 
                          Filt_X=as.numeric(SANGER_GSVA_NF[1, ]))
Ex_Filt_GSVA_QN = data.frame(Filt_O=as.numeric(SANGER_GSVA_QN[1, ]), 
                             Filt_X=as.numeric(SANGER_GSVA_QN_NF[1, ]))
Ex_Filt_SING = data.frame(Filt_O=as.numeric(SANGER_SING[1, ]), 
                          Filt_X=as.numeric(SANGER_SING_NF[1, ]))
Ex_Filt_SING_ST = data.frame(Filt_O=as.numeric(SANGER_SING_ST[1, ]), 
                             Filt_X=as.numeric(SANGER_SING_ST_NF[1, ]))

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Example Scatter Plot [Gene Filtering]")
main = sprintf("%s/Example Scatter [%s, SIDM00001, Genes filtering]", dir, methods)

plot_ex_filt = function(Filt_Ex, main=NULL, save=T, ...) {
  xlab = "SIDM00001 [Pathway Genes]"
  ylab = "SIDM00001 [Total Genes]"
  Filt_Ex %>% plot_def(Filt_O, Filt_X, main=main, xlab=xlab, ylab=ylab, xy_line=T, save=save, ...)
}

Ex_Filt_ZNorm %>% plot_ex_filt(main=main[1])
Ex_Filt_ZNormR %>% plot_ex_filt(main=main[2])
Ex_Filt_ssGSEA %>% plot_ex_filt(main=main[3])
Ex_Filt_ssGSEA_NN %>% plot_ex_filt(main=main[4])
Ex_Filt_GSVA %>% plot_ex_filt(main=main[5])
Ex_Filt_GSVA_QN %>% plot_ex_filt(main=main[6])
Ex_Filt_SING %>% plot_ex_filt(main=main[7])
Ex_Filt_SING_ST %>% plot_ex_filt(main=main[8])



### 2-2. GSVA can calculate the pathway scores differently in different cells? [O]
# Except GSVA, all methods process cells in similar values...

# examine_cell = function(Omics_Path, method) {
#   
#   nrmse = function(x, y) RMSE_Norm(x, y, na.rm=T)
#   Corr_Info = Omics_Path %>% corr_self(by_row=T, into_wide=F)
#   RMSE_Info = Omics_Path %>% stat_self(by_row=T, into_wide=F, method=nrmse)
#   
#   by = c("Cell1"="Cell1", "Cell2"="Cell2")
#   colnames(Corr_Info) = c("Cell1", "Cell2", "Corr")
#   colnames(RMSE_Info) = c("Cell1", "Cell2", "NRMSE")
#   
#   Cell_Info = merge(RMSE_Info, Corr_Info, by=by) %>% mutate(Method=method)
#   return(Cell_Info)
# }
# 
# Corr_Diff_ZNorm = SANGER_ZNorm %>% examine_cell(methods[1])
# Corr_Diff_ZNormR = SANGER_ZNormR %>% examine_cell(methods[2])
# Corr_Diff_ssGSEA = SANGER_ssGSEA %>% examine_cell(methods[3])
# Corr_Diff_ssGSEA_NN = SANGER_ssGSEA_NN %>% examine_cell(methods[4])
# Corr_Diff_GSVA = SANGER_GSVA %>% examine_cell(methods[5])
# Corr_Diff_GSVA_QN = SANGER_GSVA_QN %>% examine_cell(methods[6])
# Corr_Diff_SING = SANGER_SING %>% examine_cell(methods[7])
# Corr_Diff_SING_ST = SANGER_SING_ST %>% examine_cell(methods[8])
# 
# 
# Corr_Diff = Reduce(rbind, list(Corr_Diff_ZNorm, Corr_Diff_ZNormR, Corr_Diff_ssGSEA, Corr_Diff_ssGSEA_NN, 
#                                Corr_Diff_GSVA, Corr_Diff_GSVA_QN, Corr_Diff_SING, Corr_Diff_SING_ST))
# 
# Corr_Diff = Corr_Diff %>% 
#   mutate(Method = Method %>% factor(levels=methods))
# 
# 
# stat = c("NRMSE", "Corr")
# ylab = sprintf("%s [Different cell-pairs]", stat)
# dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Cell Distinction")
# file = sprintf("%s/Cell-Pair Comparison [%s, SANGER]", dir, stat)
# 
# Corr_Diff %>% subset(Cell1!=Cell2) %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Method, main=file[1], dpi=1500, 
#               ylab=ylab[1], point=F, violin=F, legend=F, save=T, save_svg=T, raster=T)
# 
# Corr_Diff %>% subset(Cell1!=Cell2) %>% as.data.frame %>% 
#   boxplot_def(Method, Corr, fill=Method, main=file[2], dpi=1500, 
#               ylab=ylab[2], point=F, violin=F, legend=F, save=T, save_svg=T, raster=T)


# Example
plot_ex_diff = function(Omics, main=NULL, width=15, height=15, save=T, ...) {
  
  xlab = "SIDM00001 [MEC-1, CLL]"
  ylab = "SIDM00002 [NBsusSR, NB]"
  pcc = Omics %>% t %>% as.data.frame %>% with(cor(SIDM00001, SIDM00002)) %>% round(3)
  sprintf("PCC [SIDM00001-SIDM00002] : %s", pcc) %>% print
  
  Omics %>% t %>% as.data.frame %>% 
    plot_def(SIDM00001, SIDM00002, main=main, alpha=0.5, size=2,
             xlab=xlab, ylab=ylab, width=width, height=height, 
             text_ratio=1.35, xy_line=T, save=save, save_svg=save, raster=T)
}

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Example Scatter Plot [Cell Distinction]")
main = sprintf("%s/Example Scatter [%s, SIDM00001 & SIDM00002]", dir, methods)

SANGER_ZNorm %>% plot_ex_diff(main=main[1])
SANGER_ZNormR %>% plot_ex_diff(main=main[2])
SANGER_ssGSEA %>% plot_ex_diff(main=main[3])
SANGER_ssGSEA_NN %>% plot_ex_diff(main=main[4])
SANGER_GSVA %>% plot_ex_diff(main=main[5])
SANGER_GSVA_QN %>% plot_ex_diff(main=main[6])
SANGER_SING %>% plot_ex_diff(main=main[7])
SANGER_SING_ST %>% plot_ex_diff(main=main[8])

# file = sprintf("%s/Cell_Distinction.csv", dir)
# fwrite(Corr_Diff, file=file)



### 2-3-1. GSVA can alleviate the batch-effects? [O]
# SANGER RNA-Seq vs CCLE RNA-Seq

all(stable_100 %in% colnames(GDSC_RNA))     # 93 [from 100]
all(stable_100 %in% colnames(CCLE_RNA))     # 100 [from 100]
all(common_genes %in% colnames(GDSC_RNA))   # 1410 [from 1506]
all(common_genes %in% colnames(CCLE_RNA))   # 1498 [from 1506]

common_gdsc = Reduce(intersect, list(path_genes, colnames(SANGER_RNA), colnames(GDSC_RNA)))
common_ccle = Reduce(intersect, list(path_genes, colnames(SANGER_RNA), colnames(CCLE_RNA)))

# GDSC Microarray 90.92% [4879/5366]
GDSC_ZNorm = GDSC_RNA[, common_gdsc] %>% t %>% scale %>% t %>% as.data.frame
GDSC_ZNormR = GDSC_RNA[, common_gdsc] %>% apply(1, scale_robust) %>% t %>% as.data.frame
GDSC_ssGSEA = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
GDSC_ssGSEA_NN = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
GDSC_GSVA = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
GDSC_GSVA_QN = GDSC_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", do_qnorm=T)
GDSC_SING = GDSC_RNA %>% singscore_def(Path_List, filt_genes=T)
GDSC_SING_ST = GDSC_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)

# CCLE RNA-Seq 99.20% [5323/5366]
CCLE_ZNorm = CCLE_RNA[, common_ccle] %>% t %>% scale %>% t %>% as.data.frame
CCLE_ZNormR = CCLE_RNA[, common_ccle] %>% apply(1, scale_robust) %>% t %>% as.data.frame
CCLE_ssGSEA = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea")
CCLE_ssGSEA_NN = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="ssgsea", ssgsea.norm=F)
CCLE_GSVA = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva")
CCLE_GSVA_QN = CCLE_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", do_qnorm=T)
CCLE_SING = CCLE_RNA %>% singscore_def(Path_List, filt_genes=T)
CCLE_SING_ST = CCLE_RNA %>% singscore_def(Path_List, filt_genes=T, stableGenes=stable_100)


methods = c("No_Norm", "Z_Sample", "ZRobust_Sample", "ssGSEA", "ssGSEA_Norm_X", 
            "GSVA", "GSVA_QN", "singscore", "stingscore")

Corr_CCLE_Raw = examine_cell_pair(SANGER_RNA[, common_ccle], CCLE_RNA[, common_ccle], methods[1])
Corr_CCLE_ZNorm = examine_cell_pair(SANGER_ZNorm[, common_ccle], CCLE_ZNorm[, common_ccle], methods[2])
Corr_CCLE_ZNormR = examine_cell_pair(SANGER_ZNormR[, common_ccle], CCLE_ZNormR[, common_ccle], methods[3])
Corr_CCLE_ssGSEA = examine_cell_pair(SANGER_ssGSEA, CCLE_ssGSEA, methods[4])
Corr_CCLE_ssGSEA_NN = examine_cell_pair(SANGER_ssGSEA_NN, CCLE_ssGSEA_NN, methods[5])
Corr_CCLE_GSVA = examine_cell_pair(SANGER_GSVA, CCLE_GSVA, methods[6])
Corr_CCLE_GSVA_QN = examine_cell_pair(SANGER_GSVA_QN, CCLE_GSVA_QN, methods[7])
Corr_CCLE_SING = examine_cell_pair(SANGER_SING, CCLE_SING, methods[8])
Corr_CCLE_SING_ST = examine_cell_pair(SANGER_SING_ST, CCLE_SING_ST, methods[9])

Corr_GDSC_Raw = examine_cell_pair(SANGER_RNA[, common_gdsc], GDSC_RNA[, common_gdsc], methods[1])
Corr_GDSC_ZNorm = examine_cell_pair(SANGER_ZNorm[, common_gdsc], GDSC_ZNorm[, common_gdsc], methods[2])
Corr_GDSC_ZNormR = examine_cell_pair(SANGER_ZNormR[, common_gdsc], GDSC_ZNormR[, common_gdsc], methods[3])
Corr_GDSC_ssGSEA = examine_cell_pair(SANGER_ssGSEA, GDSC_ssGSEA, methods[4])
Corr_GDSC_ssGSEA_NN = examine_cell_pair(SANGER_ssGSEA_NN, GDSC_ssGSEA_NN, methods[5])
Corr_GDSC_GSVA = examine_cell_pair(SANGER_GSVA, GDSC_GSVA, methods[6])
Corr_GDSC_GSVA_QN = examine_cell_pair(SANGER_GSVA_QN, GDSC_GSVA_QN, methods[7])
Corr_GDSC_SING = examine_cell_pair(SANGER_SING, GDSC_SING, methods[8])
Corr_GDSC_SING_ST = examine_cell_pair(SANGER_SING_ST, GDSC_SING_ST, methods[9])

col = c("Cell1", "Cell2", "Corr", "Method")

Corr_SANGER_CCLE = Reduce(rbind, 
  list(Corr_CCLE_Raw, Corr_CCLE_ZNorm, Corr_CCLE_ZNormR, Corr_CCLE_ssGSEA, Corr_CCLE_ssGSEA_NN, 
       Corr_CCLE_GSVA, Corr_CCLE_GSVA_QN, Corr_CCLE_SING, Corr_CCLE_SING_ST)
)   # 16871490

Corr_SANGER_GDSC = Reduce(rbind, 
  list(Corr_GDSC_Raw, Corr_GDSC_ZNorm, Corr_GDSC_ZNormR, Corr_GDSC_ssGSEA, Corr_GDSC_ssGSEA_NN, 
       Corr_GDSC_GSVA, Corr_GDSC_GSVA_QN, Corr_GDSC_SING, Corr_GDSC_SING_ST)
)   # 12956274

lvl_cpair = c("Same", "Different")
colnames(Corr_SANGER_GDSC)[1:2] = c("SANGER", "GDSC")
colnames(Corr_SANGER_CCLE)[1:2] = c("SANGER", "CCLE")

# idx = match(Corr_SANGER_CCLE$CCLE, Anno_Cells$BROAD_ID)
# Anno_Cells$BROAD_ID[idx] %>% is.na %>% sum
# Corr_SANGER_CCLE$CCLE_SANGER_ID = Anno_Cells$BROAD_ID[idx]

Corr_SANGER_GDSC = Corr_SANGER_GDSC %>% 
  mutate(SANGER=as.character(SANGER), GDSC=as.character(GDSC),
         Cell_Pair=ifelse(SANGER==GDSC, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=factor(Cell_Pair, levels=lvl_cpair))

Corr_SANGER_CCLE = Corr_SANGER_CCLE %>% 
  mutate(SANGER=as.character(SANGER), CCLE=as.character(CCLE),
         Cell_Pair=ifelse(SANGER==CCLE, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=factor(Cell_Pair, levels=lvl_cpair))

Corr_SANGER_GDSC$Cell_Pair %>% table   # Same 8766 [974*9]
Corr_SANGER_CCLE$Cell_Pair %>% table   # Same 9063 [1007*9]


width = 24
height = 16
legend = "Cell-Pair"

dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Batch Correction")
main = c("SANGER TPM & GDSC Array", "SANGER TPM & CCLE TPM")

file1 = sprintf("%s/Cell-Pair Comparison [%s, PCC]", dir, main)
file2 = sprintf("%s/Cell-Pair Comparison [%s, NRMSE]", dir, main)

pos = position_dodge(width=0.85)
ylab = c("Cell-Pair PCC", "Cell-Pair NRMSE [IQR]")
color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
add = list(scale_fill_manual(values=color))

# margin = margin(b=0.5, l=0.3, unit="cm")
# margin_lg = theme(legend.text=element_text(size=22.5, margin=margin))
# add = list(scale_fill_manual(values=color), margin_lg)

Corr_SANGER_GDSC %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[1], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, alpha_point=0.25, hjust=1, vjust=1, dpi=600, margin_lg=0.5,
              axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
              point=F, violin=F, raster=T, save=T, save_svg=T)

Corr_SANGER_CCLE %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[2], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, alpha_point=0.25, hjust=1, vjust=1, dpi=600, margin_lg=0.5,
              axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
              point=F, violin=F, raster=T, save=T, save_svg=T)

# Corr_SANGER_GDSC %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[1], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=600, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)
# 
# Corr_SANGER_CCLE %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[2], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, text_ratio=1.65, hjust=1, vjust=1, dpi=600, 
#               point=F, violin=F, raster=T, save=T, save_svg=T)

main = c("SANGER TPM & GSDC Array", "SANGER TPM & CCLE TPM")
file = sprintf("%s/Batch Correction [%s].csv", dir, main)
fwrite(Corr_SANGER_GDSC, file=file[1])
fwrite(Corr_SANGER_CCLE, file=file[2])



### 2-3-2. GSVA can alleviate the batch-effects? [O, Standardization]

# Scaler fitted on SANGER RNA-Seq
Scaler_Raw_C = SANGER_RNA[, common_ccle] %>% caret::preProcess()
Scaler_ZNorm_C = SANGER_ZNorm[, common_ccle] %>% caret::preProcess()
Scaler_ZNormR_C = SANGER_ZNormR[, common_ccle] %>% caret::preProcess()

Scaler_Raw_G = SANGER_RNA[, common_gdsc] %>% caret::preProcess()
Scaler_ZNorm_G = SANGER_ZNorm[, common_gdsc] %>% caret::preProcess()
Scaler_ZNormR_G = SANGER_ZNormR[, common_gdsc] %>% caret::preProcess()

Scaler_ssGSEA = SANGER_ssGSEA %>% caret::preProcess()
Scaler_ssGSEA_NN = SANGER_ssGSEA_NN %>% caret::preProcess()
Scaler_GSVA = SANGER_GSVA %>% caret::preProcess()
Scaler_GSVA_QN = SANGER_GSVA_QN %>% caret::preProcess()
Scaler_SING = SANGER_SING %>% caret::preProcess()
Scaler_SING_ST = SANGER_SING_ST %>% caret::preProcess()

# Transform SANGER RNA-Seq
scale_list = function(..., Scaler_List=NULL) {
  
  Omics_List = list(...)
  Omics_List_Scaled = list()
  
  if (length(Omics_List)!=length(Scaler_List)) {
    stop("The number of omics and scalers is not identical...")
  }
  
  for (i in 1:length(Omics_List)) {
    Omics_List_Scaled[[i]] = stats::predict(Scaler_List[[i]], Omics_List[[i]])
  }
  return(Omics_List_Scaled)
}

Scaler_List_G = list(Scaler_Raw_G, Scaler_ZNorm_G, Scaler_ZNormR_G, Scaler_ssGSEA, Scaler_ssGSEA_NN, 
                     Scaler_GSVA, Scaler_GSVA_QN, Scaler_SING, Scaler_SING_ST)

Scaler_List_C = list(Scaler_Raw_C, Scaler_ZNorm_C, Scaler_ZNormR_C, Scaler_ssGSEA, Scaler_ssGSEA_NN, 
                     Scaler_GSVA, Scaler_GSVA_QN, Scaler_SING, Scaler_SING_ST)

SANGER_List_G = scale_list(SANGER_RNA[, common_gdsc], SANGER_ZNorm[, common_gdsc], SANGER_ZNormR[, common_gdsc], 
                           SANGER_ssGSEA, SANGER_ssGSEA_NN, SANGER_GSVA, SANGER_GSVA_QN, SANGER_SING, SANGER_SING_ST, Scaler_List=Scaler_List_G)

SANGER_List_C = scale_list(SANGER_RNA[, common_ccle], SANGER_ZNorm[, common_ccle], SANGER_ZNormR[, common_ccle], 
                           SANGER_ssGSEA, SANGER_ssGSEA_NN, SANGER_GSVA, SANGER_GSVA_QN, SANGER_SING, SANGER_SING_ST, Scaler_List=Scaler_List_C)

GDSC_List = scale_list(GDSC_RNA[, common_gdsc], GDSC_ZNorm[, common_gdsc], GDSC_ZNormR[, common_gdsc], 
                       GDSC_ssGSEA, GDSC_ssGSEA_NN, GDSC_GSVA, GDSC_GSVA_QN, GDSC_SING, GDSC_SING_ST, Scaler_List=Scaler_List_G)

CCLE_List = scale_list(CCLE_RNA[, common_ccle], CCLE_ZNorm[, common_ccle], CCLE_ZNormR[, common_ccle], 
                       CCLE_ssGSEA, CCLE_ssGSEA_NN, CCLE_GSVA, CCLE_GSVA_QN, CCLE_SING, CCLE_SING_ST, Scaler_List=Scaler_List_C)


# Recalculate cell-pair correlation
methods_ = c("Raw", "ZNorm", "ZNormR", "ssGSEA", "ssGSEA_NN", "GSVA", "GSVA_QN", "SING", "SING_ST")
corr_ccle = sprintf("Corr_CCLE_%s_", methods_)
corr_gdsc = sprintf("Corr_GDSC_%s_", methods_)

for (i in 1:length(methods_)) {
  assign(corr_ccle[i], examine_cell_pair(SANGER_List_C[[i]], CCLE_List[[i]], methods[i]))
  assign(corr_gdsc[i], examine_cell_pair(SANGER_List_G[[i]], GDSC_List[[i]], methods[i]))
}

Corr_SANGER_GDSC_Norm = Reduce(rbind, mget(corr_gdsc))
Corr_SANGER_CCLE_Norm = Reduce(rbind, mget(corr_ccle))
colnames(Corr_SANGER_GDSC_Norm)[1:2] = c("SANGER", "GDSC")
colnames(Corr_SANGER_CCLE_Norm)[1:2] = c("SANGER", "CCLE")

# idx = match(Corr_SANGER_CCLE$CCLE, Anno_Cells$BROAD_ID)
# Anno_Cells$BROAD_ID[idx] %>% is.na %>% sum
# Corr_SANGER_CCLE$CCLE_SANGER_ID = Anno_Cells$BROAD_ID[idx]

Corr_SANGER_GDSC_Norm = Corr_SANGER_GDSC_Norm %>% 
  mutate(SANGER=as.character(SANGER), GDSC=as.character(GDSC),
         Cell_Pair=ifelse(SANGER==GDSC, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=factor(Cell_Pair, levels=lvl_cpair))

Corr_SANGER_CCLE_Norm = Corr_SANGER_CCLE_Norm %>% 
  mutate(SANGER=as.character(SANGER), CCLE=as.character(CCLE),
         Cell_Pair=ifelse(SANGER==CCLE, "Same", "Different"), 
         Method=factor(Method, level=methods)) %>% 
  mutate(Cell_Pair=factor(Cell_Pair, levels=lvl_cpair))

Corr_SANGER_CCLE_Norm$Cell_Pair %>% table   # Same 4028 [1007*4]
Corr_SANGER_GDSC_Norm$Cell_Pair %>% table   # Same 3896 [974*4]


dir = mkdir("../../processed_data/cell_data/BIOCARTA/GSVA_Analysis/Batch Correction")
main = c("SANGER TPM & GDSC Array", "SANGER TPM & CCLE TPM")

legend = "Cell-Pair"
file1 = sprintf("%s/Cell-Pair Comparison [%s, PCC (Scaled)]", dir, main)
file2 = sprintf("%s/Cell-Pair Comparison [%s, NRMSE (Scaled)]", dir, main)

pos = position_dodge(width=0.85)
ylab = c("Cell-Pair PCC", "Cell-Pair NRMSE [IQR]")
color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
add = list(scale_fill_manual(values=color))

# margin = margin(b=0.5, l=0.3, unit="cm")
# margin_lg = theme(legend.text=element_text(size=22.5, margin=margin))
# add = list(scale_fill_manual(values=color), margin_lg)

Corr_SANGER_GDSC_Norm %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[1], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, alpha_point=0.25, hjust=1, vjust=1, dpi=600, margin_lg=0.5,
              axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
              point=F, violin=F, raster=T, save=T, save_svg=T)

Corr_SANGER_CCLE_Norm %>% as.data.frame %>% 
  boxplot_def(Method, Corr, fill=Cell_Pair, main=file1[2], add=add, pos=pos, 
              ylab=ylab[1], legend=legend, width=width, height=height, 
              alpha=0.9, alpha_point=0.25, hjust=1, vjust=1, dpi=600, margin_lg=0.5,
              axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
              point=F, violin=F, raster=T, save=T, save_svg=T)

# Corr_SANGER_GDSC_Norm %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[1], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, hjust=1, vjust=1, dpi=600, 
#               axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
#               point=F, violin=F, raster=T, save=T, save_svg=T)
# 
# Corr_SANGER_CCLE_Norm %>% as.data.frame %>% 
#   boxplot_def(Method, NRMSE, fill=Cell_Pair, main=file2[2], add=add, pos=pos, 
#               ylab=ylab[2], legend=legend, width=width, height=height, 
#               alpha=0.9, hjust=1, vjust=1, dpi=600, 
#               axis_tl=30, axis_tx=20, legend_tl=22.5, legend_tx=22.5,
#               point=F, violin=F, raster=T, save=T, save_svg=T)


main = c("SANGER TPM & GSDC Array", "SANGER TPM & CCLE TPM")
file = sprintf("%s/Batch Correction [%s (Scaled)].csv", dir, main)
fwrite(Corr_SANGER_GDSC_Norm, file=file[1])
fwrite(Corr_SANGER_CCLE_Norm, file=file[2])



### 2-4. Complete Pathway-Pathway correlation network [RNA, k=5]

knn_graph = function(Network, k=5, mode="min", col_stat="Dist", filt_na=F) {
  
  # Complete KNN graph from Network
  # Network should have the following columns
  # Pathway1, Pathway2, col_stat [Dist or Corr]
  
  Network_Rev = Network %>% dplyr::rename(Pathway1=Pathway2, Pathway2=Pathway1)
  Network = Network %>% rbind(Network_Rev) %>% distinct(Pathway1, Pathway2, .keep_all=T)
  
  filter_net = ifelse(mode=="max", slice_max, slice_min)
  if (filt_na) Network = Network %>% subset(!is.na(object(col_stat)))
  
  Network_KNN = Network %>% group_by(Pathway1) %>% 
    filter_net(object(col_stat), n=k) %>% as.data.frame   # n x (n-1) > n x k
  
  return(Network_KNN)
}

RNA_Corr = SANGER_GSVA %>% cor %>% reshape2::melt() %>% as.data.frame
colnames(RNA_Corr) = c("Pathway1", "Pathway2", "Corr")
RNA_Corr = RNA_Corr %>% subset(Pathway1!=Pathway2)   # 36672 [192 x 191]

RNA_Corr_KNN3 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=3)   # 576
RNA_Corr_KNN5 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=5)   # 960
RNA_Corr_KNN7 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=7)   # 1344

# RNA_Corr$Corr %>% hist
# RNA_Corr_KNN3$Corr %>% hist
# RNA_Corr_KNN5$Corr %>% hist


# Save files
main = c("SANGER", "GDSC", "CCLE")
dir = mkdir("../../processed_data/cell_data/BIOCARTA")

file = sprintf("%s/%s_RNA_GSVA.csv", dir, main)
write.csv(SANGER_GSVA, file=file[1], row.names=T)
write.csv(GDSC_RNA_GSVA, file=file[2], row.names=T)
write.csv(CCLE_RNA_GSVA, file=file[3], row.names=T)

file = sprintf("%s/%s_RNA_GSVA_QN.csv", dir, main)
write.csv(SANGER_GSVA_QN, file=file[1], row.names=T)
write.csv(GDSC_RNA_GSVA_QN, file=file[2], row.names=T)
write.csv(CCLE_RNA_GSVA_QN, file=file[3], row.names=T)

file = sprintf("%s/%s_RNA_ssGSEA.csv", dir, main)
write.csv(SANGER_ssGSEA, file=file[1], row.names=T)
write.csv(GDSC_RNA_ssGSEA, file=file[2], row.names=T)
write.csv(CCLE_RNA_ssGSEA, file=file[3], row.names=T)

file = sprintf("%s/%s_RNA_SING.csv", dir, main)
write.csv(SANGER_SING, file=file[1], row.names=T)
write.csv(GDSC_RNA_SING, file=file[2], row.names=T)
write.csv(CCLE_RNA_SING, file=file[3], row.names=T)

dir = mkdir("../../processed_data/net_data/BIOCARTA")
file = sprintf("%s/RNA_Corr_KNN%s.csv", dir, c(3, 5, 7))
write.csv(RNA_Corr_KNN3, file=file[1], row.names=F)
write.csv(RNA_Corr_KNN5, file=file[2], row.names=F)
write.csv(RNA_Corr_KNN7, file=file[3], row.names=F)

file = sprintf("%s/RNA_Corr.csv", dir)
write.csv(RNA_Corr, file=file, row.names=F)



### 2-5. Compress STRING & RegNetwork into Pathway-Pathway Association network [k=5]
# [2015] Uncovering disease-disease relationships through the incomplete interactome

calc_sep_score = function(Network, geneset_A, geneset_B) {
  
  # Calculate separation score of a pair of pathways
  # This function requires the installation of dnet, igraph package
  # Users must use the function separation_score, rather than this function
  
  inf_omit = function(x) x[!is.infinite(x)]
  dist_AA = distances(Network, geneset_A, geneset_A, mode="out") %>% as.numeric %>% inf_omit %>% mean
  dist_BB = distances(Network, geneset_B, geneset_B, mode="out") %>% as.numeric %>% inf_omit %>% mean
  dist_AB = distances(Network, geneset_A, geneset_B, mode="out") %>% as.numeric %>% inf_omit %>% mean
  
  if (dist_AA==0 | is.na(dist_AA)) dist_AA = NA
  if (dist_BB==0 | is.na(dist_BB)) dist_BB = NA
  if (is.nan(dist_AB)) dist_AB = NA
  
  # Separation score is only calculated when all dists are calculated
  # If geneset_A are not connected, dist_AA = 0 > NA
  # If geneset_B are not connected, dist_BB = 0 > NA
  # If geneset_AB are not connected, dist_AB = NaN > NA
  
  ss_AB = dist_AB - (dist_AA + dist_BB) / 2
  ss_info = c(ss_AB, dist_AA, dist_BB, dist_AB)
  return(ss_info)
}

separation_score = function(Network, Path_List, min_genes=NULL, verbose=T, cores=10) {
  
  # Calculate separation score of all pairwise pathways
  # This function requires the installation of dnet, igraph package
  # The Network should contain two columns named Node1, Node2
  # Otherwise, the first two columns are recognized as Nodes
  
  suppressMessages(library(dnet))
  cond1 = "Node1" %in% colnames(Network)
  cond2 = "Node2" %in% colnames(Network)
  
  if (!cond1 | !cond2) {
    Network = Network[, 1:2]
    colnames(Network) = c("Node1", "Node2")
  }
  
  genes_net = union(Network$Node1, Network$Node2)
  genes_path = Path_List %>% unlist %>% unique
  genes_int = intersect(genes_net, genes_path)
  
  Path_List_Ori = Path_List
  int_ratio = length(genes_int) / length(genes_path)
  Path_List = Path_List %>% sapply(function(x) x[x %in% genes_int])
  
  sprintf("# Pathway size range : [%s, %s]", 
          min(sapply(Path_List, length)), max(sapply(Path_List, length))) %>% print
  
  sprintf("# Total Genes utilized : %.2f%% [%s/%s]", 
          int_ratio*100, length(genes_int), length(genes_path)) %>% print
  
  if (is.numeric(min_genes)) {
    idx = sapply(Path_List, length)>=min_genes
    Path_List = Path_List[idx]
    sprintf("# Total Pathway utilized : %.2f%% [%s/%s]", 
            100 * sum(idx) / length(idx), sum(idx), length(idx)) %>% print
  }
  
  Network = Network %>% dplyr::select(from=Node1, to=Node2)
  Network = Network %>% graph_from_data_frame(directed=T)
  nC2 = combn(length(Path_List), 2) %>% t   
  # [n_path x (n_path-1)/2, 2]
  
  if (cores<2) {
    Network_SS = data.frame()
    for (i in 1:nrow(nC2)) {
      path_A = names(Path_List)[nC2[i, 1]]
      path_B = names(Path_List)[nC2[i, 2]]
      geneset_A = Path_List[[path_A]]
      geneset_B = Path_List[[path_B]]
      
      ss_info = calc_sep_score(Network, geneset_A, geneset_B)
      Network_SS_ = c(path_A, path_B, ss_info)
      Network_SS = Network_SS %>% rbind(Network_SS_)
    }
  } else {
    cluster = multicores(cores=cores)
    Network_SS = foreach (i=1:nrow(nC2), .combine=rbind) %dopar% {
      path_A = names(Path_List)[nC2[i, 1]]
      path_B = names(Path_List)[nC2[i, 2]]
      geneset_A = Path_List[[path_A]]
      geneset_B = Path_List[[path_B]]
      
      ss_info = calc_sep_score(Network, geneset_A, geneset_B)
      Network_SS_ = c(path_A, path_B, ss_info)
      return(Network_SS_)
    }
    Network_SS = Network_SS %>% as.data.frame
    stopCluster(cluster)
  }
  
  rownames(Network_SS) = NULL
  Network_SS[, 3:6] = Network_SS[, 3:6] %>% sapply(as.numeric)
  colnames(Network_SS) = c("Pathway1", "Pathway2", "Dist", "Dist_P1", "Dist_P2", "Dist_P1_P2")
  
  inlen = function(x, y) intersect(x, y) %>% length
  num_overlap = nC2 %>% apply(1, function(x) Reduce(inlen, Path_List_Ori[x]))
  Network_SS$Num_Overlap = num_overlap
  
  Path_Nums = sapply(Path_List_Ori, length)
  idx1 = match(Network_SS$Pathway1, names(Path_Nums))
  idx2 = match(Network_SS$Pathway2, names(Path_Nums))
  
  Network_SS = Network_SS %>% 
    mutate(Num1 = Path_Nums[idx1], Num2 = Path_Nums[idx2],
           Ratio_Overlap = Num_Overlap / (Num1+Num2-Num_Overlap)) %>% 
    relocate(Num_Overlap, .before=Ratio_Overlap)
  
  # Undirect the Pathway-Pathway Interaction Graph
  Network_SS_Rev = Network_SS
  colnames(Network_SS_Rev)[1:2] = c("Pathway2", "Pathway1")
  Network_SS = Network_SS %>% rbind(Network_SS_Rev)
  
  return(Network_SS)
}

ggnet2_def = function(Network, col_from=NULL, col_to=NULL, main=NULL, nodes=NULL, add=NULL,
                      color="steelblue1", color_edge="grey50", size=4, size_edge=0.4, 
                      alpha=1, alpha_edge=0.8, width=18, height=15, dpi=400,
                      label=NULL, input_net=F, verbose=T, save=F) {
  
  # suppressMessages(library(ggnetwork))
  # suppressMessages(library(GGally))
  # suppressMessages(library(intergraph))
  
  if (!input_net) {
    if (!is.null(col_from) & !is.null(col_to)) {
      col_from = deparse(substitute(col_from))
      col_to = deparse(substitute(col_to))
      Network = Network[, c(col_from, col_to)]
    } else Network = Network[, 1:2]
    
    # Funciton graph_from_data_frame automatically takes edges from first 2 columns
    # Function ggnet2 cannot plot multiple node-node relations [should set directed=T]
    Network = Network %>% graph_from_data_frame(directed=T)
  }
  
  if (verbose) {
    net_info = igraph::components(Network)
    sprintf("# Nodes Number : %s", length(V(Network))) %>% print
    sprintf("# Edges Number : %s", length(E(Network))) %>% print
    sprintf("# Subgraph Number : %s", net_info$no) %>% print
    for (i in 1:length(net_info$csize)) {
      sprintf("# Subgraph Size [%sth] : %s", i, net_info$csize[i]) %>% print
    }
  }
  
  V(Network)$color = color
  if (!is.null(nodes)) {
    Network = delete_vertices(Network, !(names(V(Network)) %in% nodes))
  }
  
  Network = ggnetwork(Network)
  if (!is.null(label)) {
    Network$name = ifelse(Network$name %in% label, Network$name, NA)
  }
  
  pl = ggplot(Network, aes(x=x, y=y, xend=xend, yend=yend)) + 
    theme_blank() + 
    geom_edges(color=color_edge, alpha=alpha_edge, size=size_edge) + 
    geom_nodes(color=color, alpha=alpha, size=size)
  
  # geom_nodetext_repel(aes(label=name), max.overlaps=10)
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  
  if (save) {
    save_fig(pl, main=main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
}

knn_graph_ovl = function(Network_SS, k=5, inverse=F) {
  
  # Complete a KNN graph from the separation score
  # Network_SS should have the following columns
  # [Pathway1, Pathway2, Dist, Ratio_Overlap]
  
  # If tied with or do not have separation scores [Dist],
  # Filter the pairs by pathway overlap ratios [Ratio_Overlap]
  # Some pathways might not have KNN by utilizing either Dist or Ratio_Overlap
  
  if (!inverse) {
    # Network_KNN = Network_SS %>% mutate(Pathway=as.factor(Pathway1)) %>%
    #   group_by(Pathway) %>% arrange(Dist, desc(Ratio_Overlap), .by_group=T) %>%
    #   do(head(., n=k)) %>% subset(!is.na(Dist) | !Ratio_Overlap==0) %>% as.data.frame
    
    Network_KNN = Network_SS %>% mutate(Pathway=as.factor(Pathway1)) %>%
      group_by(Pathway) %>% arrange(Dist, desc(Ratio_Overlap), .by_group=T) %>%
      do(head(., n=k)) %>% as.data.frame
  } else {
    Network_KNN = Network_SS %>% mutate(Pathway=as.factor(Pathway1)) %>% 
      group_by(Pathway) %>% arrange(desc(Dist), Ratio_Overlap, .by_group=T) %>% 
      do(head(., n=k)) %>% as.data.frame
  }
  
  info_knn = Network_KNN$Pathway %>% table
  idx = which(info_knn != k)
  
  if (length(idx)!=0) {
    for (i in idx) {
      path = names(info_knn)[i]
      num_nn = info_knn[i]
      sprintf("# Pathway with lower than k-edges : %s [%s]", path, num_nn) %>% print
    }
  }
  
  num_full = length(info_knn) * k
  sprintf("# The num of KNN edges : %s / %s", nrow(Network_KNN), num_full) %>% print
  return(Network_KNN)
}

knn_graph_rand = function(Path_List, k=5, seed=2021) {
  
  set.seed(seed)
  Network_KNN = data.frame()
  pathways = Path_List %>% names
  
  for (i in 1:length(pathways)) {
    path_query = rep(pathways[i], k)
    path_rand = sample(pathways[-i], k)
    Network_KNN_ = data.frame(Pathway1=path_query, Pathway2=path_rand)
    Network_KNN = Network_KNN %>% rbind(Network_KNN_)
  }
  
  sprintf("KNN Edges : %s", nrow(Network_KNN)) %>% print
  return(Network_KNN)
}

row_overlap = function(df1, df2, col=c("Pathway1", "Pathway2"), verbose=T) {
  
  nrow_min = min(nrow(df1), nrow(df2))
  if (is.null(col)) col = intersect(colnames(df1), colnames(df2))
  nrow_overlap = generics::intersect(df1[, col], df2[, col]) %>% nrow
  per_overlap = 100 * nrow_overlap / nrow_min
  
  if (verbose) {
    sprintf("Row Intersect : %s / %s [%.2f%%]", nrow_overlap, nrow_min, per_overlap) %>% print
  } else return(per_overlap)
}

row_overlap_heatmap = function(Net_List, col=c("Pathway1", "Pathway2"), 
                               net_names=NULL, main=NULL, size=6.4, round=2, 
                               width=22.5, height=16, dpi=400, save=T) {
  
  # Net_List = list(...)
  if (is.null(net_names)) net_names = sprintf("Network%s", 1:length(Net_List))
  
  row_ovl = c()
  names(Net_List) = net_names
  Net_Overlap_Info = expand.grid(Net1=net_names, Net2=net_names)
  
  for (i in 1:nrow(Net_Overlap_Info)) {
    net1 = Net_Overlap_Info$Net1[i]
    net2 = Net_Overlap_Info$Net2[i]
    row_ovl_temp = row_overlap(Net_List[[net1]], Net_List[[net2]], col=col, verbose=F)
    row_ovl = row_ovl %>% c(row_ovl_temp)
  }
  
  legend = "Similarity [%]"
  Net_Overlap_Info$Similarity = row_ovl
  color = scale_fill_gradient(low="white", high="firebrick3")
  add = list(theme(legend.key.width=unit(0.8, "cm"), 
                   legend.key.height=unit(1, "cm")))
  
  Net_Overlap_Info %>% 
    grid_def(Net1, Net2, fill=Similarity, main=main, color=color,
             mean_summary=F, legend=legend, dpi=dpi, size=size, add=add,
             axis_tx=18, legend_tl=18, legend_tx=18, margin_lg=0.75, round=round, 
             axis_face="plain", legend_face="plain", width=width, height=height, save=save)
  
  return(Net_Overlap_Info)
}

dir = "../../processed_data/net_data/STRING"
file= sprintf("%s/STRING_Filt_Ent.csv", dir)
STRING_700 = read.csv(file)

dir = "../../processed_data/net_data/RegNetwork"
file= sprintf("%s/RegNetwork_Ent.csv", dir)
RegNetwork = read.csv(file)

STRING_700 = STRING_700 %>% mutate(Node1=as.character(Node1), 
                                   Node2=as.character(Node2))
RegNetwork = RegNetwork %>% mutate(Node1=as.character(Node1), 
                                   Node2=as.character(Node2))
STRING_900 = STRING_700 %>% subset(Combined_Score>=900)   # 479908 > 236392

col = c("Node1", "Node2")
row_overlap(STRING_700, RegNetwork, col=col)   # 10036 / 371910 [2.70%]

STRING_700_SS = STRING_700 %>% separation_score(Path_List)
# Pathway size range : [10, 351]
# Total genes utilized : 92.06% [4939/5365]

STRING_900_SS = STRING_900 %>% separation_score(Path_List)
# Pathway size range : [9, 345]
# Total Genes utilized : 85.83% [4605/5365]

RegNetwork_SS = RegNetwork %>% separation_score(Path_List)
# Pathway size range : [10, 349]
# Total Genes utilized : 94.28% [5058/5365]


# [Confirmation] Do separation scores truly represent the distance of pathway-pathway?
net_type = c("STRING_700", "STRING_900", "RegNetwork")
dir = mkdir("../../processed_data/net_data/BIOCARTA/Network_Analysis")
main = sprintf("%s/Separation Score & Overlap Ratio [%s]", dir, net_type)

xlab = "Separation Score"
ylab = "Overlap Ratio"

# STRING_700_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[1], xlab=xlab, ylab=ylab, 
#            alpha=0.5, axis_tl=20, save=T)   # [△]
# STRING_900_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[2], xlab=xlab, ylab=ylab, 
#            alpha=0.5, axis_tl=20, save=T)   # [△]
# RegNetwork_SS %>% 
#   plot_def(Dist, Ratio_Overlap, main=main[3], xlab=xlab, ylab=ylab, 
#            alpha=0.5,axis_tl=20, save=T)   # [O]

col = list(fill="lightgray")
margin = margin(10, 10, 10, 10, unit="pt")
font_label = font("xylab", size=30, margin=margin)
font_text = font("xy.text", size=22.5, color="grey30", margin=margin)

font_ = font_label + font_text
main = sprintf("%s/Separation Score & Overlap Ratio [%s, with histogram]", dir, net_type)

pl1 <- STRING_700_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)
pl2 = STRING_900_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)
pl3 = RegNetwork_SS %>% 
  ggscatterhist("Dist", "Ratio_Overlap", xlab=xlab, ylab=ylab, 
                alpha=0.25, size=1.5, margin.params=col)

pl1$sp = pl1$sp + font_
pl2$sp = pl2$sp + font_
pl3$sp = pl3$sp + font_

suppressMessages(library(ggrastr))
pl1$sp$layers[[1]] = pl1$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")
pl2$sp$layers[[1]] = pl2$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")
pl3$sp$layers[[1]] = pl3$sp$layers[[1]] %>% rasterise(dpi=1200, dev="ragg_png")

pl1 %>% save_fig_ggpubr(main=main[1], width=15, height=15, dpi=1200, svg=T)
pl2 %>% save_fig_ggpubr(main=main[2], width=15, height=15, dpi=1200, svg=T)
pl3 %>% save_fig_ggpubr(main=main[3], width=15, height=15, dpi=1200, svg=T)

# Some pathways in RegNetwork have negative separative scores with no gene-overlap ratios...
STRING_700_SS$Dist %>% na.omit %>% range %>% round(3)   # 0.005   2.630
STRING_900_SS$Dist %>% na.omit %>% range %>% round(3)   # 0.006   3.164
RegNetwork_SS$Dist %>% na.omit %>% range %>% round(3)   # -0.342  3.017

# RegNetwork_SS %>% subset(Dist<0) %>% 
#   plot_def(Dist, Dist_P1_P2, color=Ratio_Overlap, size=Ratio_Overlap, alpha=0.5, save=F)
# 
# RegNetwork_SS %>% subset(Dist<0) %>% 
#   plot_def(Dist, Ratio_Overlap, color=Dist_P1_P2, size=Dist_P1_P2, alpha=0.25, save=F)
# 
# RegNetwork_SS %>% subset(Dist<0) %>% 
#   plot_def(Dist_P1_P2, Ratio_Overlap, color=Dist, size=Dist, alpha=0.25, save=F)
# 
# RegNetwork_SS %>% subset(Dist<0) %>% 
#   plot_def(Dist_P1_P2, Dist_P1+Dist_P2, color=Dist, size=Ratio_Overlap, alpha=0.25, save=F)


# KNN Graphs
STRING_700_KNN3 = STRING_700_SS %>% knn_graph_ovl(k=3)   # 876 [292*3]
STRING_900_KNN3 = STRING_900_SS %>% knn_graph_ovl(k=3)   # 876 [292*3]
RegNetwork_KNN3 = RegNetwork_SS %>% knn_graph_ovl(k=3)   # 876 [292*3]

STRING_700_KNN5 = STRING_700_SS %>% knn_graph_ovl(k=5)   # 1460 [292*5]
STRING_900_KNN5 = STRING_900_SS %>% knn_graph_ovl(k=5)   # 1460 [292*5]
RegNetwork_KNN5 = RegNetwork_SS %>% knn_graph_ovl(k=5)   # 1460 [292*5]

STRING_700_KNN7 = STRING_700_SS %>% knn_graph_ovl(k=7)   # 2044 [292*7]
STRING_900_KNN7 = STRING_900_SS %>% knn_graph_ovl(k=7)   # 2044 [292*7]
RegNetwork_KNN7 = RegNetwork_SS %>% knn_graph_ovl(k=7)   # 2044 [292*7]

# Cf. KNN Graphs with Maximum Separation Scores [Ablation Test]
STRING_700_KNN5_Inv = STRING_700_SS %>% knn_graph_ovl(k=5, inverse=T)
STRING_900_KNN5_Inv = STRING_900_SS %>% knn_graph_ovl(k=5, inverse=T)
RegNetwork_KNN5_Inv = RegNetwork_SS %>% knn_graph_ovl(k=5, inverse=T)


# KNN Graphs [Random]
seed = 2021
Random_KNN3 = knn_graph_rand(Path_List, k=3, seed=seed)
Random_KNN5 = knn_graph_rand(Path_List, k=5, seed=seed)
Random_KNN7 = knn_graph_rand(Path_List, k=7, seed=seed)

# Random_KNN3_ = knn_graph_rand(Path_List, k=3, seed=seed)
# Random_KNN5_ = knn_graph_rand(Path_List, k=5, seed=seed)
# Random_KNN7_ = knn_graph_rand(Path_List, k=7, seed=seed)
# 
# identical(Random_KNN3, Random_KNN3_)   # T
# identical(Random_KNN5, Random_KNN5_)   # T
# identical(Random_KNN7, Random_KNN7_)   # T

min_max_mean = function(x) {
  c(mean(x, na.rm=T), min(x, na.rm=T), max(x, na.rm=T))
}

examine_graph = function(Network, col=c("Pathway1", "Pathway2")) {
  
  num_edge = Network %>% nrow
  num_node = union(Network[, col[1]], Network[, col[2]]) %>% length
  
  Network = Network[, col] %>% graph_from_data_frame(directed=T)
  degree_list = Network %>% igraph::degree()
  path_list = Network %>% igraph::distances() %>% as.data.frame %>% unlist
  
  num_subgraph = components(Network)$no
  info_degree = degree_list %>% min_max_mean
  info_path = path_list[path_list!=0] %>% min_max_mean
  
  info_graph = c(num_node, num_edge, info_degree, info_path, num_subgraph)
  return(info_graph)
}

examine_graph_list = function(..., net_names=NULL, col=c("Pathway1", "Pathway2")) {
  
  Network_Info = data.frame()
  for (Network in list(...)) {
    network_info = Network %>% examine_graph(col=col)
    Network_Info = Network_Info %>% rbind(network_info)
  }
  
  col_ = c("Mean", "Min", "Max")
  col1 = sprintf("Degree_%s", col_)
  col2 = sprintf("Path_%s", col_)
  
  col_info = c("Nodes", "Edges", col1, col2, "Num_Component")
  colnames(Network_Info) = col_info
  Network_Info = Network_Info %>% mutate(Network=net_names) %>% 
    relocate(Network, .before=everything()) %>% as.data.frame
  
  return(Network_Info)
}

choose_freq_path = function(Network, col=c("Pathway1", "Pathway2"), top=10) {
  col_rev = col %>% setNames(rev(col))
  Network_Rev = Network %>% rename(col_rev)
  Network_UD = Network[, col] %>% rbind(Network_Rev[, col]) %>% distinct
  path_degree = Network_UD[[col[1]]] %>% table %>% sort(decreasing=T)
  pathways = names(path_degree)[1:top]
  return(pathways)
}

net_names = c("STRING_700", "STRING_900", "RegNetwork", "RNA_Corr", "Random")
PPA_Info = examine_graph_list(STRING_700_KNN5, STRING_900_KNN5, RegNetwork_KNN5, 
                              RNA_Corr_KNN5, Random_KNN5, net_names=net_names)


# All graphs have no multiple subgraphs
suppressMessages(library(ggnetwork))
main = c("GSVA_Corr", "RegNetwork", "STRING_700", "STRING_900")
dir = mkdir("../../processed_data/net_data/BIOCARTA/Network_Analysis")
file = sprintf("%s/KNN5 Graph [%s]", dir, main)

label_corr = RNA_Corr_KNN5 %>% choose_freq_path(top=20)
label_regnet = RegNetwork_KNN5 %>% choose_freq_path(top=20)
label_str_700 = STRING_700_KNN5 %>% choose_freq_path(top=20)
label_str_900 = STRING_900_KNN5 %>% choose_freq_path(top=20)

color = "royalblue1"
add = list(geom_nodetext_repel(aes(label=name), max.overlaps=8, force=2, size=3.2))
RNA_Corr_KNN5 %>% ggnet2_def(main=file[1], color=color, label=label_corr, add=add, dpi=1500, save=T)
RegNetwork_KNN5 %>% ggnet2_def(main=file[2], color=color, label=label_regnet, add=add, dpi=1500, save=T)
STRING_700_KNN5 %>% ggnet2_def(main=file[3], color=color, label=label_str_700, add=add, dpi=1500, save=T)
STRING_900_KNN5 %>% ggnet2_def(main=file[4], color=color, label=label_str_900, add=add, dpi=1500, save=T)

# Similarity between KNN Networks [a little]
row_overlap(STRING_700_KNN5, STRING_900_KNN5)   # 1108 / 1460 [75.89%]
row_overlap(STRING_700_KNN5, RegNetwork_KNN5)   # 315 / 1460 [21.58%]
row_overlap(STRING_700_KNN5, RNA_Corr_KNN5)     # 552 / 1460 [37.81%]

net_names = c("Random", "GSVA_Corr", "RegNetwork", "STRING_700", "STRING_900")
KNN3_List = list(Random_KNN3, RNA_Corr_KNN3, RegNetwork_KNN3, STRING_700_KNN3, STRING_900_KNN3)
KNN5_List = list(Random_KNN5, RNA_Corr_KNN5, RegNetwork_KNN5, STRING_700_KNN5, STRING_900_KNN5)
KNN7_List = list(Random_KNN7, RNA_Corr_KNN7, RegNetwork_KNN7, STRING_700_KNN7, STRING_900_KNN7)

knn = c("KNN3", "KNN5", "KNN7")
dir = "../../processed_data/net_data/BIOCARTA/Network_Analysis"
main = sprintf("%s/KNN Graph Similarity [%s]", dir, knn)

KNN3_Sim = row_overlap_heatmap(KNN3_List, net_names=net_names, main=main[1], save=T)
KNN5_Sim = row_overlap_heatmap(KNN5_List, net_names=net_names, main=main[2], save=T)
KNN7_Sim = row_overlap_heatmap(KNN7_List, net_names=net_names, main=main[3], save=T)

file = sprintf("%s/PPA_Information.csv", dir)
write.csv(PPA_Info, file=file, row.names=F)



##### 4. Save Files

dir = "../../processed_data/net_data/BIOCARTA"

main = c("STRING_700", "STRING_900", "RegNetwork")
file = sprintf("%s/SS_%s.csv", dir, main)
write.csv(STRING_700_SS, file=file[1], row.names=F)
write.csv(STRING_900_SS, file=file[2], row.names=F)
write.csv(RegNetwork_SS, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STRING_700.csv", dir, c(3, 5, 7))
write.csv(STRING_700_KNN3, file=file[1], row.names=F)
write.csv(STRING_700_KNN5, file=file[2], row.names=F)
write.csv(STRING_700_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_STRING_900.csv", dir, c(3, 5, 7))
write.csv(STRING_900_KNN3, file=file[1], row.names=F)
write.csv(STRING_900_KNN5, file=file[2], row.names=F)
write.csv(STRING_900_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_RegNetwork.csv", dir, c(3, 5, 7))
write.csv(RegNetwork_KNN3, file=file[1], row.names=F)
write.csv(RegNetwork_KNN5, file=file[2], row.names=F)
write.csv(RegNetwork_KNN7, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_Random.csv", dir, c(3, 5, 7))
write.csv(Random_KNN3, file=file[1], row.names=F)
write.csv(Random_KNN5, file=file[2], row.names=F)
write.csv(Random_KNN7, file=file[3], row.names=F)

main = c("STRING_700", "STRING_900", "RegNetwork")
file = sprintf("%s/KNN5_%s_Inv.csv", dir, main)
write.csv(STRING_700_KNN5_Inv, file=file[1], row.names=F)
write.csv(STRING_900_KNN5_Inv, file=file[2], row.names=F)
write.csv(RegNetwork_KNN5_Inv, file=file[3], row.names=F)

# Merge three graphs for RGCN

col = c("Pathway1", "Pathway2")
Edge_KNN3 = Reduce(rbind, list(STRING_900_KNN3[, col] %>% mutate(Edge_Type="STRING_900"), 
                               RegNetwork_KNN3[, col] %>% mutate(Edge_Type="RegNetwork"), 
                               RNA_Corr_KNN3[, col] %>% mutate(Edge_Type="RNA_Corr"))) 

Edge_KNN5 = Reduce(rbind, list(STRING_900_KNN5[, col] %>% mutate(Edge_Type="STRING_900"), 
                               RegNetwork_KNN5[, col] %>% mutate(Edge_Type="RegNetwork"), 
                               RNA_Corr_KNN5[, col] %>% mutate(Edge_Type="RNA_Corr"))) 

Edge_KNN7 = Reduce(rbind, list(STRING_900_KNN7[, col] %>% mutate(Edge_Type="STRING_900"), 
                               RegNetwork_KNN7[, col] %>% mutate(Edge_Type="RegNetwork"), 
                               RNA_Corr_KNN7[, col] %>% mutate(Edge_Type="RNA_Corr"))) 

RNA_Corr_KNN5_Inv = RNA_Corr %>% knn_graph(col_stat="Corr", mode="min", k=5)   # 960
Edge_Inv = Reduce(rbind, list(STRING_900_KNN5_Inv[, col] %>% mutate(Edge_Type="STRING_900"), 
                              RegNetwork_KNN5_Inv[, col] %>% mutate(Edge_Type="RegNetwork"), 
                              RNA_Corr_KNN5_Inv[, col] %>% mutate(Edge_Type="RNA_Corr")))

file = sprintf("%s/KNN5_RNA_Corr_Inv.csv", dir)
write.csv(RNA_Corr_KNN5_Inv, file=file, row.names=F)

file = sprintf("%s/KNN3_STR9_Reg_Corr.csv", dir)
write.csv(Edge_KNN3, file=file, row.names=F)

file = sprintf("%s/KNN5_STR9_Reg_Corr.csv", dir)
write.csv(Edge_KNN5, file=file, row.names=F)

file = sprintf("%s/KNN7_STR9_Reg_Corr.csv", dir)
write.csv(Edge_KNN7, file=file, row.names=F)

file = sprintf("%s/KNN5_STR9_Reg_Corr_Inv.csv", dir)
write.csv(Edge_Inv, file=file, row.names=F)


file = sprintf("%s/Pathway_Analysis.RData", dir)
save(RegNetwork_SS, STRING_700_SS, STRING_900_SS, RNA_Corr,
     Random_KNN3, Random_KNN5, Random_KNN7,
     STRING_700_KNN3, STRING_900_KNN3, RegNetwork_KNN3, RNA_Corr_KNN3, 
     STRING_700_KNN5, STRING_900_KNN5, RegNetwork_KNN5, RNA_Corr_KNN5, 
     STRING_700_KNN7, STRING_900_KNN7, RegNetwork_KNN7, RNA_Corr_KNN7, file=file)


supplementary = T
if (supplementary) {
  
  save_for_nc = function(df_list, dir=".", num=1, num_fig=NULL, rowNames=F, suppl=T) {
    
    suppressMessages(library(openxlsx))
    is_list = inherits(df_list, "list")
    if (is_list & is.null(num_fig)) num_fig = letters[1:length(df_list)]
    
    if (!suppl) {
      sheets = sprintf("Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Fig. %s", num)
      file = sprintf("%s/SourceData_Fig%s.xlsx", dir, num)
    } else {
      sheets = sprintf("Supplementary Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Supplementary Fig. %s", num)
      file = sprintf("%s/SourceData_SupplementaryFig%s.xlsx", dir, num)
    }
    
    write.xlsx(df_list, file=file, sheetName=sheets, rowNames=rowNames)
  }
  
  
  ### [Supplementar Data] Supplementar Data 2
  Edge_KNN5_ = Edge_KNN5
  Edge_KNN5_$Edge_Type[Edge_KNN5_$Edge_Type=="RNA_Corr"] = "GSVA_Corr"
  save_for_nc(Edge_KNN5_, num=2, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 1
  STRING_700_SS_ = STRING_700_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  STRING_900_SS_ = STRING_900_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  RegNetwork_SS_ = RegNetwork_SS %>% 
    rename(Separation_Score=Dist, Dist_Pathway1=Dist_P1, Dist_Pathway2=Dist_P2, 
           Dist_Pathway12=Dist_P1_P2, Num_Genes_Pathway1=Num1, Num_Genes_Pathway2=Num2, 
           Num_Genes_Overlap=Num_Overlap, Ratio_Genes_Overlap=Ratio_Overlap)
  
  KNN5_Sim_ = KNN5_Sim %>% rename(Network_X=Net1, Network_Y=Net2)
  PPA_Info_ = list(STRING_700_SS_, STRING_900_SS_, RegNetwork_SS_, KNN5_Sim_)
  
  # dir = "../../processed_data/net_data/BIOCARTA"
  PPA_Info_ %>% save_for_nc(num=1, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 2
  STRING_700_KNN5_ = STRING_700_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  STRING_900_KNN5_ = STRING_900_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  RegNetwork_KNN5_ = RegNetwork_KNN5 %>% 
    subset(select=c(Pathway1, Pathway2, Dist, Ratio_Overlap)) %>% 
    rename(Separation_Score=Dist, Ratio_Genes_Overlap=Ratio_Overlap)
  
  KNN5_List_ = list(STRING_700_KNN5_, STRING_900_KNN5_, RegNetwork_KNN5_, RNA_Corr_KNN5)
  
  # dir = "../../processed_data/net_data/BIOCARTA"
  KNN5_List_ %>% save_for_nc(num=2, suppl=T)
}
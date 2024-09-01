#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
source("../functions.R")
loadings()

file = "2128_genes.txt"
geneset_ori = scan(file, what="")

dir = "../../do_better/PaccMann/data/gene_expression"
file = sprintf("%s/gdsc-rnaseq_gene-expression.csv", dir)
EXP_PaccMann_Ori = fread_def(file)   # 457 x 19689
cells_ori = EXP_PaccMann_Ori %>% rownames

### Not all genes are included in RNA-Seq data
# https://github.com/PaccMann/paccmann_rl/issues/9
# "Regarding missing genes, we do zero-imputation..."

dir = "../../raw_data/GDSC"
file1 = sprintf("%s/E-MTAB-3983-query-results.tpms.tsv", dir)
file2 = sprintf("%s/E-MTAB-3983-query-results.fpkms.tsv", dir)
EXP_TPM = fread(file1)    # 57503 x 459
EXP_FPKM = fread(file2)   # 54352 x 459

# Both data are raw counts, not Log2-scale
EXP_TPM_Ori[, -c(1:2)] %>% colSums(na.rm=T) %>% range    # 999762.3 1000066.8
EXP_FPKM_Ori[, -c(1:2)] %>% colSums(na.rm=T) %>% range   # 238894.6 2013142.7

names(EXP_TPM)[1:2] = c("Gene_ID", "Gene_Name")
names(EXP_FPKM)[1:2] = c("Gene_ID", "Gene_Name")

# All cell-lines in original data are already in processed data in PaccMann
# Not need to process the raw data from the start...
identical(colnames(EXP_TPM), colnames(EXP_FPKM))
cells = colnames(EXP_TPM)[-c(1:2)]
cells = cells %>% strsplit(",") %>% sapply(function(x) x[1])

all(cells %in% cells_ori)   # T
all(cells_ori %in% cells)   # T
identical(sort(cells), sort(cells_ori))   # T

sum(geneset_ori %in% EXP_TPM$Gene_Name)            # 2084/2128
sum(geneset_ori %in% EXP_FPKM$Gene_Name)           # 2075/2128
sum(geneset_ori %in% colnames(EXP_PaccMann_Ori))   # 2102/2128


# Calculate the average of columns excluding the specified columns
avg_exp = T
if (avg_exp) {
  col_except = c("Gene_Name", "Gene_ID")
  EXP_TPM = EXP_TPM[, lapply(.SD, mean), by=Gene_Name, .SDcols=setdiff(names(EXP_TPM), col_except)]      # 56450
  EXP_FPKM = EXP_FPKM[, lapply(.SD, mean), by=Gene_Name, .SDcols=setdiff(names(EXP_FPKM), col_except)]   # 53320
  EXP_TPM = EXP_TPM %>% as.data.frame
  EXP_FPKM = EXP_FPKM %>% as.data.frame
  
  rownames(EXP_TPM) = EXP_TPM$Gene_Name
  rownames(EXP_FPKM) = EXP_FPKM$Gene_Name
  EXP_TPM = EXP_TPM[, -1]
  EXP_FPKM = EXP_FPKM[, -1]
  
} else {
  EXP_TPM = EXP_TPM %>% as.data.frame
  EXP_FPKM = EXP_FPKM %>% as.data.frame
  EXP_TPM = EXP_TPM[!duplicated(EXP_TPM$Gene_Name), ]      # 56450
  EXP_FPKM = EXP_FPKM[!duplicated(EXP_RPKM$Gene_Name), ]   # 53320
  
  rownames(EXP_TPM) = EXP_TPM$Gene_Name
  rownames(EXP_FPKM) = EXP_FPKM$Gene_Name
  EXP_TPM = EXP_TPM[, -c(1:2)]
  EXP_FPKM = EXP_FPKM[, -c(1:2)]
}


colnames(EXP_TPM) = colnames(EXP_TPM) %>% strsplit(",") %>% sapply(function(x) x[1])
colnames(EXP_FPKM) = colnames(EXP_FPKM) %>% strsplit(",") %>% sapply(function(x) x[1])
EXP_TPM = EXP_TPM %>% t %>% as.data.frame
EXP_FPKM = EXP_FPKM %>% t %>% as.data.frame

EXP_TPM = EXP_TPM[cells_ori, colnames(EXP_TPM) %in% geneset_ori]      # 457 x 2084
EXP_FPKM = EXP_FPKM[cells_ori, colnames(EXP_FPKM) %in% geneset_ori]   # 457 x 2075
EXP_PaccMann_Sub = EXP_PaccMann_Ori[, colnames(EXP_PaccMann_Ori) %in% geneset_ori]   # 457 x 2102

EXP_TPM[is.na(EXP_TPM)] = 0
EXP_FPKM[is.na(EXP_FPKM)] = 0
EXP_PaccMann_Sub[is.na(EXP_PaccMann_Sub)] = 0

EXP_TPM = EXP_TPM %>% asinh
EXP_FPKM = EXP_FPKM %>% asinh

create_dist_exp = function(EXP, EXP_Ori, mut=F, levels_mut=c(0, 1)) {
  cells_int = intersect(rownames(EXP_Ori), rownames(EXP))
  genes_int = intersect(colnames(EXP_Ori), colnames(EXP))
  sprintf("Cells : %s / %s", length(cells_int), nrow(EXP_Ori)) %>% print
  sprintf("Genes : %s / %s", length(genes_int), ncol(EXP_Ori)) %>% print
  
  Dist_EXP = data.frame(EXP = EXP[cells_int, genes_int] %>% unlist, 
                        EXP_Ori = EXP_Ori[cells_int, genes_int] %>% unlist)
  
  if (mut) {
    Dist_EXP = Dist_EXP %>% 
      mutate(EXP = EXP %>% factor(levels=levels_mut), 
             EXP_Ori = EXP_Ori %>% factor(levels=levels_mut)) %>% 
      with(caret::confusionMatrix(data=EXP, reference=EXP_Ori))
    Dist_EXP = Dist_EXP$table %>% as.data.frame
    
    Dist_EXP = Dist_EXP %>% 
      mutate(EXP = Prediction %>% factor(levels=levels_mut), 
             EXP_Ori = Reference %>% factor(levels=levels_mut)) %>% 
      subset(select=c(EXP, EXP_Ori, Freq))
  }
  return(Dist_EXP)
}

Dist_EXP_TPM = create_dist_exp(EXP_TPM, EXP_PaccMann_Ori)
Dist_EXP_FPKM = create_dist_exp(EXP_FPKM, EXP_PaccMann_Ori)

width = 15
height = 15
omics = c("TPM", "FPKM")
main = sprintf("Distribution of EXP [PaccMann, %s]", omics)
type = c("GDSC", "GDSC, Github")

labs = sprintf("Expression [%s]", type)
Dist_EXP_TPM %>% plot_def(EXP, EXP_Ori, main=main[1], xlab=labs[1], ylab=labs[2],
                          alpha=0.5, axis_tl=25, axis_tx=20, 
                          width=width, height=height, dpi=1200, 
                          xy_line=T, raster=T, save=T, save_svg=T)

labs = sprintf("Expression [%s]", type)
Dist_EXP_FPKM %>% plot_def(EXP, EXP_Ori, main=main[2], xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20, 
                           width=width, height=height, dpi=1200, 
                           xy_line=T, raster=T, save=T, save_svg=T)


# IC50 data must contain cell-line names?
# Or transform cell-line names into sanger ids?

dir = "../../processed_data/ic50_data/GDSC"

file = sprintf("%s/IC50_GDSC.csv", dir)
IC50_GDSC = read.csv(file)

file = sprintf("%s/IC50_GDSC1.csv", dir)
IC50_GDSC1 = read.csv(file)

file = sprintf("%s/IC50_GDSC2.csv", dir)
IC50_GDSC2 = read.csv(file)

# Cell-line names
sum(cells_ori %in% IC50_GDSC$CELL_LINE_NAME)    # 395 [from 457]
sum(cells_ori %in% IC50_GDSC1$CELL_LINE_NAME)   # 391 [from 457]
sum(cells_ori %in% IC50_GDSC2$CELL_LINE_NAME)   # 389 [from 457]


col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "CELL_LINE_NAME", "DRUG_CID", "LN_IC50")
IC50_GDSC = IC50_GDSC[, col] %>% subset(CELL_LINE_NAME %in% cells)     # 373681 > 150204
IC50_GDSC1 = IC50_GDSC1[, col] %>% subset(CELL_LINE_NAME %in% cells)   # 266573 > 107879
IC50_GDSC2 = IC50_GDSC2[, col] %>% subset(CELL_LINE_NAME %in% cells)   # 199687 > 79704

colnames(IC50_GDSC)[c(1, 5)] = c("Cell", "Drug")
colnames(IC50_GDSC1)[c(1, 5)] = c("Cell", "Drug")
colnames(IC50_GDSC2)[c(1, 5)] = c("Cell", "Drug")

dir = mkdir("../../do_better/PaccMann/_data")
file1 = sprintf("%s/IC50_GDSC.csv", dir)
file2 = sprintf("%s/IC50_GDSC1.csv", dir)
file3 = sprintf("%s/IC50_GDSC2.csv", dir)

write.csv(IC50_GDSC, file=file1, row.names=T)
write.csv(IC50_GDSC1, file=file2, row.names=T)
write.csv(IC50_GDSC2, file=file3, row.names=T)

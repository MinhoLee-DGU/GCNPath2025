#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(openxlsx))

source("../functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)




##### 2. Process data for TGSA

file = "Supplementary_Tables.xlsx"
Genes = read.xlsx(file, sheet="Table S4")
geneset_ori = Genes$ENTREZID %>% unique   # 1014

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
TPM_Ent = fread_def(file, col_numeric=T)

geneset = colnames(TPM_Ent)[colnames(TPM_Ent) %in% geneset_ori]   # 1001
TPM_Ent = TPM_Ent[, col]   # 1431 x 1001

dir = mkdir("../../benchmark_test/RF/_data")
file = sprintf("%s/SANGER_RNA_TPM.csv", dir)
write.csv(TPM_Ent, file=file, row.names=T)




##### 3. Process CCLE TPM Data

# Aggregate Omics [Union]
fill_col_na = function(df, col, value=0) {
  df = df %>% as.data.frame
  col_na = col[!(col %in% colnames(df))]
  sprintf("The number of columns missing : %s", len(col_na)) %>% print
  df[, col_na] = value
  df = df[, col]
  return(df)
}

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/OmicsExpressionProteinCodingGenesTPMLogp1.csv", dir)
EXP_CCLE = fread_def(file, na.strings="", header=T)   # 1450 x 19193

col_ent = gsub("(.*) \\((.*)\\)", "\\2", colnames(EXP_CCLE))
EXP_CCLE = EXP_CCLE %>% as.data.frame %>% setNames(col_ent)
EXP_CCLE[is.na(EXP_CCLE)] = 0

# Rename Cells [BROAD > SANGER]
idx = match(rownames(EXP_CCLE), Anno_Cells$BROAD_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_CCLE = EXP_CCLE[!is.na(cells), ]
rownames(EXP_CCLE) = na.omit(cells)   # 1310 x 37566 [-140]

# Filter Genes
EXP_RF_CCLE = EXP_CCLE %>% fill_col_na(geneset)   # 1310 x 1001 [NA 21]

dir = "../../benchmark_test/RF/_data"
file = sprintf("%s/CCLE_RNA_TPM.csv", dir)
write.csv(EXP_RF_CCLE, row.names=T, file=file)




##### 4. Process GDSC Microarray Data

dir = "../../benchmark_test/HiDRA/_data"
file = sprintf("%s/GDSC_RNA_Array.csv", dir)
EXP_GDSC = fread_def(file, header=T, na.strings="")   # 968 x 17419

# Rename Cells [COSMIC > SANGER]
idx = match(rownames(EXP_GDSC), Anno_Cells$COSMIC_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_GDSC = EXP_GDSC[!is.na(idx), ]
rownames(EXP_GDSC) = na.omit(cells)   # 966 x 706 [-2]

# Rename Genes [Symbol > Entrez]
idx = match(colnames(EXP_GDSC), Anno_Genes$HGNC_SYMBOL)
genes = Anno_Genes$ENTREZ_ID[idx]
EXP_GDSC_Ent = EXP_GDSC[, !is.na(idx)]
colnames(EXP_GDSC_Ent) = na.omit(genes)   # 968 x 16669 [-753]

# Filter Genes
EXP_RF_GDSC = EXP_GDSC_Ent %>% fill_col_na(geneset)       # 966 x 1001 [NA 99]

dir = "../../benchmark_test/RF/_data"
file = sprintf("%s/GDSC_RNA_Array.csv", dir)
write.csv(EXP_RF_GDSC, row.names=T, file=file)

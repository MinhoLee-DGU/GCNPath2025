#!/usr/bin/env Rscript

##### 1. Packages and Data

source("../functions.R")
loadings()

dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/gene_identifiers_20191101.csv", dir)
Anno_Genes = read.csv(file, header=T, na.strings="")
colnames(Anno_Genes) = colnames(Anno_Genes) %>% toupper

file = sprintf("%s/model_list_20230517.csv", dir)
Anno_Cells = read.csv(file, header=T, na.strings="")

dir = "../../raw_data/GDSC"
file = sprintf("%s/Cell_listThu Jun  8 07_41_54 2023.csv", dir)
Anno_Cells_GDSC = read.csv(file, header=T, na.strings="")

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/Cell_lines_annotations_20181226.txt", dir)
Anno_Cells_CCLE = read.csv(file, sep="\t")


# Cf. Be careful when using the function match
# Set parameter incomlarables=NA in order not to match NA values
# 
# a = c(NA, NA, 3, 4)
# b = c(1, 2, NA, 4)
# idx = match(a, b)   # 3, 3, NA, 4
# 
# a = c(NA, NA, 3, 4)
# b = c(1, 2, NA, 4)
# idx = match(a, b, incomparables=NA)   # NA, NA, NA, 4

match_multi = function(df1, df2, col1, col2=NULL, values_na=NA) {
  
  idx = NULL
  if (is.null(col2)) col2 = col1
  if (length(col1)!=length(col2)){
    stop("The lengths of col1 & col2 must be same!")
  }
  
  for (i in 1:length(col1)) {
    idx_temp = match(df1[[col1[i]]], df2[[col2[i]]], incomparables=values_na)
    if (is.null(idx)) {
      idx = idx_temp
    } else {
      idx = ifelse(!is.na(idx), idx, idx_temp)
    }
  }
  
  return(idx)
}

# Confirmed that match_multi properly work!
# a = 1:12 %>% matrix(ncol=2) %>% as.data.frame
# b = a
# a$V1[2:3] = NA
# match_multi(a, b, col1="V1", col2="V1")
# # 1 NA NA 4 5 6
# match_multi(b, b, col1="V1", col2="V1")
# # 1 2 3 4 5 6


##### 2. Integrate cell annotations [SANGER Passport, GDSC, CCLE]
# To get TCGA code, SANGER IDs and BROAD IDs of each cell lines

# SANGER Passport & GDSC
col = c("model_id", "BROAD_ID", "COSMIC_ID", "CCLE_ID", "RRID", 
        "model_name", "synonyms", "tissue", "cancer_type", "cancer_type_detail")

Anno_Cells = Anno_Cells[, col]
colnames(Anno_Cells)[1] = "SANGER_MODEL_ID"
colnames(Anno_Cells) = colnames(Anno_Cells) %>% toupper

colnames(Anno_Cells_GDSC)[c(2, 4)] = c("SANGER_MODEL_ID", "TCGA_CODE")
colnames(Anno_Cells_GDSC) = gsub("\\.", "_", colnames(Anno_Cells_GDSC)) %>% toupper
Anno_Cells_GDSC$TCGA_CODE %>% is.na %>% sum   # 8 [1939]


col1 = c("SANGER_MODEL_ID", "COSMIC_ID", "MODEL_NAME")
col2 = c("SANGER_MODEL_ID", "COSMIC_ID", "CELL_LINE_NAME")
idx = match_multi(Anno_Cells, Anno_Cells_GDSC, col1=col1, col2=col2)
tcga_code = Anno_Cells_GDSC$TCGA_CODE[idx]

Anno_Cells = Anno_Cells %>% 
  mutate(TCGA_CODE = tcga_code) %>% 
  relocate(TCGA_CODE, .after=SYNONYMS)

Anno_Cells$TCGA_CODE %>% is.na %>% sum   # 1126


# SANGER Passport & CCLE_DepMap
col = c("SangerModelID", "ModelID", "COSMICID", "CCLEName", "RRID", "CellLineName", 
        "StrippedCellLineName", "OncotreeCode", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype")

Anno_Cells_CCLE = Anno_Cells_CCLE[, col]

col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "CCLE_ID", "RRID", "MODEL_NAME", 
        "SYNONYMS", "TCGA_CODE", "TISSUE", "CANCER_TYPE", "CANCER_TYPE_DETAIL")

colnames(Anno_Cells_CCLE) = col

# [Anno_Cells, SANGER] Fill BROAD_ID using other unique cell ids
col = c("SANGER_MODEL_ID", "COSMIC_ID", "CCLE_ID", "RRID")
idx = match_multi(Anno_Cells, Anno_Cells_CCLE, col1=col)
broad_id = Anno_Cells_CCLE$BROAD_ID[idx]
Anno_Cells$BROAD_ID %>% is.na %>% sum   # 384

Anno_Cells = Anno_Cells %>% 
  mutate(BROAD_ID = ifelse(!is.na(BROAD_ID), BROAD_ID, broad_id))
Anno_Cells$BROAD_ID %>% is.na %>% sum   # 373

# [Anno_Cells, SANGER] Fill TCGA_Code using unique cell ids
col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "CCLE_ID", "RRID")
idx = match_multi(Anno_Cells, Anno_Cells_CCLE, col1=col)
tcga_code = Anno_Cells_CCLE$TCGA_CODE[idx]
Anno_Cells$TCGA_CODE %>% is.na %>% sum   # 1126

Anno_Cells = Anno_Cells %>% 
  mutate(TCGA_CODE = ifelse(!is.na(TCGA_CODE), TCGA_CODE, tcga_code))
Anno_Cells$TCGA_CODE %>% is.na %>% sum   # 487
# Anno_Cells$TCGA_CODE %>% table


dir = mkdir("../../processed_data/cell_data/SANGER_Passports")

file = sprintf("%s/Anno_Cells.csv", dir)
write.csv(Anno_Cells, row.names=F, file=file)

file = sprintf("%s/Anno_Genes.csv", dir)
write.csv(Anno_Genes, row.names=F, file=file)

dir = mkdir("../../processed_data/cell_data/CCLE_DepMap")
file = sprintf("%s/Anno_Cells.csv", dir)
write.csv(Anno_Cells_CCLE, row.names=F, file=file)



##### 2. RNA Processing [TPM, SANGER Passports]

# fpkm_to_tpm = function(RNA, do_log=T) {
#   
#   # Input : FPKM [Genes x Samples]
#   # Output : TPM [Genes x Samples]
#   # TPM = FPKM x 10^6 / sum(FPKM)
#   
#   genes = RNA %>% rownames
#   fpkm_per_sample = RNA %>% colSums
#   RNA = RNA %>% apply(1, function(x) x/fpkm_per_sample) %>% t %>% as.data.frame
#   if (do_log) RNA = log2(RNA*10^6 + 1) %>% as.data.frame
#   rownames(RNA) = genes
#   return(RNA)
# }

dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/rnaseq_tpm_20220624.csv", dir)
TPM = read.csv(file, header=T, row.names=1, na.strings="")

symbols = TPM$X[5:nrow(TPM)]
cells = TPM[1, 2:ncol(TPM)] %>% as.character
TPM = TPM[5:nrow(TPM), 2:ncol(TPM)]

# Fill NAs in TPM Data
TPM %>% is.na %>% sum   # 473793
TPM[is.na(TPM)] = 0

TPM %>% is.na %>% sum       # 0
TPM = TPM %>% mutate_if(is.character, as.numeric)
TPM %>% colSums %>% range   # [969887.4, 1000032.7]

# Log2 Normalization with pseudo count 1
cells %>% unique %>% length     # 1431 [from 1431]
symbols %>% unique %>% length   # 37600 [from 37602]
TPM = log(TPM+1, 2) %>% as.data.frame

# Drop duplicated genes in symbol
idx = symbols %>% duplicated
TPM_Sym = TPM[!idx, ]
rownames(TPM_Sym) = symbols[!idx]
TPM_Sym = TPM_Sym %>% t %>% as.data.frame
# Genes x Cells  >  Cells x Genes

TPM_Ent = TPM
idx = match(rownames(TPM_Ent), Anno_Genes$GENE_ID)
entrez_ids = Anno_Genes$ENTREZ_ID[idx] %>% as.character
entrez_ids %>% is.na %>% sum   # 36

idx = entrez_ids %>% is.na
TPM_Ent = TPM_Ent[!idx, ]
rownames(TPM_Ent) = entrez_ids[!idx]
TPM_Ent = TPM_Ent %>% t %>% as.data.frame
# Genes x Cells  >  Cells x Genes


dir = "../../processed_data/cell_data/SANGER_Passports"
file1 = sprintf("%s/TPM_Sym.csv", dir)
file2 = sprintf("%s/TPM_Ent.csv", dir)
file = sprintf("%s/TPM.RData", dir)

write.csv(TPM_Sym, file=file1, quote=F)
write.csv(TPM_Ent, file=file2, quote=F)
save(Anno_Cells, Anno_Genes, TPM_Sym, TPM_Ent, file=file)



##### 3. RNA Processing [Microarray, GDSC]

dir = "../../raw_data/GDSC"
file = sprintf("%s/Cell_line_RMA_proc_basalExp.txt", dir)
RNA_Array_GDSC = read.csv(file, na.strings="", header=T, sep="\t")

RNA_Array_GDSC[, 3:1020] %>% is.na %>% sum                     # 0
RNA_Array_GDSC$GENE_SYMBOLS %>% is.na %>% sum                  # 318
RNA_Array_GDSC = subset(RNA_Array_GDSC, !is.na(GENE_SYMBOLS))  # 17419

sources = sub("(.*)Source:(.*);Acc:(.*)]", "\\2", RNA_Array_GDSC$GENE_title)
symbols = sub("(.*)Source:(.*);Acc:(.*)]", "\\3", RNA_Array_GDSC$GENE_title)
sources %>% unique   # HGNC

hgnc_ids = sprintf("HGNC:%s", symbols)
intersect(hgnc_ids, Anno_Genes$HGNC_ID) %>% length   # 17408 [17419]

idx1 = match(hgnc_ids, Anno_Genes$HGNC_ID)
idx2 = match(RNA_Array_GDSC$GENE_SYMBOLS, Anno_Genes$HGNC_SYMBOL)
idx3 = match(RNA_Array_GDSC$GENE_SYMBOLS, Anno_Genes$COSMIC_GENE_SYMBOL)

Anno_Genes$ENTREZ_ID[idx1] %>% na   # [NA] 14
Anno_Genes$ENTREZ_ID[idx2] %>% na   # [NA] 753
Anno_Genes$ENTREZ_ID[idx3] %>% na   # [NA] 1472


RNA_Array_GDSC$GENE_title = hgnc_ids
RNA_Array_GDSC$ENTREZ_ID = Anno_Genes$ENTREZ_ID[idx1]
RNA_Array_GDSC = RNA_Array_GDSC %>% relocate(ENTREZ_ID, .after=GENE_title)

cells = colnames(RNA_Array_GDSC)[4:ncol(RNA_Array_GDSC)]
cells = sub("DATA\\.", "", cells)
col_dot = cells[grep("\\.", cells)]
col_dot_ = sub("\\.(.*)", "", col_dot)

# Those cells are duplicated...
col_dot %in% (Anno_Cells$COSMIC_ID %>% unique)    # F, F, F, F
col_dot_ %in% (Anno_Cells$COSMIC_ID %>% unique)   # T, T, T, T
# "1503362.1" "1330983.1" "909976.1"  "905954.1"

colnames(RNA_Array_GDSC)[4:ncol(RNA_Array_GDSC)] = cells
idx = (colnames(RNA_Array_GDSC) %in% col_dot)   # 4
RNA_Array_GDSC = RNA_Array_GDSC[, !idx]         # 1021 > 1017


RNA_Array_GDSC_Sym = RNA_Array_GDSC
RNA_Array_GDSC_Ent = RNA_Array_GDSC %>% subset(!is.na(ENTREZ_ID))
rownames(RNA_Array_GDSC_Sym) = RNA_Array_GDSC_Sym$GENE_SYMBOLS
rownames(RNA_Array_GDSC_Ent) = RNA_Array_GDSC_Ent$ENTREZ_ID

RNA_Array_GDSC_Sym = RNA_Array_GDSC_Sym[, -c(1:3)]
RNA_Array_GDSC_Ent = RNA_Array_GDSC_Ent[, -c(1:3)]
RNA_Array_GDSC_Sym = RNA_Array_GDSC_Sym %>% t %>% as.data.frame   # 1014 x 17419
RNA_Array_GDSC_Ent = RNA_Array_GDSC_Ent %>% t %>% as.data.frame   # 1014 x 17405

idx = match(rownames(RNA_Array_GDSC_Sym), Anno_Cells$COSMIC_ID)
cells_sanger = Anno_Cells$SANGER_MODEL_ID[idx]

idx = cells_sanger %>% is.na
RNA_Array_GDSC_Sym = RNA_Array_GDSC_Sym[!idx, ]   # 1014 > 1006
RNA_Array_GDSC_Ent = RNA_Array_GDSC_Ent[!idx, ]   # 1014 > 1006
rownames(RNA_Array_GDSC_Sym) = cells_sanger[!idx]
rownames(RNA_Array_GDSC_Ent) = cells_sanger[!idx]


dir = mkdir("../../processed_data/cell_data/GDSC")
file1 = sprintf("%s/RNA_Array_Sym.csv", dir)
file2 = sprintf("%s/RNA_Array_Ent.csv", dir)
file3 = sprintf("%s/Anno_Cells.csv", dir)

write.csv(RNA_Array_GDSC_Sym, file=file1, quote=F)
write.csv(RNA_Array_GDSC_Ent, file=file2, quote=F)
write.csv(Anno_Cells_GDSC, file=file3, row.names=F)



##### 4. RNA Processing [RNA-Seq, CCLE Depmap 23Q2]

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/OmicsExpressionProteinCodingGenesTPMLogp1.csv", dir)
TPM_CCLE = read.csv(file, row.names=1, check.names=F)   # 1450 x 19193
colnames(TPM_CCLE) = gsub("(.*) \\((.*)\\)", "\\1", colnames(TPM_CCLE))

TPM_CCLE_BROAD_Sym = TPM_CCLE
TPM_CCLE_BROAD_Ent = TPM_CCLE
idx = match(colnames(TPM_CCLE_BROAD_Ent), Anno_Genes$HGNC_SYMBOL)

TPM_CCLE_BROAD_Ent = TPM_CCLE_BROAD_Ent[, !is.na(idx)]   # 19193 > 18782
colnames(TPM_CCLE_BROAD_Ent) = Anno_Genes$ENTREZ_ID[na.omit(idx)]

TPM_CCLE_Sym = TPM_CCLE_BROAD_Sym
TPM_CCLE_Ent = TPM_CCLE_BROAD_Ent

idx = match(rownames(TPM_CCLE_Sym), Anno_Cells$BROAD_ID)
TPM_CCLE_Sym = TPM_CCLE_Sym[!is.na(idx), ]   # 1450 > 1310
TPM_CCLE_Ent = TPM_CCLE_Ent[!is.na(idx), ]   # 1450 > 1310
rownames(TPM_CCLE_Sym) = Anno_Cells$SANGER_MODEL_ID[na.omit(idx)]
rownames(TPM_CCLE_Ent) = Anno_Cells$SANGER_MODEL_ID[na.omit(idx)]


dir = mkdir("../../processed_data/cell_data/CCLE_DepMap")
file1 = sprintf("%s/TPM_BROAD_ID_Sym.csv", dir)
file2 = sprintf("%s/TPM_BROAD_ID_Ent.csv", dir)
file3 = sprintf("%s/TPM_Sym.csv", dir)
file4 = sprintf("%s/TPM_Ent.csv", dir)

write.csv(TPM_CCLE_BROAD_Sym, file=file1, quote=F)
write.csv(TPM_CCLE_BROAD_Ent, file=file2, quote=F)
write.csv(TPM_CCLE_Sym, file=file3, quote=F)
write.csv(TPM_CCLE_Ent, file=file4, quote=F)

rm(TPM, TPM_CCLE, RNA_Array_GDSC)


supplementary = F
if (supplementary) {
  suppressMessages(library(openxlsx))
  
  ### [Supplementary Data] Supplementary Data 1
  dir = "../../processed_data/ic50_data/GDSC"
  file = sprintf("%s/IC50_GDSC.csv", dir)
  IC50_GDSC = read.csv(file)
  
  Anno_Cells_ = Anno_Cells %>% subset(SANGER_MODEL_ID %in% rownames(TPM_Ent))
  Anno_Cells_$GDSC_Cell = Anno_Cells_$SANGER_MODEL_ID %in% unique(IC50_GDSC$SANGER_MODEL_ID)
  Anno_Cells_ = Anno_Cells_ %>% arrange(SANGER_MODEL_ID) %>% as.data.frame
  Anno_Cells_$GDSC_Cell %>% sum   # 972
  
  dir = "../../processed_data/cell_data/SANGER_Passports"
  sheets = "Supplementary Data 1"
  file = sprintf("%s/%s.xlsx", dir, sheets)
  write.xlsx(Anno_Cells_, file=file, sheetName=sheets, rowNames=F)
}

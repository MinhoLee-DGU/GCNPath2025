#!/usr/bin/env Rscript

##### 1. Packages and Data

source("../functions.R")
loadings()

# Gene Annotation from SANGER Cell Passport
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

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


# This code is for preparing test dataset only...
# Train dataset of GDSC microarray should be created first
# Run first "HiDRA_FeatureGeneration.py"
dir = "../../do_better/HiDRA/Prediction/input_dir"
file_list = sprintf("%s/%s.csv", dir, 0:185)

EXP_List_Ori = list()
Gene_List_Ori = list()
for (i in 1:length(file_list)) {
  EXP_List_Ori[[i]] = fread_def(file_list[i])
  Gene_List_Ori[[i]] = EXP_List_Ori[[i]] %>% colnames
}

dir = "../../do_better/HiDRA/_data/ProcessedFile"
file_list = sprintf("%s/%s.csv", dir, 0:185)

EXP_List = list()
Gene_List = list()
for (i in 1:length(file_list)) {
  EXP_List[[i]] = fread_def(file_list[i])
  Gene_List[[i]] = EXP_List[[i]] %>% colnames
}

EXP = Reduce(cbind, EXP_List)
EXP_Ori = Reduce(cbind, EXP_List_Ori)

EXP = EXP[, !duplicated(colnames(EXP))]
EXP_Ori = EXP_Ori[, !duplicated(colnames(EXP_Ori))]

geneset = Gene_List %>% unlist %>% unique           # 4590
geneset_ori = Gene_List_Ori %>% unlist %>% unique   # 4592

rownames(EXP_Ori)[!(rownames(EXP_Ori) %in% Anno_Cells$MODEL_NAME)]
# [Name] H2731, H2369, H2804
# [COSMIC] 1240134, 1290808, 1240136

idx = match(rownames(EXP_Ori), Anno_Cells$MODEL_NAME)
rownames(EXP_Ori)[!is.na(idx)] = Anno_Cells$COSMIC_ID[na.omit(idx)]
rownames(EXP_Ori)[rownames(EXP_Ori)=="H2731"] = 1240134
rownames(EXP_Ori)[rownames(EXP_Ori)=="H2369"] = 1290808
rownames(EXP_Ori)[rownames(EXP_Ori)=="H2804"] = 1240136

Dist_EXP = create_dist_exp(EXP, EXP_Ori)
# Cells : 10 / 10
# Genes : 4590 / 4592


width = 15
height = 15
main = "Distribution of EXP [HiDRA]"
type = c("GDSC", "GDSC, Github")

labs = sprintf("Expression [%s]", type)
Dist_EXP %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                      alpha=0.5, axis_tl=25, axis_tx=20, 
                      width=width, height=height, dpi=1200, 
                      xy_line=T, raster=T, save=T, save_svg=T)


# Does GDSC microarray data lack two genes?
# dir = "../../do_better/HiDRA/_data"
# file = sprintf("%s/Cell_line_RMA_proc_basalExp.txt", dir)
# a = fread(file, header=T)   # 968 x 17419
# 
# sum(geneset_ori %in% EXP$GENE_SYMBOLS)   # 4590
# geneset_ori[!(geneset_ori %in% EXP$GENE_SYMBOLS)]   # UQCRHL, OR56A5



##### 2. Process SANGER RNA-Seq Data
# In HiDRA, z-socre in each sample is applied
# This ranking-based method might be sensitive to the total number of genes
# Conduct z-score normalization after remaining these genes only...

fill_col_na = function(df, col, value=0) {
  df = df %>% as.data.frame
  col_na = col[!(col %in% colnames(df))]
  info = sum(col %in% colnames(df))
  info_pct = round(100*info/length(col), 2)
  sprintf("# Columns : %s / %s [%s%%]", info, length(col), info_pct) %>% print
  if (length(col_na)!=0) df[, col_na] = value
  df = df[, col]
  return(df)
}

# Retrieve the total number of genes
# First run HiDRA_FeatureGeneration.py
dir = "../../do_better/HiDRA/_data"
file = sprintf("%s/GDSC_RNA_Array.csv", dir)
geneset_total = read.csv(file, row.names=1, nrow=1, check.names=F) %>% colnames   # 17419

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Sym.csv", dir)
SANGER_TPM = fread_def(file, header=T, check_names=F)
SANGER_TPM = SANGER_TPM %>% fill_col_na(geneset_total)
# Columns : 16316 / 17419 [93.67%]

SANGER_TPM = SANGER_TPM %>% t %>% scale %>% t %>% as.data.frame
SANGER_TPM %>% rowMeans %>% range        # 0
SANGER_TPM %>% apply(1, sd) %>% range    # 1
sum(geneset %in% colnames(SANGER_TPM))   # 4590 [4590]

# check_cols = c()
# Check_Cols = list()
dir = mkdir("../../do_better/HiDRA/_data/ProcessedFile_SANGER")

for (i in 1:length(file_list)) {
  file = sprintf("%s/%s.csv", dir, i-1)
  SANGER_TPM_Temp = SANGER_TPM[, Gene_List[[i]]]
  fwrite(SANGER_TPM_Temp, file=file, row.names=T, col.names=T)
  
  # Do specific gene names are conserved? [Ex. HLA-E > HLA.E]
  # Check_Cols[[i]] = fread_def(file)
  # check_cols_ = identical(colnames(SANGER_TPM_Temp), colnames(Check_Cols[[i]]))
  # check_cols = check_cols %>% c(check_cols_)   # all T
}



# GDSC vs SANGER
# Rename Cells [SANGER > COSMIC]
SANGER_TPM_ = SANGER_TPM
idx = match(rownames(SANGER_TPM_), Anno_Cells$SANGER_MODEL_ID)
cells = Anno_Cells$COSMIC_ID[idx]
SANGER_TPM_ = SANGER_TPM_[!is.na(cells), geneset]
rownames(SANGER_TPM_) = na.omit(cells)

Dist_EXP_SANGER = create_dist_exp(EXP, SANGER_TPM_)
# Cells : 959 / 1045
# Genes : 4590 / 4590

width = 15
height = 15
main = "Distribution of EXP [HiDRA, GDSC & SANGER]"
type = c("GDSC", "SANGER")

labs = sprintf("Expression [%s]", type)
Dist_EXP_SANGER %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                             alpha=0.1, axis_tl=25, axis_tx=20, 
                             width=width, height=height, dpi=600, 
                             xy_line=T, raster=T, save=T, save_svg=T)



##### 3. Process CCLE RNA-Seq Data

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/OmicsExpressionProteinCodingGenesTPMLogp1.csv", dir)
CCLE_TPM = fread_def(file)   # 1450 x 19193
colnames(CCLE_TPM) = gsub("(.*) \\((.*)\\)", "\\1", colnames(CCLE_TPM))
CCLE_TPM = CCLE_TPM %>% fill_col_na(geneset_total)
# Columns : 15963 / 17419 [91.64%]

CCLE_TPM = CCLE_TPM %>% t %>% scale %>% t %>% as.data.frame
CCLE_TPM %>% rowMeans %>% range %>% round(3)       # 0
CCLE_TPM %>% apply(1, sd) %>% range %>% round(3)   # 1
sum(geneset %in% colnames(CCLE_TPM))               # 4590 [4590]

dir = mkdir("../../do_better/HiDRA/_data/ProcessedFile_CCLE")
for (i in 1:length(file_list)) {
  file = sprintf("%s/%s.csv", dir, i-1)
  CCLE_TPM_Temp = CCLE_TPM[, Gene_List[[i]]]
  fwrite(CCLE_TPM_Temp, file=file, row.names=T, col.names=T)
}

# GDSC vs CCLE
CCLE_TPM_ = CCLE_TPM
idx = match(rownames(CCLE_TPM_), Anno_Cells$BROAD_ID)
cells = Anno_Cells$COSMIC_ID[idx]
CCLE_TPM_ = CCLE_TPM_[!is.na(cells), geneset]
rownames(CCLE_TPM_) = na.omit(cells)

Dist_EXP_CCLE = create_dist_exp(EXP, CCLE_TPM_)
# Cells : 700 / 774
# Genes : 4590 / 4590

width = 15
height = 15
main = "Distribution of EXP [HiDRA, GDSC & CCLE]"
type = c("GDSC", "CCLE")

labs = sprintf("Expression [%s]", type)
Dist_EXP_CCLE %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.1, axis_tl=25, axis_tx=20, 
                           width=width, height=height, dpi=600, 
                           xy_line=T, raster=T, save=T, save_svg=T)


# ##### Cf. Cell-line intersection
# 
# file1 = "../../do_better/HiDRA/_data/Cell_line_RMA_proc_basalExp.txt"
# file2 = "../../do_better/HiDRA/_data/GDSC_RNA_Array.csv"
# file3 = "../../do_better/HiDRA/_data/IC50_GDSC.txt"
# 
# EXP_Before = fread(file1, sep="\t", check.names=F)
# EXP_After = fread_def(file2, check_names=F)
# IC50_GDSC = read.csv(file3, sep="\t")
# 
# cells_before = gsub("DATA.", "", colnames(EXP_Before)[3:ncol(EXP_Before)])
# cells_after = rownames(EXP_After)
# cells_ic50 = IC50_GDSC$Cell_COSMIC %>% unique
# 
# intersect(cells_before, cells_ic50) %>% length   # 947
# intersect(cells_after, cells_ic50) %>% length    # 947

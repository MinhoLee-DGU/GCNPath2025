#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
source("../functions.R")
loadings()

# Gene list from 2128_genes.pkl
file = "2128_genes.txt"
geneset_ori = scan(file, what=character())   # 2128

dir = "../_prepare_paccmann"
file = sprintf("%s/gdsc-rnaseq_gene-expression.csv", dir)
EXP_PaccMann_Ori = fread_def(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)


##### 2. Process SANGER TPM data
dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/rnaseq_tpm_20220624.csv", dir)
EXP = fread(file, header=T, na.strings="")

cells_name = EXP[1, 3:ncol(EXP)] %>% as.character
EXP = EXP[5:nrow(EXP), 2:ncol(EXP)]

# Fill NAs in EXP Data
EXP %>% is.na %>% sum   # 473793
EXP[is.na(EXP)] = 0
EXP %>% is.na %>% sum       # 0

col = names(EXP)[-1]
EXP = EXP[, (col):=lapply(.SD, as.numeric, as.is=T), .SDcols=col]
EXP = EXP[, lapply(.SD, mean), by=V2, .SDcols=col]   # 37602 > 37600

symbols = EXP$V2
EXP = EXP[, -"V2", with=F] %>% as.data.frame
rownames(EXP) = symbols

EXP = EXP %>% t %>% asinh %>% as.data.frame
EXP_PaccMann = EXP[, colnames(EXP) %in% geneset_ori]    # 1431 x 2093

EXP_PaccMann_ = EXP_PaccMann
EXP_PaccMann_Ori_ = EXP_PaccMann_Ori
rownames(EXP_PaccMann_) = toupper(cells_name)
rownames(EXP_PaccMann_Ori_) = rownames(EXP_PaccMann_Ori_) %>% toupper
sum(rownames(EXP_PaccMann_Ori_) %in% rownames(EXP_PaccMann_))   # 419/457


# Distributions of EXP between SAGNER & GDSC are similar...
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

Dist_EXP = create_dist_exp(EXP_PaccMann_, EXP_PaccMann_Ori_)

width = 15
height = 15
main = "Distribution of EXP [PaccMann_SANGER]"
type = c("SANGER", "GDSC, Github")

labs = sprintf("Expression [%s]", type)
Dist_EXP %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                      alpha=0.5, axis_tl=25, axis_tx=20, 
                      width=width, height=height, dpi=1200, 
                      xy_line=T, raster=T, save=T, save_svg=T)


dir = mkdir("../../do_better/PaccMann_SANGER/_data")

file = sprintf("%s/2093_genes.txt", dir)
geneset = colnames(EXP_PaccMann)
cat(geneset, file=file, sep="\n")

file = sprintf("%s/EXP.csv", dir)
write.csv(EXP_PaccMann, file=file, row.names=T)



##### Cf. Compare data distribution [Github_CCLE, Github_GDSC]

dir = "../../do_better/PaccMann/data/gene_expression"

file = sprintf("%s/ccle-rnaseq_gene-expression.csv", dir)
EXP_CCLE_Ori = fread_def(file, header=T)
EXP_CCLE_Ori_ = EXP_CCLE_Ori

file = sprintf("%s/gdsc-rma_gene-expression.csv", dir)
EXP_GDSC_Ori = fread_def(file, header=T)
EXP_GDSC_Ori_ = EXP_GDSC_Ori

# SANGER vs CCLE
idx = match(rownames(EXP_CCLE_Ori_), Anno_Cells$CCLE_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_CCLE_Ori_ = EXP_CCLE_Ori_[!is.na(cells), ]
rownames(EXP_CCLE_Ori_) = na.omit(cells)

labs = c("Expression [SANGER]", "Expression [CCLE]")
main = "Distribution of Expression [SANGER & CCLE]"
Dist_EXP_CCLE = create_dist_exp(EXP_PaccMann, EXP_CCLE_Ori_)
Dist_EXP_CCLE %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=20, axis_tx=18, margin=0.4,
                           width=width, height=height, xy_line=T, save=T, save_svg=F)

# SANGER vs GDSC_Github
idx = match(rownames(EXP_GDSC_Ori_), Anno_Cells$COSMIC_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_GDSC_Ori_ = EXP_GDSC_Ori_[!is.na(cells), ]
rownames(EXP_GDSC_Ori_) = na.omit(cells)

labs = c("RNA [Github, GDSC]", "RNA [SANGER]")
main = "Distribution of RNA [Github_GDSC & SANGER]"
Dist_EXP_GDSC = create_dist_exp(EXP_PaccMann, EXP_GDSC_Ori_)
Dist_EXP_GDSC %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=20, axis_tx=18, margin=0.4,
                           width=width, height=height, xy_line=T, save=T, save_svg=F)



##### 3. Process CCLE TPM data

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/OmicsExpressionProteinCodingGenesTPMLogp1.csv", dir)
EXP_CCLE = fread_def(file)   # 1450 x 19193
colnames(EXP_CCLE) = gsub("(.*) \\((.*)\\)", "\\1", colnames(EXP_CCLE))

# Fill NAs in EXP Data
EXP_CCLE %>% is.na %>% sum   # 0
# EXP[is.na(EXP)] = 0
# EXP %>% is.na %>% sum       # 0

# Duplicated genes not exist...
EXP_CCLE %>% colnames %>% length              # 19193
EXP_CCLE %>% colnames %>% unique %>% length   # 19193

# Arcsine Transformation
EXP_CCLE = 2**EXP_CCLE - 1

# EXP from CCLE was truly log2(TPM+1) normalized...
EXP_CCLE %>% rowSums %>% hist
EXP_CCLE %>% rowSums %>% median %>% round(3)   # 948157.7
EXP_CCLE %>% rowSums %>% range %>% round(3)    # 662892.2 971245.6

EXP_CCLE = EXP_CCLE %>% asinh %>% as.data.frame
EXP_CCLE = EXP_CCLE[, colnames(EXP_CCLE) %in% geneset]   # 1450 x 2083 [2093]

dir = mkdir("../../do_better/PaccMann_SANGER/_data")
file = sprintf("%s/EXP_CCLE.csv", dir)
write.csv(EXP_CCLE, file=file, row.names=T)



##### 4. Process GDSC Microarray data

dir = "../../do_better/HiDRA/_data"
file = sprintf("%s/GDSC_RNA_Array.csv", dir)
EXP_GDSC = fread_def(file, header=T, na.strings="")   # 968 x 17419

# Fill NAs in EXP Data
EXP_GDSC %>% is.na %>% sum   # 0
# EXP[is.na(EXP)] = 0
# EXP %>% is.na %>% sum       # 0

# Duplicated genes not exist...
EXP_GDSC %>% colnames %>% length              # 17419
EXP_GDSC %>% colnames %>% unique %>% length   # 17419

EXP_GDSC = EXP_GDSC %>% asinh %>% as.data.frame
EXP_GDSC = EXP_GDSC[, colnames(EXP_GDSC) %in% geneset]   # 968 x 2055 [2093]

dir = mkdir("../../do_better/PaccMann_SANGER/_data")
file = sprintf("%s/EXP_GDSC.csv", dir)
write.csv(EXP_GDSC, file=file, row.names=T)



##### Cf. Compare data distribution

# SANGER vs CCLE
EXP_CCLE_ = EXP_CCLE
idx = match(rownames(EXP_CCLE_), Anno_Cells$BROAD_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_CCLE_ = EXP_CCLE_[!is.na(cells), ]
rownames(EXP_CCLE_) = na.omit(cells)

labs = c("Expression [SANGER]", "Expression [CCLE]")
main = "Distribution of EXP [PaccMann_SANGER, SANGER & CCLE]"
Dist_EXP_CCLE = create_dist_exp(EXP_PaccMann, EXP_CCLE_)
# Cells : 1007 / 1310
# Genes : 2083 / 2083

Dist_EXP_CCLE %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20,
                           width=width, height=height, xy_line=T, raster=T, save=T, save_svg=T)

# SANGER vs GDSC
EXP_GDSC_ = EXP_GDSC
idx = match(rownames(EXP_GDSC_), Anno_Cells$COSMIC_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_GDSC_ = EXP_GDSC_[!is.na(cells), ]
rownames(EXP_GDSC_) = na.omit(cells)

labs = c("Expression [GDSC]", "Expression [SANGER]")
main = "Distribution of EXP [PaccMann_SANGER, SANGER & GDSC]"
Dist_EXP_GDSC = create_dist_exp(EXP_PaccMann, EXP_GDSC_)
# Cells : 959 / 966
# Genes : 2055 / 2055

Dist_EXP_GDSC %>% plot_def(EXP_Ori, EXP, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20,
                           width=width, height=height, xy_line=T, raster=T, save=T, save_svg=T)

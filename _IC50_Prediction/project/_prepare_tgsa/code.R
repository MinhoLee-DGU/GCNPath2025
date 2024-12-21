#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
source("../functions.R")
loadings()

# Cell Annotation from SANGER Cell Passport
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

# Gene Annotation from SANGER Cell Passport
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)



##### 2. Process data [TGSA]

# CGC Genes from the TGSA original data [706]
dir = "../../do_better/TGSA/data/CellLines_DepMap/CCLE_580_18281/census_706"
file_exp = sprintf("%s/exp.csv", dir)
file_mut = sprintf("%s/mu.csv", dir)
file_cnv = sprintf("%s/cn.csv", dir)

EXP_TGSA_Ori = read.csv(file_exp, row.names=1, check.names=F)
MUT_TGSA_Ori = read.csv(file_mut, row.names=1, check.names=F)
CNV_TGSA_Ori = read.csv(file_cnv, row.names=1, check.names=F)

identical(colnames(EXP_TGSA_Ori), colnames(MUT_TGSA_Ori))   # T
identical(colnames(EXP_TGSA_Ori), colnames(CNV_TGSA_Ori))   # T
geneset_ori = gsub("\\(([0-9]*)\\)", "\\1", colnames(EXP_TGSA_Ori)) %>% as.character

EXP_TGSA_Ori = EXP_TGSA_Ori %>% setNames(geneset_ori)
MUT_TGSA_Ori = MUT_TGSA_Ori %>% setNames(geneset_ori)
CNV_TGSA_Ori = CNV_TGSA_Ori %>% setNames(geneset_ori)

Anno_Genes$ENTREZ_ID = Anno_Genes$ENTREZ_ID %>% as.character
sum(geneset_ori %in% Anno_Genes$ENTREZ_ID)                            # 706
Anno_Genes_TGSA = Anno_Genes %>% subset(ENTREZ_ID %in% geneset_ori)   # 706


dir = "../../raw_data/CCLE_DepMap"
file_exp = sprintf("%s/OmicsExpressionProteinCodingGenesTPMLogp1.csv", dir)
file_mut = sprintf("%s/OmicsSomaticMutations.csv", dir)
file_cnv = sprintf("%s/OmicsCNGene.csv", dir)

MUT = fread(file_mut, na.strings="", header=T)       # 1408099 x 56
EXP = fread_def(file_exp, na.strings="", header=T)   # 1450 x 19193
CNV = fread_def(file_cnv, na.strings="", header=T)   # 1814 x 25369


# EXP Data
EXP = EXP %>% as.data.frame
col_sym = gsub("(.*) \\((.*)\\)", "\\1", colnames(EXP))
col_ent = gsub("(.*) \\((.*)\\)", "\\2", colnames(EXP))

EXP_Sym = EXP %>% setNames(col_sym)
EXP_Ent = EXP %>% setNames(col_ent)

colnames(EXP_Ent) %>% class                # character
all(geneset_ori %in% colnames(EXP_Ent))   # T
idx = (colnames(EXP_Ent) %in% geneset_ori)
EXP_TGSA = EXP_Ent[, idx]   # 19193 > 706


# MUT Data
MUT$EntrezGeneID = gsub("\\.0", "", MUT$EntrezGeneID)
MUT_TGSA = MUT[EntrezGeneID %in% geneset_ori]   # 1408099 > 81634
MUT_TGSA$EntrezGeneID %>% unique %>% length      # 698

rm(MUT)
any_ = function(x) ifelse(length(x)>=1, 1, 0)
MUT_TGSA = dcast(MUT_TGSA, ModelID~EntrezGeneID, value.var="EntrezGeneID", fun.agg=any_)
MUT_TGSA = MUT_TGSA %>% as.data.frame

rownames(MUT_TGSA) = MUT_TGSA$ModelID
MUT_TGSA = MUT_TGSA[, -1]
MUT_TGSA %>% unlist %>% table   # 0, 1 only


# CNV Data
CNV = CNV %>% as.data.frame
col_sym = gsub("(.*) \\((.*)\\)", "\\1", colnames(CNV))
col_ent = gsub("(.*) \\((.*)\\)", "\\2", colnames(CNV))

CNV_Sym = CNV %>% setNames(col_sym)
CNV_Ent = CNV %>% setNames(col_ent)

colnames(CNV_Ent) %>% class                # character
all(geneset_ori %in% colnames(CNV_Ent))   # T
idx = (colnames(CNV_Ent) %in% geneset_ori)
CNV_TGSA = CNV_Ent[, idx]   # 25368 > 706
CNV_TGSA %>% unlist %>% hist


# Aggregate Omics [Union]
fill_col_na = function(df, col, value=0) {
  df = df %>% as.data.frame
  col_na = col[!(col %in% colnames(df))]
  df[, col_na] = value
  df = df[, col]
  return(df)
}

MUT_TGSA = MUT_TGSA %>% fill_col_na(geneset_ori, value=0)   # 698 > 706
# EXP_TGSA = EXP_TGSA %>% fill_col_na(geneset_ori, value=0)   # 706 > 706
# CNV_TGSA = CNV_TGSA %>% fill_col_na(geneset_ori, value=1)   # 706 > 706

identical(colnames(MUT_TGSA), colnames(EXP_TGSA))   # T
identical(colnames(MUT_TGSA), colnames(CNV_TGSA))   # T

cells = rownames(MUT_TGSA) %>% intersect(rownames(EXP_TGSA)) %>% intersect(rownames(CNV_TGSA))   # 1399
MUT_TGSA = MUT_TGSA[cells, ]   # 1738 > 1399
EXP_TGSA = EXP_TGSA[cells, ]   # 1450 > 1399
CNV_TGSA = CNV_TGSA[cells, ]   # 1814 > 1399


# Confirmed that values are very similar as original
# Some values can vary due to data processing and/or database update

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

Dist_MUT = create_dist_exp(MUT_TGSA, MUT_TGSA_Ori, mut=T)
Dist_EXP = create_dist_exp(EXP_TGSA, EXP_TGSA_Ori)
Dist_CNV = create_dist_exp(CNV_TGSA, CNV_TGSA_Ori)

width = 15
height = 15
omics = c("MUT", "EXP", "CNV")
main = sprintf("Distribution of %s [TGSA]", omics)
type = c("CCLE", "CCLE, Github")

labs = sprintf("Mutation [%s]", type)
color = scale_fill_gradient(low="grey90", high="grey90")
Dist_MUT %>% grid_def(EXP, EXP_Ori, fill=Freq, main=main[1], color=color,
                      xlab=labs[1], ylab=labs[2], axis_tl=25, axis_tx=20, 
                      margin=0.4, angle=0, width=width, height=height, 
                      mean_summary=F, legend=F, save=T, save_svg=T)

labs = sprintf("Expression [%s]", type)
Dist_EXP %>% plot_def(EXP, EXP_Ori, main=main[2], xlab=labs[1], ylab=labs[2],
                      alpha=0.5, axis_tl=25, axis_tx=20, 
                      width=width, height=height, dpi=1200, 
                      xy_line=T, raster=T, save=T, save_svg=T)

labs = sprintf("CNV [%s]", type)
Dist_CNV %>% plot_def(EXP, EXP_Ori, main=main[3], xlab=labs[1], ylab=labs[2],
                      alpha=0.5, axis_tl=25, axis_tx=20, 
                      width=width, height=height, dpi=1200, 
                      xy_line=T, raster=T, save=T, save_svg=T)


dir = mkdir("../../do_better/TGSA/_data")
file_exp = sprintf("%s/EXP.csv", dir)
file_mut = sprintf("%s/MUT.csv", dir)
file_cnv = sprintf("%s/CNV.csv", dir)

write.csv(EXP_TGSA, row.names=T, file=file_exp)
write.csv(MUT_TGSA, row.names=T, file=file_mut)
write.csv(CNV_TGSA, row.names=T, file=file_cnv)

EXP_Total_Scaled = EXP_Ent[rownames(EXP_TGSA), ] %>% scale %>% as.data.frame
EXP_Total_Scaled %>% sapply(mean) %>% range %>% round(3)   # [0, 0]
EXP_Total_Scaled %>% sapply(sd) %>% range %>% round(3)     # [1, 1]

file = sprintf("%s/EXP_Total_Scaled.csv", dir)
write.csv(EXP_Total_Scaled, row.names=T, file=file)



##### 3. Process GDSC Microarray Data

dir = "../../do_better/HiDRA/_data"
file = sprintf("%s/GDSC_RNA_Array.csv", dir)
EXP_GDSC = fread_def(file, header=T, na.strings="")   # 968 x 17419

# Rename Cells [COSMIC > BROAD]
idx = match(rownames(EXP_GDSC), Anno_Cells$COSMIC_ID)
cells = Anno_Cells$BROAD_ID[idx]
EXP_GDSC = EXP_GDSC[!is.na(cells), ]
rownames(EXP_GDSC) = na.omit(cells)   # 966 x 17419 [-2]

# Rename Genes [Symbol > Entrez]
idx = match(colnames(EXP_GDSC), Anno_Genes$HGNC_SYMBOL)
genes = Anno_Genes$ENTREZ_ID[idx]
EXP_GDSC = EXP_GDSC[, !is.na(genes)]
colnames(EXP_GDSC) = na.omit(genes)   # 966 x 16669 [-750]

# Filter Genes [CGC 706]
EXP_TGSA_GDSC = EXP_GDSC %>% fill_col_na(geneset_ori)   # 966 x 706

# Retrieve Cells [Intersecting in MUT, EXP, CNV]
cells = Reduce(intersect, list(rownames(MUT_TGSA), rownames(EXP_TGSA_GDSC), rownames(CNV_TGSA)))   # 697
EXP_TGSA_GDSC = EXP_TGSA_GDSC[cells, ]
MUT_TGSA_GDSC = MUT_TGSA[cells, ]
CNV_TGSA_GDSC = CNV_TGSA[cells, ]

identical(colnames(EXP_TGSA_GDSC), colnames(MUT_TGSA_GDSC))   # T
identical(colnames(EXP_TGSA_GDSC), colnames(MUT_TGSA_GDSC))   # T

dir = mkdir("../../do_better/TGSA/_data_gdsc")
file1 = sprintf("%s/MUT.csv", dir)
file2 = sprintf("%s/EXP.csv", dir)
file3 = sprintf("%s/CNV.csv", dir)

write.csv(MUT_TGSA_GDSC, file=file1, row.names=T)
write.csv(EXP_TGSA_GDSC, file=file2, row.names=T)
write.csv(CNV_TGSA_GDSC, file=file3, row.names=T)


# Expression for the cell similarity graph
EXP_Total_Scaled_GDSC = EXP_GDSC[rownames(EXP_TGSA_GDSC), ] %>% scale %>% as.data.frame   # 697 x 16669
EXP_Total_Scaled_GDSC %>% sapply(mean) %>% range %>% round(3)   # [0, 0]
EXP_Total_Scaled_GDSC %>% sapply(sd) %>% range %>% round(3)     # [1, 1]

dir = "../../do_better/TGSA/_data_gdsc"
file = sprintf("%s/EXP_Total_Scaled.csv", dir)
write.csv(EXP_Total_Scaled_GDSC, row.names=T, file=file)

# GDSC vs CCLE
labs = c("Expression [CCLE]", "Expression [GDSC]")
main = "Distribution of EXP [TGSA, CCLE & GDSC]"
Dist_EXP_GDSC = create_dist_exp(EXP_TGSA, EXP_TGSA_GDSC)
Dist_EXP_GDSC %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20, 
                           width=width, height=height, dpi=1200, 
                           xy_line=T, raster=T, save=T, save_svg=T)



##### 4. Process SANGER TPM Data
# For TGSA, intersect cells from MUT, EXP, CNV
# For DRPreter, we can utilize data from DRPreter_SANGER

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
EXP_SANGER = fread_def(file, header=T, na.strings="")   # 1431 x 37566

# Rename Cells [COSMIC > BROAD]
idx = match(rownames(EXP_SANGER), Anno_Cells$SANGER_MODEL_ID)
cells = Anno_Cells$BROAD_ID[idx]
EXP_SANGER = EXP_SANGER[!is.na(cells), ]
rownames(EXP_SANGER) = na.omit(cells)   # 1293 x 37566 [-138]

# Filter Genes [CGC 706]
EXP_TGSA_SANGER = EXP_SANGER %>% fill_col_na(geneset_ori)   # 1293 x 706

# Retrieve Cells [Intersecting in MUT, EXP, CNV]
cells = Reduce(intersect, list(rownames(MUT_TGSA), rownames(EXP_TGSA_SANGER), rownames(CNV_TGSA)))   # 984
EXP_TGSA_SANGER = EXP_TGSA_SANGER[cells, colnames(MUT_TGSA_SANGER)]
MUT_TGSA_SANGER = MUT_TGSA[cells, ]
CNV_TGSA_SANGER = CNV_TGSA[cells, ]

identical(colnames(EXP_TGSA_SANGER), colnames(MUT_TGSA_SANGER))   # T
identical(colnames(EXP_TGSA_SANGER), colnames(CNV_TGSA_SANGER))   # T

dir = mkdir("../../do_better/TGSA/_data_sanger")
file1 = sprintf("%s/MUT.csv", dir)
file2 = sprintf("%s/EXP.csv", dir)
file3 = sprintf("%s/CNV.csv", dir)

write.csv(MUT_TGSA_SANGER, file=file1, row.names=T)
write.csv(EXP_TGSA_SANGER, file=file2, row.names=T)
write.csv(CNV_TGSA_SANGER, file=file3, row.names=T)


# Expression for the cell similarity graph
EXP_Total_Scaled_SANGER = EXP_SANGER[rownames(EXP_TGSA_SANGER), ] %>% scale %>% as.data.frame   # 984 x 37566
EXP_Total_Scaled_SANGER %>% sapply(mean) %>% range %>% round(3)   # [NaN, NaN]
EXP_Total_Scaled_SANGER %>% sapply(sd) %>% range %>% round(3)     # [NA, NA]

col = sapply(EXP_Total_Scaled_SANGER, function(x) sum(is.na(x))==0)   # 36397
EXP_Total_Scaled_SANGER = EXP_Total_Scaled_SANGER[, col]   # 37566 > 36397 [-1169]
EXP_Total_Scaled_SANGER %>% sapply(mean) %>% range %>% round(3)   # [0, 0]
EXP_Total_Scaled_SANGER %>% sapply(sd) %>% range %>% round(3)     # [1, 1]

dir = "../../do_better/TGSA/_data_sanger"
file = sprintf("%s/EXP_Total_Scaled.csv", dir)
write.csv(EXP_Total_Scaled_SANGER, row.names=T, file=file)

# SANGER vs CCLE
labs = c("Expression [CCLE]", "Expression [SANGER]")
main = "Distribution of EXP [TGSA, CCLE & SANGER]"
Dist_EXP_SANGER = create_dist_exp(EXP_TGSA, EXP_TGSA_SANGER)
Dist_EXP_SANGER %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                             alpha=0.5, axis_tl=25, axis_tx=20, 
                             width=width, height=height, dpi=1200, 
                             xy_line=T, raster=T, save=F, save_svg=T)

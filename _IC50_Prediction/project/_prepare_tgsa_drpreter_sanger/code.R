#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
source("../functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)



##### 2. Process data for TGSA

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
geneset_tgsa = gsub("\\(([0-9]*)\\)", "\\1", colnames(EXP_TGSA_Ori)) %>% as.character

EXP_TGSA_Ori = EXP_TGSA_Ori %>% setNames(geneset_tgsa)
MUT_TGSA_Ori = MUT_TGSA_Ori %>% setNames(geneset_tgsa)
CNV_TGSA_Ori = CNV_TGSA_Ori %>% setNames(geneset_tgsa)

Anno_Genes$ENTREZ_ID = Anno_Genes$ENTREZ_ID %>% as.character
sum(geneset_tgsa %in% Anno_Genes$ENTREZ_ID)                            # 706
Anno_Genes_TGSA = Anno_Genes %>% subset(ENTREZ_ID %in% geneset_tgsa)   # 706



### 2. Process Omics Data

dir = "../../processed_data/cell_data/SANGER_Passports"
file_exp = sprintf("%s/TPM_Ent.csv", dir)

dir = "../../raw_data/SANGER_Passports"
file_mut = sprintf("%s/mutations_all_20230202.csv", dir)
file_cnv = sprintf("%s/WES_pureCN_CNV_genes_20221213.csv", dir)

EXP = fread(file_exp, header=T)   # 1432 x 37567
MUT = fread(file_mut)             # 10050692 x 13
CNV = fread(file_cnv)             # 24488104 x 24


# EXP Data
EXP = EXP %>% data.frame(row.names=EXP$V1, check.names=F)
EXP = EXP[, -1]

colnames(EXP) %>% class                # character
all(geneset_tgsa %in% colnames(EXP))   # T
idx = (colnames(EXP) %in% geneset_tgsa)
EXP_TGSA = EXP[, idx]   # 37263 > 706


# MUT Data
# Follow the data processing/filtering in DepMap Omics if possible
# https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=OmicsSomaticMutations.csv
# https://storage.googleapis.com/shared-portal-files/Tools/%5BDMC%20Communication%5D%2022Q4%20Mutation%20Pipeline%20Update.pdf

# Filter out Non-coding mutation
# Filter out Low allele frequency (<0.15)
# SANGER mutation data are already somatic mutations
# https://depmap.sanger.ac.uk/documentation/datasets/mutation

MUT_TGSA = MUT[gene_id %in% Anno_Genes_TGSA$GENE_ID]   # 10050692 > 642036
MUT_TGSA = MUT_TGSA[vaf>=0.15]                         # 642036 > 609340
MUT_TGSA = MUT_TGSA[(coding)]                          # 609340 > 73221
MUT_TGSA$gene_symbol %>% unique %>% length             # 704

# One-hot mutation matrix
any_ = function(x) ifelse(length(x)>=1, 1, 0)
MUT_TGSA = dcast(MUT_TGSA, model_id~gene_id, value.var="gene_id", fun.agg=any_, fill=0)
MUT_TGSA = MUT_TGSA %>% data.frame(row.names=MUT_TGSA$model_id)
MUT_TGSA = MUT_TGSA[, -1]

MUT_TGSA %>% is.na %>% sum      # 0
MUT_TGSA %>% unlist %>% table   # 0, 1 only

all(colnames(MUT_TGSA) %in% Anno_Genes_TGSA$GENE_ID)   # T
idx = match(colnames(MUT_TGSA), Anno_Genes_TGSA$GENE_ID)
col = Anno_Genes_TGSA$ENTREZ_ID[idx]
col %>% is.na %>% sum   # 0
colnames(MUT_TGSA) = col


# CNV Data [CNV Gene Data (WES)]
# Total copy number and categorical CNA calls derived from WES data processed through PureCN
# cell-line & organoid data were processed in SANGER Passports with PureCN [WES] & PURPLE [WGS], respectively
# https://depmap.sanger.ac.uk/documentation/datasets/copy-number
# 
# CCLE DepMap calculates CNV as log2(CN ratio + 1) [OmicsCNGene.csv]
# "Values are calculated by mapping genes onto the segment level calls and computing a weighted average along the genomic coordinate"
# https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=OmicsCNGene.csv
# https://github.com/broadinstitute/depmap_omics/blob/master/depmapomics/copynumbers.py
# 
# In PureCN, CN ratio is given in log2
# Therefore, we converted the CNV into log2(2**gene_mean+1) values (see Table 2)
# 
# Table 2: callAlterations output columns.
# "Gene copy number log2-ratios (not adjusted for purity/ploidy).
# gene.mean is weighted by interval weights when these are available."
# https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/PureCN.pdf

xlab = "Gene Mean\n[gene_mean]"
ylab = "Gene Mean\n[gatk_mean_log2_copy_ratio]"
main = "CNV [706 CGC Genes, Gene Mean & Gene Mean GATK]"
CNV_TGSA = CNV[gene_id %in% Anno_Genes_TGSA$GENE_ID]   # 24488104 > 938543

CNV_TGSA[!is.na(gene_mean) & !is.na(gatk_mean_log2_copy_ratio)] %>% 
  plot_def(gene_mean, gatk_mean_log2_copy_ratio, main=main, size=1.2, 
           alpha=0.5, xlab=xlab, ylab=ylab, xy_line=T, raster=T, save=T)

CNV_TGSA$gene_mean %>% mean(na.rm=T)   # -0.04359442
CNV_TGSA = dcast(CNV_TGSA, model_id~gene_id, value.var="gene_mean", fun.agg=mean, fill=0)
CNV_TGSA = CNV_TGSA %>% data.frame(row.names=CNV_TGSA$model_id)
CNV_TGSA = CNV_TGSA[, -1]   # 1253 x 698

CNV_TGSA %>% is.na %>% sum   # 20094
CNV_TGSA[is.na(CNV_TGSA)] = 0
CNV_TGSA %>% unlist %>% hist
CNV_TGSA %>% unlist %>% mean %>% round(3)   # -0.041

CNV_TGSA = CNV_TGSA %>% apply(2, function(x) log(2**x+1, base=2))
CNV_TGSA %>% unlist %>% hist
CNV_TGSA %>% unlist %>% mean %>% round(3)   # 1.002

idx = match(colnames(CNV_TGSA), Anno_Genes_TGSA$GENE_ID)
colnames(CNV_TGSA) = Anno_Genes_TGSA$ENTREZ_ID[idx]

CNV_TGSA_Ori %>% unlist %>% hist
CNV_TGSA_Ori %>% unlist %>% mean %>% round(3)   # 1.001


# Aggregate Omics [Union]
fill_col_na = function(df, col, value=0) {
  df = df %>% as.data.frame
  col_na = col[!(col %in% colnames(df))]
  sprintf("The number of columns missing : %s", len(col_na)) %>% print
  df[, col_na] = value
  df = df[, col]
  return(df)
}

MUT_TGSA = MUT_TGSA %>% fill_col_na(geneset_tgsa, value=0)   # 705 > 706
EXP_TGSA = EXP_TGSA %>% fill_col_na(geneset_tgsa, value=0)   # 706 > 706
CNV_TGSA = CNV_TGSA %>% fill_col_na(geneset_tgsa, value=1)   # 698 > 706

identical(colnames(MUT_TGSA), colnames(EXP_TGSA))   # T
identical(colnames(MUT_TGSA), colnames(CNV_TGSA))   # T

cells = Reduce(intersect, list(rownames(MUT_TGSA), rownames(EXP_TGSA), rownames(CNV_TGSA)))
MUT_TGSA = MUT_TGSA[cells, ]   # 1357 > 1183
EXP_TGSA = EXP_TGSA[cells, ]   # 1293 > 1183
CNV_TGSA = CNV_TGSA[cells, ]   # 1253 > 1183

dir = mkdir("../../do_better/TGSA_SANGER/_data")
file_mut = sprintf("%s/MUT.csv", dir)
file_exp = sprintf("%s/EXP.csv", dir)
file_cnv = sprintf("%s/CNV.csv", dir)

write.csv(MUT_TGSA, file=file_mut, row.names=T)
write.csv(EXP_TGSA, file=file_exp, row.names=T)
write.csv(CNV_TGSA, file=file_cnv, row.names=T)


# Expression for the cell similarity graph [TGSA only]
EXP_TGSA_Total = EXP[rownames(EXP_TGSA), ] %>% scale %>% as.data.frame
EXP_TGSA_Total %>% sapply(mean) %>% range   # [NaN NaN]
EXP_TGSA_Total %>% sapply(sd) %>% range     # [NA NA]

file_exp_total = sprintf("%s/EXP_Total_Scaled.csv", dir)
col = sapply(EXP_TGSA_Total, function(x) sum(is.na(x))==0)
EXP_TGSA_Total = EXP_TGSA_Total[, col]   # 37566 > 36462
fwrite(EXP_TGSA_Total, row.names=T, file=file_exp_total)


# Confirmed that EXP & CNV are very similar as original [not MUT]
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

identical(rownames(EXP_TGSA_Ori), rownames(MUT_TGSA_Ori))   # T
identical(rownames(EXP_TGSA_Ori), rownames(CNV_TGSA_Ori))   # T
idx = match(rownames(EXP_TGSA_Ori), Anno_Cells$BROAD_ID)
cells_sanger = Anno_Cells$SANGER_MODEL_ID[idx]

rownames(EXP_TGSA_Ori) = cells_sanger
rownames(MUT_TGSA_Ori) = cells_sanger
rownames(CNV_TGSA_Ori) = cells_sanger

Dist_MUT_TGSA = create_dist_exp(MUT_TGSA, MUT_TGSA_Ori, mut=T)
Dist_EXP_TGSA = create_dist_exp(EXP_TGSA, EXP_TGSA_Ori)
Dist_CNV_TGSA = create_dist_exp(CNV_TGSA, CNV_TGSA_Ori)

width = 15
height = 15

omics = c("MUT", "EXP", "CNV")
main = sprintf("Distribution of %s [TGSA_SANGER]", omics)
type = c("SANGER", "CCLE, Github")

labs = sprintf("Mutation [%s]", type)
color = scale_fill_gradient(low="grey90", high="grey90")
Dist_MUT_TGSA %>% grid_def(EXP, EXP_Ori, fill=Freq, main=main[1], color=color,
                           xlab=labs[1], ylab=labs[2], axis_tl=25, axis_tx=20, 
                           margin=0.4, angle=0, width=width, height=height, 
                           mean_summary=F, legend=F, save=T, save_svg=T)

labs = sprintf("Expression [%s]", type)
Dist_EXP_TGSA %>% plot_def(EXP, EXP_Ori, main=main[2], xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20, 
                           width=width, height=height, dpi=1200, 
                           xy_line=T, raster=T, save=T, save_svg=T)

labs = sprintf("CNV [%s]", type)
Dist_CNV_TGSA %>% plot_def(EXP, EXP_Ori, main=main[3], xlab=labs[1], ylab=labs[2],
                           alpha=0.5, axis_tl=25, axis_tx=20, 
                           width=width, height=height, dpi=1200, 
                           xy_line=T, raster=T, save=T, save_svg=T)

rm(EXP, CNV, MUT)



##### 3. Process data for DRPreter

# KEGG genes from the DRPreter original data [2369]
dir = "../../do_better/DRPreter/Data/Cell"
file = sprintf("%s/CCLE_2369_EXP.csv", dir)
EXP_DRPreter_Ori = read.csv(file, row.names=1, check.names=F)
geneset_drpreter = colnames(EXP_DRPreter_Ori)

# SANGER Cell-Passport Data
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Sym.csv", dir)
EXP = fread(file, header=T, check.names=F)
EXP = EXP %>% data.frame(row.names=EXP$V1, check.names=F)
EXP = EXP[, -1]

# Filter out genes not in original data
all(geneset_drpreter %in% colnames(EXP))   # T
EXP_DRPreter = EXP[, geneset_drpreter]     # 2369
# EXP = EXP[, colnames(EXP) %in% geneset_drpreter]

dir = mkdir("../../do_better/DRPreter_SANGER/_data")
file = sprintf("%s/EXP.csv", dir)
write.csv(EXP_DRPreter, file=file, row.names=T)


# # Expression for the cell similarity graph
# EXP_DRPreter_Total = EXP[rownames(EXP_DRPreter), ] %>% scale %>% as.data.frame
# EXP_DRPreter_Total %>% sapply(mean) %>% range   # [NaN NaN]
# EXP_DRPreter_Total %>% sapply(sd) %>% range     # [NA NA]
# 
# file_exp_total = sprintf("%s/EXP_Total_Scaled.csv", dir)
# col = sapply(EXP_DRPreter_Total, function(x) sum(is.na(x))==0)   # 36773 [from 37600]
# EXP_DRPreter_Total = EXP_DRPreter_Total[, col]
# fwrite(EXP_DRPreter_Total, row.names=T, file=file_exp_total)


# Confirmed that values are very similar as original
# Some values can vary due to data processing and/or database update
idx = match(rownames(EXP_DRPreter_Ori), Anno_Cells$BROAD_ID)
cells_sanger = Anno_Cells$SANGER_MODEL_ID[idx]
rownames(EXP_DRPreter_Ori) = cells_sanger
Dist_EXP_DRP = create_dist_exp(EXP_DRPreter, EXP_DRPreter_Ori)

main = "Distribution of EXP [DRPreter_SANGER]"
type = c("SANGER", "CCLE, Github")

labs = sprintf("Expression [%s]", type)
Dist_EXP_DRP %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                          alpha=0.5, axis_tl=25, axis_tx=20, 
                          width=width, height=height, dpi=1200, 
                          xy_line=T, raster=T, save=T, save_svg=T)

rm(EXP)



##### 4. Process GDSC Microarray Data

dir = "../../do_better/HiDRA/_data"
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

# Filter Genes [CGC 706, KEGG 2369]
EXP_TGSA_GDSC = EXP_GDSC_Ent %>% fill_col_na(geneset_tgsa)       # 966 x 706 [NA 41]
EXP_DRPreter_GDSC = EXP_GDSC %>% fill_col_na(geneset_drpreter)   # 968 x 2369 [NA 181]

# Retrieve Cells [Intersecting in MUT, EXP, CNV (TGSA only)]
cells = Reduce(intersect, list(rownames(MUT_TGSA), rownames(EXP_TGSA_GDSC), rownames(CNV_TGSA)))   # 939
EXP_TGSA_GDSC = EXP_TGSA_GDSC[cells, ]
MUT_TGSA_GDSC = MUT_TGSA[cells, ]
CNV_TGSA_GDSC = CNV_TGSA[cells, ]

dir = mkdir("../../do_better/TGSA_SANGER/_data_gdsc")
file1 = sprintf("%s/MUT.csv", dir)
file2 = sprintf("%s/EXP.csv", dir)
file3 = sprintf("%s/CNV.csv", dir)

write.csv(MUT_TGSA_GDSC, file=file1, row.names=T)
write.csv(EXP_TGSA_GDSC, file=file2, row.names=T)
write.csv(CNV_TGSA_GDSC, file=file3, row.names=T)

# dir = mkdir("../../do_better/DRPreter_SANGER/_data_gdsc")
dir = mkdir("../../do_better/DRPreter_SANGER/_data")
file = sprintf("%s/EXP_GDSC.csv", dir)
write.csv(EXP_DRPreter_GDSC, file=file, row.names=T)

# dir = mkdir("../../do_better/DRPreter/_data_gdsc")
dir = mkdir("../../do_better/DRPreter/_data")
file = sprintf("%s/EXP_GDSC.csv", dir)
write.csv(EXP_DRPreter_GDSC, file=file, row.names=T)


# Expression for the cell similarity graph [TGSA only]
EXP_TGSA_GDSC_Total = EXP_GDSC_Ent[rownames(EXP_TGSA_GDSC), ] %>% scale %>% as.data.frame    # 968 x 17419
EXP_TGSA_GDSC_Total %>% sapply(mean) %>% range %>% round(3)   # [0, 0]
EXP_TGSA_GDSC_Total %>% sapply(sd) %>% range %>% round(3)     # [1, 1]

dir = "../../do_better/TGSA_SANGER/_data_gdsc"
file_exp_total = sprintf("%s/EXP_Total_Scaled.csv", dir)
col = sapply(EXP_TGSA_GDSC_Total, function(x) sum(is.na(x))==0)
EXP_TGSA_GDSC_Total = EXP_TGSA_GDSC_Total[, col]   # 16669 > 16669
fwrite(EXP_TGSA_GDSC_Total, row.names=T, file=file_exp_total)

# SANGER vs GDSC [TGSA]
labs = c("Expression [SANGER]", "Expression [GDSC]")
main = "Distribution of EXP [TGSA_SANGER, SANGER & GDSC]"
Dist_EXP_TGSA_GDSC = create_dist_exp(EXP_TGSA, EXP_TGSA_GDSC)
Dist_EXP_TGSA_GDSC %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                                alpha=0.5, axis_tl=25, axis_tx=20,
                                width=width, height=height, dpi=1200,
                                xy_line=T, raster=T, save=T, save_svg=T)

# SANGER vs GDSC [DRPreter]
EXP_DRPreter_GDSC_ = EXP_DRPreter_GDSC
idx = match(rownames(EXP_DRPreter_GDSC_), Anno_Cells$COSMIC_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_DRPreter_GDSC_ = EXP_DRPreter_GDSC_[!is.na(cells), ]
rownames(EXP_DRPreter_GDSC_) = na.omit(cells)

labs = c("Expression [SANGER]", "Expression [GDSC]")
main = "Distribution of EXP [DRPreter_SANGER, SANGER & GDSC]"
Dist_EXP_DRP_GDSC = create_dist_exp(EXP_DRPreter, EXP_DRPreter_GDSC_)
Dist_EXP_DRP_GDSC %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                               alpha=0.5, axis_tl=25, axis_tx=20,
                               width=width, height=height, dpi=1200,
                               xy_line=T, raster=T, save=T, save_svg=T)



##### 5. Process CCLE TPM Data [TGSA only]

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

# Filter Genes [CGC 706]
EXP_TGSA_CCLE = EXP_CCLE %>% fill_col_na(geneset_tgsa)   # 1310 x 706

# Retrieve Cells [Intersecting in MUT, EXP, CNV]
cells = Reduce(intersect, list(rownames(MUT_TGSA), rownames(EXP_TGSA_CCLE), rownames(CNV_TGSA)))   # 910
EXP_TGSA_CCLE = EXP_TGSA_CCLE[cells, colnames(MUT_TGSA_SANGER)]
MUT_TGSA_CCLE = MUT_TGSA[cells, ]
CNV_TGSA_CCLE = CNV_TGSA[cells, ]

identical(colnames(EXP_TGSA_CCLE), colnames(MUT_TGSA_CCLE))   # T
identical(colnames(EXP_TGSA_CCLE), colnames(CNV_TGSA_CCLE))   # T

dir = mkdir("../../do_better/TGSA_SANGER/_data_ccle")
file1 = sprintf("%s/MUT.csv", dir)
file2 = sprintf("%s/EXP.csv", dir)
file3 = sprintf("%s/CNV.csv", dir)

write.csv(MUT_TGSA_CCLE, file=file1, row.names=T)
write.csv(EXP_TGSA_CCLE, file=file2, row.names=T)
write.csv(CNV_TGSA_CCLE, file=file3, row.names=T)


# Expression for the cell similarity graph
EXP_TGSA_Total_CCLE = EXP_CCLE[rownames(EXP_TGSA_CCLE), ] %>% scale %>% as.data.frame   # 910 x 19193
EXP_TGSA_Total_CCLE %>% sapply(mean) %>% range %>% round(3)   # [NaN, NaN]
EXP_TGSA_Total_CCLE %>% sapply(sd) %>% range %>% round(3)     # [NA, NA]

col = sapply(EXP_TGSA_Total_CCLE, function(x) sum(is.na(x))==0)   # 19189
EXP_TGSA_Total_CCLE = EXP_TGSA_Total_CCLE[, col]   # 19193 > 19189 [-4]

dir = "../../do_better/TGSA_SANGER/_data_ccle"
file = sprintf("%s/EXP_Total_Scaled.csv", dir)
write.csv(EXP_TGSA_Total_CCLE, row.names=T, file=file)

# SANGER vs CCLE [TGSA]
labs = c("Expression [SANGER]", "Expression [CCLE]")
main = "Distribution of EXP [TGSA_SANGER, SANGER & CCLE]"
Dist_EXP_TGSA_CCLE = create_dist_exp(EXP_TGSA, EXP_TGSA_CCLE)
Dist_EXP_TGSA_CCLE %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                                alpha=0.5, axis_tl=25, axis_tx=20,
                                width=width, height=height, dpi=1200,
                                xy_line=T, raster=T, save=T, save_svg=T)


# SANGER vs CCLE [DRPreter]

dir = "../../do_better/DRPreter/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_DRP_CCLE = read.csv(file, row.names=1)   # 1389 x 2369 [-140]

# Rename Cells [BROAD > SANGER]
idx = match(rownames(EXP_DRP_CCLE), Anno_Cells$BROAD_ID)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
EXP_DRP_CCLE = EXP_DRP_CCLE[!is.na(cells), ]
rownames(EXP_DRP_CCLE) = na.omit(cells)   # 1304 x 2369 [-140]

labs = c("Expression [SANGER]", "Expression [CCLE]")
main = "Distribution of EXP [DRPreter_SANGER, SANGER & CCLE]"
Dist_EXP_DRP_CCLE = create_dist_exp(EXP_DRPreter, EXP_DRP_CCLE)
Dist_EXP_DRP_CCLE %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                                alpha=0.5, axis_tl=25, axis_tx=20,
                                width=width, height=height, dpi=1200,
                                xy_line=T, raster=T, save=T, save_svg=T)

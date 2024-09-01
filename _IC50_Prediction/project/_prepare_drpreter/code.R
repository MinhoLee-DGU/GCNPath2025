#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
source("../functions.R")
loadings()

# Gene Annotation from SANGER Cell Passport
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)



##### 2. Process data [DRPreter]

# KEGG genes from the DRPreter original data [2369]
dir = "../../do_better/DRPreter/Data/Cell"
file = sprintf("%s/CCLE_2369_EXP.csv", dir)
EXP_DRPreter_Ori = read.csv(file, row.names=1, check.names=F)
geneset_ori = EXP_DRPreter_Ori %>% colnames %>% as.character

# # CCLE TPM data [DepMap 21Q4]
# dir = "../../raw_data/CCLE"
# file = sprintf("%s/CCLE_expression.csv", dir)
# EXP = fread_def(file, check_names=F)   # 1389 x 19177
# 
# col_exp = colnames(EXP) %>% col_to_symbol(geneset_ori)   # 19177 [from 19177]
# EXP_Sym = EXP
# colnames(EXP_Sym) = col_exp
# 
# all(geneset_ori %in% colnames(EXP_Sym))   # F
# (genes_no = setdiff(geneset_ori, colnames(EXP_Sym)))
# # CASP12  MICOS10-NBL1  UBE2NL  BUB1B-PAK6

col_to_symbol = function(col, col_interest=NULL) {
  col = gsub("(.*)\\((.*)\\)", "\\1", col) %>% as.character
  col = col %>% trimws
  sprintf("Unique genes : %s [from %s]", length(unique(col)), length(col)) %>% print
  
  if (!is.null(col_interest)) {
    col_dup = col[duplicated(col)] %>% unique
    col_dup_in = intersect(col_interest, col_dup)
    sprintf("Any column interested is duplicated? [%s]", length(col_dup_in)) %>% print
    if (length(col_dup_in)!=0) sprintf("# %s", paste(col_dup_in, collapse=",")) %>% print
  }
  
  return(col)
}

# CCLE TPM data [DepMap 21Q4]
dir = "../../raw_data/CCLE"
file = sprintf("%s/CCLE_expression_full.csv", dir)
EXP = fread_def(file, check_names=F)   # 1389 x 52063

EXP_Total = EXP
col_exp_total = colnames(EXP_Total) %>% col_to_symbol
EXP_Total = EXP_Total[, !duplicated(col_exp_total)]
colnames(EXP_Total) = col_exp_total[!duplicated(col_exp_total)]   # 1389 x 52054

col_exp = colnames(EXP) %>% col_to_symbol(geneset_ori)   # 52054 [from 52063]
colnames(EXP) = col_exp
EXP_Sym = EXP

all(geneset_ori %in% colnames(EXP_Sym))   # T
EXP_DRPreter = EXP_Sym[, geneset_ori]     # 52063 > 2369
# EXP_DRPreter = EXP_Sym[, colnames(EXP_Sym) %in% geneset_ori]

dir = mkdir("../../do_better/DRPreter/_data")
file_exp = sprintf("%s/EXP.csv", dir)
write.csv(EXP_DRPreter, row.names=T, file=file_exp)

file_exp_total = sprintf("%s/EXP_Total.csv", dir)
write.csv(EXP_Total, row.names=T, file=file_exp_total)


# # Expression for the cell similarity graph
# file_exp_scaled = sprintf("%s/EXP_Scaled.csv", dir)
# EXP_DRPreter_Scaled = EXP_DRPreter %>% scale %>% as.data.frame
# EXP_DRPreter_Scaled %>% sapply(mean) %>% range   # [-1.241301e-15  1.216673e-15]
# EXP_DRPreter_Scaled %>% sapply(sd) %>% range     # [1, 1]
# write.csv(EXP_DRPreter_Scaled, row.names=T, file=file_exp_scaled)

# Confirmed that values are the same as original
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

Dist_EXP = create_dist_exp(EXP_DRPreter, EXP_DRPreter_Ori)

width = 15
height = 15

omics = c("MUT", "EXP", "CNV")
main = "Distribution of EXP [DRPreter]"
type = c("CCLE", "CCLE, Github")

labs = sprintf("Expression [%s]", type)
Dist_EXP %>% plot_def(EXP, EXP_Ori, main=main, xlab=labs[1], ylab=labs[2],
                      alpha=0.5, axis_tl=25, axis_tx=20, 
                      width=width, height=height, dpi=1200, 
                      xy_line=T, raster=T, save=T, save_svg=T)

# cells_int = intersect(rownames(EXP_DRPreter), rownames(EXP_DRPreter_Ori))   # 580
# diff_exp_drpreter = (EXP_DRPreter[cells_int, ]-EXP_DRPreter_Ori[cells_int, ]) %>% unlist
# range(diff_exp_drpreter)   # -1.776357e-15  8.881784e-16
# 
# main = sprintf("%s/Difference in EXP [DRPreter]", dir)
# diff_exp_drpreter %>% hist_def(main=main, width=16.5, height=15, dist=F, save=T)



# ##### 3. Process data [TGSA]
# 
# # CGC Genes from the TGSA original data [706]
# file_exp = "exp.csv"
# file_mut = "mu.csv"
# file_cnv = "cn.csv"
# 
# EXP_TGSA_Ori = read.csv(file_exp, row.names=1, check.names=F)
# MUT_TGSA_Ori = read.csv(file_mut, row.names=1, check.names=F)
# CNV_TGSA_Ori = read.csv(file_cnv, row.names=1, check.names=F)
# 
# identical(colnames(EXP_TGSA_Ori), colnames(MUT_TGSA_Ori))   # T
# identical(colnames(EXP_TGSA_Ori), colnames(CNV_TGSA_Ori))   # T
# geneset_tgsa = gsub("\\(([0-9]*)\\)", "\\1", colnames(EXP_TGSA_Ori)) %>% as.character
# 
# Anno_Genes$ENTREZ_ID = Anno_Genes$ENTREZ_ID %>% as.character
# sum(geneset_tgsa %in% Anno_Genes$ENTREZ_ID)                            # 706
# Anno_Genes_TGSA = Anno_Genes %>% subset(ENTREZ_ID %in% geneset_tgsa)   # 706
# 
# colnames(EXP_TGSA_Ori) = geneset_tgsa
# colnames(MUT_TGSA_Ori) = geneset_tgsa
# colnames(CNV_TGSA_Ori) = geneset_tgsa
# 
# 
# dir = "../../raw_data/CCLE"
# file_exp = sprintf("%s/CCLE_expression.csv", dir)
# 
# file_mut1 = sprintf("%s/CCLE_mutations_bool_hotspot.csv", dir)
# file_mut2 = sprintf("%s/CCLE_mutations_bool_damaging.csv", dir)
# file_mut3 = sprintf("%s/CCLE_mutations_bool_nonconserving.csv", dir)
# file_mut4 = sprintf("%s/CCLE_mutations_bool_otherconserving.csv", dir)
# 
# file_cnv1 = sprintf("%s/CCLE_gene_cn.csv", dir)
# file_cnv2 = sprintf("%s/CCLE_wes_gene_cn.csv", dir)
# 
# EXP = fread(file_exp, na.strings="")     # 1450 x 19194
# MUT1 = fread(file_mut1, na.strings="")   # 1408099 x 56
# MUT2 = fread(file_mut2, na.strings="")   # 1408099 x 56
# MUT3 = fread(file_mut3, na.strings="")   # 1408099 x 56
# MUT4 = fread(file_mut4, na.strings="")   # 1408099 x 56
# CNV1 = fread(file_cnv1, na.strings="")   # 1814 x 25369
# CNV2 = fread(file_cnv2, na.strings="")   # 1814 x 25369
# 
# to_data_frame = function(df) {
#   df = df %>% data.frame(row.names=df$V1, check.names=F)
#   df = df[, -1]
#   return(df)
# }
# 
# to_entrez = function(Omics) {
#   colnames(Omics) = gsub(".*\\(([0-9]*)\\)", "\\1", colnames(Omics)) %>% as.character
#   return(Omics)
# }
# EXP_Ent = EXP %>% to_data_frame %>% to_entrez
# MUT1_Ent = MUT1 %>% to_data_frame %>% to_entrez
# MUT2_Ent = MUT2 %>% to_data_frame %>% to_entrez
# MUT3_Ent = MUT3 %>% to_data_frame %>% to_entrez
# MUT4_Ent = MUT4 %>% to_data_frame %>% to_entrez
# CNV1_Ent = CNV1 %>% to_data_frame %>% to_entrez
# CNV2_Ent = CNV2 %>% to_data_frame %>% to_entrez
# 
# all(geneset_tgsa %in% colnames(EXP_Ent))    # T
# all(geneset_tgsa %in% colnames(MUT1_Ent))   # F
# all(geneset_tgsa %in% colnames(MUT2_Ent))   # F
# all(geneset_tgsa %in% colnames(MUT3_Ent))   # T
# all(geneset_tgsa %in% colnames(MUT4_Ent))   # F
# all(geneset_tgsa %in% colnames(CNV1_Ent))   # T
# all(geneset_tgsa %in% colnames(CNV2_Ent))   # T
# 
# EXP_Ent = EXP_Ent[, colnames(EXP_Ent) %in% geneset_tgsa]
# MUT1_Ent = MUT1_Ent[, colnames(MUT1_Ent) %in% geneset_tgsa]   # 461 
# MUT2_Ent = MUT2_Ent[, colnames(MUT2_Ent) %in% geneset_tgsa]   # 689
# MUT3_Ent = MUT3_Ent[, colnames(MUT3_Ent) %in% geneset_tgsa]   # 706
# MUT4_Ent = MUT4_Ent[, colnames(MUT4_Ent) %in% geneset_tgsa]   # 20
# CNV1_Ent = CNV1_Ent[, colnames(CNV1_Ent) %in% geneset_tgsa]   # 
# CNV2_Ent = CNV2_Ent[, colnames(CNV2_Ent) %in% geneset_tgsa]   # 
# 
# 
# 
# # Confirmed that values are very similar as original
# # Some values can vary due to data processing and/or database update
# 
# confirm_mut = function(MUT, MUT_Ori) {
#   cells_int = intersect(rownames(MUT), rownames(MUT_Ori))
#   genes_int = intersect(colnames(MUT), colnames(MUT_Ori))
#   
#   lvl = c(0, 1)
#   mut = MUT[cells_int, genes_int] %>% unlist %>% factor(levels=lvl)
#   mut_ori = MUT_Ori[cells_int, genes_int] %>% unlist %>% factor(levels=lvl)
#   
#   mut_table = caret::confusionMatrix(mut, mut_ori)
#   rownames(mut_table$table) = sprintf("Latest_%s", c(0, 1))
#   colnames(mut_table$table) = sprintf("Github_%s", c(0, 1))
#   
#   sprintf("Cells : %s", length(cells_int)) %>% print
#   sprintf("Genes : %s", length(genes_int)) %>% print
#   return(mut_table$table)
# }
# 
# confirm_exp = function(EXP, EXP_Ori) {
#   cells_int = intersect(rownames(EXP), rownames(EXP_Ori))
#   genes_int = intersect(colnames(EXP), colnames(EXP_Ori))
#   exp_diff = (EXP[cells_int, genes_int]-EXP_Ori[cells_int, genes_int]) %>% unlist
#   exp_diff %>% hist_def(save=F)
#   exp_diff %>% range(na.rm=T) %>% print
# }
# 
# confirm_exp(EXP_Ent, EXP_TGSA_Ori)
# confirm_exp(CNV1_Ent, CNV_TGSA_Ori)
# confirm_exp(CNV2_Ent, CNV_TGSA_Ori)
# confirm_mut(MUT1_Ent, MUT_TGSA_Ori)   # 
# confirm_mut(MUT2_Ent, MUT_TGSA_Ori)   # 
# confirm_mut(MUT3_Ent, MUT_TGSA_Ori)   # 
# confirm_mut(MUT4_Ent, MUT_TGSA_Ori)   # 

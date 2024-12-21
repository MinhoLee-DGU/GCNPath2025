#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(openxlsx))
source("../functions.R")
loadings()

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/Anno_Drugs.csv", dir)
Anno_Drugs = read.csv(file)   # 700
Anno_Drugs = Anno_Drugs %>% subset(!is.na(Drug_CID))  # 597

dir = "../../processed_data/cell_data/GDSC"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells_GDSC = read.csv(file)   # 1939
Anno_Cells_GDSC$SANGER_MODEL_ID %>% is.na %>% sum   # 0

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)   # 2099
Anno_Cells$SANGER_MODEL_ID %>% is.na %>% sum   # 0

dir = "../../raw_data/GDSC"
file1 = sprintf("%s/PANCANCER_IC_Thu Jun  8 07_43_49 2023.csv", dir)
file2 = sprintf("%s/PANCANCER_IC_Thu Jun  8 07_44_45 2023.csv", dir)
IC50_GDSC1 = read.csv(file1) %>% as.data.frame  # 333292
IC50_GDSC2 = read.csv(file2) %>% as.data.frame  # 243466

IC50_GDSC1$Dataset.Version %>% table   # GDSC1
IC50_GDSC2$Dataset.Version %>% table   # GDSC2
colnames(IC50_GDSC1) = gsub("\\.", "_", colnames(IC50_GDSC1)) %>% toupper
colnames(IC50_GDSC2) = gsub("\\.", "_", colnames(IC50_GDSC2)) %>% toupper


# Match SANGER_MODEL_ID
idx1 = match(IC50_GDSC1$COSMIC_ID, Anno_Cells_GDSC$COSMIC_ID)
idx2 = match(IC50_GDSC2$COSMIC_ID, Anno_Cells_GDSC$COSMIC_ID)
sanger_id1 = Anno_Cells_GDSC$SANGER_MODEL_ID[idx1]
sanger_id2 = Anno_Cells_GDSC$SANGER_MODEL_ID[idx2]

sanger_id1 %>% is.na %>% sum   # 0
sanger_id2 %>% is.na %>% sum   # 0
IC50_GDSC1$SANGER_MODEL_ID = sanger_id1
IC50_GDSC2$SANGER_MODEL_ID = sanger_id2

idx1 = match(IC50_GDSC1$SANGER_MODEL_ID, Anno_Cells$SANGER_MODEL_ID)
idx2 = match(IC50_GDSC2$SANGER_MODEL_ID, Anno_Cells$SANGER_MODEL_ID)
broad_id1 = Anno_Cells$BROAD_ID[idx1]
broad_id2 = Anno_Cells$BROAD_ID[idx2]

broad_id1 %>% is.na %>% sum   # 0
broad_id2 %>% is.na %>% sum   # 0
IC50_GDSC1$BROAD_ID = broad_id1
IC50_GDSC2$BROAD_ID = broad_id2

IC50_GDSC1 = IC50_GDSC1 %>% 
  relocate(SANGER_MODEL_ID, BROAD_ID, COSMIC_ID, .before=CELL_LINE_NAME)
IC50_GDSC2 = IC50_GDSC2 %>% 
  relocate(SANGER_MODEL_ID, BROAD_ID, COSMIC_ID, .before=CELL_LINE_NAME)

# SANGER_MODEL_ID can be unique identifiers
# CELL_LINE_NAME & SANGER_MODEL_ID is not dulplicated (1:1 corresponding)

IC50_GDSC1$SANGER_MODEL_ID %>% unique %>% length   # 970
IC50_GDSC1$CELL_LINE_NAME %>% unique %>% length    # 970
IC50_GDSC2$SANGER_MODEL_ID %>% unique %>% length   # 969
IC50_GDSC2$CELL_LINE_NAME %>% unique %>% length    # 969


# Match DRUG_CID
idx1_1 = match(IC50_GDSC1$DRUG_ID, Anno_Drugs$Drug_ID)
idx2_1 = match(IC50_GDSC2$DRUG_ID, Anno_Drugs$Drug_ID)
drug_cid1_1 = Anno_Drugs$Drug_CID[idx1_1]
drug_cid2_1 = Anno_Drugs$Drug_CID[idx2_1]

idx1_2 = match(IC50_GDSC1$DRUG_NAME, Anno_Drugs$Name)
idx2_2 = match(IC50_GDSC2$DRUG_NAME, Anno_Drugs$Name)
drug_cid1_2 = Anno_Drugs$Drug_CID[idx1_2]
drug_cid2_2 = Anno_Drugs$Drug_CID[idx2_2]

drug_cid1 = ifelse(!is.na(drug_cid1_1), drug_cid1_1, drug_cid1_2)
drug_cid2 = ifelse(!is.na(drug_cid2_1), drug_cid2_1, drug_cid2_2)
IC50_GDSC1$DRUG_CID = drug_cid1
IC50_GDSC2$DRUG_CID = drug_cid2

IC50_GDSC1$IC50 %>% range   # -10.57774   12.35455
IC50_GDSC2$IC50 %>% range   # -8.769011   13.847363

IC50_GDSC1 = IC50_GDSC1 %>% mutate(LN_IC50=IC50) %>% 
  subset(select=-IC50) %>% relocate(LN_IC50, .before=AUC) %>% 
  relocate(SANGER_MODEL_ID, BROAD_ID, COSMIC_ID, CELL_LINE_NAME, 
           DRUG_CID, DRUG_NAME, DRUG_ID, .before=everything())

IC50_GDSC2 = IC50_GDSC2 %>% mutate(LN_IC50=IC50) %>% 
  subset(select=-IC50) %>% relocate(LN_IC50, .before=AUC) %>% 
  relocate(SANGER_MODEL_ID, BROAD_ID, COSMIC_ID, CELL_LINE_NAME, 
           DRUG_CID, DRUG_NAME, DRUG_ID, .before=everything())


# Filter IC50s without corresponding PubChem CIDs
IC50_GDSC1 = IC50_GDSC1 %>% subset(!is.na(DRUG_CID))    # 333292 > 286620
IC50_GDSC2 = IC50_GDSC2 %>% subset(!is.na(DRUG_CID))    # 243466 > 207350
IC50_GDSC1$TAG = paste0(IC50_GDSC1$SANGER_MODEL_ID, "@", IC50_GDSC1$DRUG_CID)
IC50_GDSC2$TAG = paste0(IC50_GDSC2$SANGER_MODEL_ID, "@", IC50_GDSC2$DRUG_CID)

rm(sanger_id1, broad_id1, drug_cid1, drug_cid1_1, drug_cid1_2, 
   sanger_id2, broad_id2, drug_cid2, drug_cid2_1, drug_cid2_2)



##### 2. Remove Duplicated IC50

average_dup = function(IC50) {
  
  n_before = IC50 %>% nrow
  num_tags = IC50$TAG %>% table
  dup_tags = names(num_tags)[num_tags>1]
  
  col = c("SANGER_MODEL_ID", "DRUG_CID", "LN_IC50", "AUC", "TAG")
  IC50_Dup_Raw = IC50[, col] %>% subset(TAG %in% dup_tags)
  IC50_Dup = IC50_Dup_Raw %>% group_by(TAG) %>% 
    summarise(LN_IC50_Mean=mean(LN_IC50), AUC_Mean=mean(AUC), N_Dup=n()) %>% as.data.frame
  
  idx = match(IC50_Dup_Raw$TAG, IC50_Dup$TAG)
  IC50_Dup_Raw = IC50_Dup_Raw %>% 
    mutate(LN_IC50_Mean=IC50_Dup$LN_IC50_Mean[idx], AUC_Mean=IC50_Dup$AUC_Mean[idx]) %>% 
    relocate(TAG, .after=everything())
  
  idx = match(IC50$TAG, IC50_Dup$TAG)
  IC50$AUC = ifelse(is.na(idx), IC50$AUC, IC50_Dup$AUC_Mean[idx])
  IC50$LN_IC50 = ifelse(is.na(idx), IC50$LN_IC50, IC50_Dup$LN_IC50_Mean[idx])
  IC50 = IC50 %>% subset(!duplicated(TAG))
  
  n_after = IC50 %>% nrow
  sprintf("# Remove Duplicates : %s > %s", n_before, n_after) %>% print
  IC50_List = list(IC50=IC50, IC50_Dup=IC50_Dup, IC50_Dup_Raw=IC50_Dup_Raw)
  return(IC50_List)
}

# Duplicates Within GDSC1 & GDSC2
IC50_GDSC1_List = IC50_GDSC1 %>% average_dup   # 286620 > 266573 [264328 in Release 8.2]
IC50_GDSC2_List = IC50_GDSC2 %>% average_dup   # 207350 > 199687 [120930 in Release 8.2]

# Duplicates Between GDSC1 & GDSC2
# Merge two datasets into IC50_GDSC
identical(colnames(IC50_GDSC1_List$IC50), colnames(IC50_GDSC2_List$IC50))   # T
IC50_GDSC = rbind(IC50_GDSC1_List$IC50, IC50_GDSC2_List$IC50)   # 466260 [266573+199687]
IC50_GDSC_List = IC50_GDSC %>% average_dup   # 466260 > 373681 [335035 in Release 8.2]



##### 3. Analysis of IC50 Duplicates

# IC50 data have their own variation
# It means that we can not theoretically decrease the prediction RMSE below the IC50 variation RMSE
# If we can make a very accurate model... it means just overfitting to train data... so bad!!!

IC50_GDSC1_List$IC50_Dup$N_Dup %>% table   # 2 [18415], 3 [816] 
IC50_GDSC2_List$IC50_Dup$N_Dup %>% table   # 2 [7663]
IC50_GDSC_List$IC50_Dup$N_Dup %>% table    # 2 [92579]

# RMSE between GDSC1 & GDSC2 increased into 0.944
rmse1 = IC50_GDSC1_List$IC50_Dup_Raw %>% with(RMSE(LN_IC50, LN_IC50_Mean))   # 0.8540166
rmse2 = IC50_GDSC2_List$IC50_Dup_Raw %>% with(RMSE(LN_IC50, LN_IC50_Mean))   # 0.8957159
rmse12 = IC50_GDSC_List$IC50_Dup_Raw %>% with(RMSE(LN_IC50, LN_IC50_Mean))   # 0.9449771

corr1 = IC50_GDSC1_List$IC50_Dup_Raw %>% with(cor(LN_IC50, LN_IC50_Mean))    # 0.9464328
corr2 = IC50_GDSC2_List$IC50_Dup_Raw %>% with(cor(LN_IC50, LN_IC50_Mean))    # 0.9617034
corr12 = IC50_GDSC_List$IC50_Dup_Raw %>% with(cor(LN_IC50, LN_IC50_Mean))    # 0.9384612

dir = mkdir("../../processed_data/ic50_data/GDSC")
file = sprintf("%s/LN_IC50 & LN_IC50_Mean [GDSC%s]", dir, c(1, 2, ""))

xlab = bquote(ln(IC["50"]))
ylab = bquote(Mean~ln(IC["50"]))
# xlab = bquote(bold(ln(IC["50"])))
# ylab = bquote(bold(Mean~ln(IC["50"])))

IC50_GDSC1_List$IC50_Dup_Raw %>% 
  plot_def(LN_IC50, LN_IC50_Mean, main=file[1], 
           xlab=xlab, ylab=ylab, size=1.5, alpha=0.5, axis_tl=30, axis_tx=25, 
           dpi=1200, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T)
IC50_GDSC2_List$IC50_Dup_Raw %>% 
  plot_def(LN_IC50, LN_IC50_Mean, main=file[2], 
           xlab=xlab, ylab=ylab, size=1.5, alpha=0.5, axis_tl=30, axis_tx=25,  
           dpi=1200, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T)
IC50_GDSC_List$IC50_Dup_Raw %>% 
  plot_def(LN_IC50, LN_IC50_Mean, main=file[3], 
           xlab=xlab, ylab=ylab, size=1.5, alpha=0.5, axis_tl=30, axis_tx=25, 
           dpi=1200, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T)


dir = mkdir("../../processed_data/ic50_data/GDSC")
file = sprintf("%s/AUC & AUC_Mean [GDSC%s]", dir, c(1, 2, ""))

xlab = "AUC"
ylab = "Mean AUC"
# xlab = bquote(bold(AUC))
# ylab = bquote(bold(Mean~AUC))

IC50_GDSC1_List$IC50_Dup_Raw %>% 
  plot_def(AUC, AUC_Mean, main=file[1], alpha=0.5, xlab=xlab, ylab=ylab, size=1.5,
           axis_tl=30, axis_tx=25, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T, dpi=1500)
IC50_GDSC2_List$IC50_Dup_Raw %>% 
  plot_def(AUC, AUC_Mean, main=file[2], alpha=0.5, xlab=xlab, ylab=ylab, size=1.5,
           axis_tl=30, axis_tx=25, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T, dpi=1500)
IC50_GDSC_List$IC50_Dup_Raw %>% 
  plot_def(AUC, AUC_Mean, main=file[3], alpha=0.5, xlab=xlab, ylab=ylab, size=1.5,
           axis_tl=30, axis_tx=25, xy_line=T, force_bold=F, raster=T, save=T, save_svg=T, dpi=1500)



##### 4. Save Data

# IC50 Data [Full Information]
IC50_GDSC1 = IC50_GDSC1_List$IC50
IC50_GDSC2 = IC50_GDSC2_List$IC50
IC50_GDSC = IC50_GDSC_List$IC50

IC50_GDSC1$SANGER_MODEL_ID %>% unique %>% length   # 970
IC50_GDSC2$SANGER_MODEL_ID %>% unique %>% length   # 969
IC50_GDSC$SANGER_MODEL_ID %>% unique %>% length    # 978
IC50_GDSC$SANGER_MODEL_ID %>% table %>% hist       # [23, 432] 

IC50_GDSC1$DRUG_CID %>% unique %>% length   # 321
IC50_GDSC2$DRUG_CID %>% unique %>% length   # 236
IC50_GDSC$DRUG_CID %>% unique %>% length    # 432
IC50_GDSC$DRUG_CID %>% table %>% hist       # [99,979]

dir = "../../processed_data/ic50_data/GDSC"
file = sprintf("%s/IC50_GDSC%s.csv", dir, c(1, 2, ""))
write.csv(IC50_GDSC1, file[1], row.names=F)
write.csv(IC50_GDSC2, file[2], row.names=F)
write.csv(IC50_GDSC, file[3], row.names=F)

file = sprintf("%s/IC50_GDSC%s.xlsx", dir, c(1, 2, ""))
write.xlsx(IC50_GDSC1_List, file=file[1], rowNames=F)
write.xlsx(IC50_GDSC2_List, file=file[2], rowNames=F)
write.xlsx(IC50_GDSC_List, file=file[3], rowNames=F)


# IC50 Data [Simplified Information]
col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "DRUG_CID", "LN_IC50")
IC50_GDSC1_Brief = IC50_GDSC1[, col]
IC50_GDSC2_Brief = IC50_GDSC2[, col]
IC50_GDSC_Brief = IC50_GDSC[, col]

col = c("Cell", "Cell_BROAD", "Cell_COSMIC", "Drug", "LN_IC50")
colnames(IC50_GDSC1_Brief) = col
colnames(IC50_GDSC2_Brief) = col
colnames(IC50_GDSC_Brief) = col

file = sprintf("%s/IC50_GDSC%s.txt", dir, c(1, 2, ""))
write.table(IC50_GDSC1_Brief, file[1], sep="\t", row.names=F)
write.table(IC50_GDSC2_Brief, file[2], sep="\t", row.names=F)
write.table(IC50_GDSC_Brief, file[3], sep="\t", row.names=F)

file = sprintf("%s/IC50_GDSC.RData", dir)
save(IC50_GDSC1, IC50_GDSC2, IC50_GDSC, 
     IC50_GDSC1_List, IC50_GDSC2_List, IC50_GDSC_List, file=file)

rm(IC50_GDSC1_Brief, IC50_GDSC2_Brief, IC50_GDSC_Brief)



supplementary = T
if (supplementary) {
  suppressMessages(library(openxlsx))
  
  ### [Supplementary Data] Supplementary Data 4-6
  col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "CELL_LINE_NAME", 
          "DRUG_CID", "DRUG_NAME", "DRUG_ID", "LN_IC50", "DATASET_VERSION")
  
  IC50_GDSC1_ = IC50_GDSC1[, col]
  IC50_GDSC2_ = IC50_GDSC2[, col]
  IC50_GDSC_ = IC50_GDSC[, col]
  
  sheets = "Supplementary Data 4"
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(IC50_GDSC1_, file=file, sheetName=sheets, rowNames=F)
  
  sheets = "Supplementary Data 5"
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(IC50_GDSC2_, file=file, sheetName=sheets, rowNames=F)
  
  sheets = "Supplementary Data 6"
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(IC50_GDSC_, file=file, sheetName=sheets, rowNames=F)
  
  
  ### [Source Data] Supplementary Fig. 37
  IC50_GDSC1_Dup_ = IC50_GDSC1_List$IC50_Dup_Raw %>% subset(select=-c(TAG, AUC, AUC_Mean))
  IC50_GDSC2_Dup_ = IC50_GDSC2_List$IC50_Dup_Raw %>% subset(select=-c(TAG, AUC, AUC_Mean))
  IC50_GDSC_Dup_ = IC50_GDSC_List$IC50_Dup_Raw %>% subset(select=-c(TAG, AUC, AUC_Mean))
  df_list = list(IC50_GDSC1_Dup_, IC50_GDSC2_Dup_, IC50_GDSC_Dup_)
  
  num_fig = letters[1:3]
  sheets = sprintf("Supplementary Fig. 37%s", num_fig)
  file = "SourceData_SupplementaryFig37.xlsx"
  write.xlsx(df_list, file=file, sheetName=sheets, rowNames=F)
}


# ##### 5. SANGER Passport or DepMap CCLE?
# # Confirmed the SANGER Passport guarantees more IC50 available than DepMap CCLE 
# 
# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/TPM_Ent.csv", dir)
# SANGER_TPM = read.csv(file, row.names=1, check.names=F)   # 1431 x 37566
# 
# dir = "../../processed_data/cell_data/CCLE_DepMap"
# file = sprintf("%s/TPM_Ent.csv", dir)
# CCLE_TPM = read.csv(file, row.names=1, check.names=F)     # 1310 x 18782
# 
# cells_sanger = rownames(SANGER_TPM)
# cells_ccle = rownames(CCLE_TPM)
# cells_ic50 = IC50_GDSC$SANGER_MODEL_ID %>% unique
# 
# intersect(cells_ic50, cells_sanger) %>% length   # 972
# intersect(cells_ic50, cells_ccle) %>% length     # 703
# 
# IC50_GDSC %>% subset(SANGER_MODEL_ID %in% cells_sanger) %>% nrow    # 371751 [from 373681]
# IC50_GDSC1 %>% subset(SANGER_MODEL_ID %in% cells_sanger) %>% nrow   # 265448 [from 266573]
# IC50_GDSC2 %>% subset(SANGER_MODEL_ID %in% cells_sanger) %>% nrow   # 198551 [from 199687]
# 
# IC50_GDSC %>% subset(SANGER_MODEL_ID %in% cells_ccle) %>% nrow    # 270717 [from 373681]
# IC50_GDSC1 %>% subset(SANGER_MODEL_ID %in% cells_ccle) %>% nrow   # 192186 [from 266573]
# IC50_GDSC2 %>% subset(SANGER_MODEL_ID %in% cells_ccle) %>% nrow   # 145421 [from 199687]

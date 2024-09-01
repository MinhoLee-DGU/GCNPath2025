#!/usr/bin/env Rscript

##### 1. Packages and Data

source("../functions.R")
loadings()

dir = "../../processed_data/ic50_data/GDSC"
file = sprintf("%s/IC50_GDSC.csv", dir)
IC50_GDSC = read.csv(file)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells_SANGER = read.csv(file)

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/Anno_Drugs.csv", dir)
Anno_Drugs_GDSC = read.csv(file)

dir = "../../raw_data/CCLE_DepMap"
file = sprintf("%s/CCLE_NP24.2009_Drug_data_2015.02.24.csv", dir)
IC50_CCLE = read.csv(file, check.names=F)

dir = "../../processed_data/cell_data/CCLE_DepMap"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells_CCLE = read.csv(file, na.strings=c("", "NA"))



##### 2. Process CCLE IC50 Data

symbols = c(" ", "\\(", "\\)")
colnames(IC50_CCLE) = colnames(IC50_CCLE) %>% standard_cols(symbols=symbols)

IC50_CCLE = IC50_CCLE %>% 
  dplyr::rename(Cell=CCLE_Cell_Line_Name, Cell_Primary=Primary_Cell_Line_Name, 
                Drug=Compound, Doses_uM=Doses_UM, EC50_uM=EC50_UM, IC50_uM=IC50_UM) %>% 
  mutate(LN_IC50=log(IC50_uM)) %>% relocate(LN_IC50, IC50_uM, EC50_uM, .after=Drug) %>% as.data.frame

ic50_ccle_min = sub(",.*", "", IC50_CCLE$Doses) %>% as.numeric
ic50_ccle_max = sub(".*,", "", IC50_CCLE$Doses) %>% as.numeric

sum(IC50_CCLE$IC50_uM==ic50_ccle_min)   # 48 [from 11670]
sum(IC50_CCLE$IC50_uM==ic50_ccle_max)   # 6451 [from 11670]
capping = (IC50_CCLE$IC50_uM==ic50_ccle_min) | (IC50_CCLE$IC50_uM==ic50_ccle_max)   # 6499 [6451+48]
IC50_CCLE = IC50_CCLE %>% mutate(Capping=capping) %>% relocate(Capping, .after=EC50_uM)



##### 3. Process CCLE Drug Annotation

cids_first_match = function(Anno_CIDs, verbose_dup=T) {
  Anno_CIDs_ = Anno_CIDs %>% subset(!is.na(Drug_CID)) %>% distinct
  drugs_dup = Anno_CIDs_$Query[duplicated(Anno_CIDs_$Query)] %>% unique
  Anno_CIDs_Dup = Anno_CIDs_[Anno_CIDs_$Query %in% drugs_dup, ]
  
  if (verbose_dup) {
    if (nrow(Anno_CIDs_Dup)!=0) {
      # Anno_CIDs_Dup %>% View
      Anno_CIDs_Dup %>% print
    } else print("No drugs with multiple CIDs...")
  }
  
  Anno_CIDs_First = Anno_CIDs_ %>% distinct(Query, .keep_all=T)
  return(Anno_CIDs_First)
}

suppressMessages(library(webchem))
drug_ccle = IC50_CCLE$Drug %>% unique   # 24
Anno_Drugs_CCLE = get_cid(drug_ccle, from="name", domain="compound", match="all", verbose=T) 

Anno_Drugs_CCLE = Anno_Drugs_CCLE %>% as.data.frame
colnames(Anno_Drugs_CCLE) = c("Query", "Drug_CID")
Anno_Drugs_CCLE = cids_first_match(Anno_Drugs_CCLE)
Anno_Drugs_CCLE = Anno_Drugs_CCLE %>% rename(Drug_Name=Query)
# Confirmed drugs with 1:multi match have similar structures... 
# [17-AAG, Nutlin-3, TKI258]

# If drugs exist in GDSC also, use those CIDs...
col = c("Name", "Synonyms", "Drug_CID")
Anno_Drugs_GDSC = Anno_Drugs_GDSC[, col] %>% 
  subset(!is.na(Drug_CID)) %>% 
  separate_rows(Name, sep=",|, ") %>% 
  separate_rows(Synonyms, sep=",|, ") %>% 
  mutate(Name=trimws(Name), Synonyms=trimws(Synonyms)) %>% as.data.frame
  
idx1 = match(tolower(Anno_Drugs_CCLE$Drug_Name), tolower(Anno_Drugs_GDSC$Name))
idx2 = match(tolower(Anno_Drugs_CCLE$Drug_Name), tolower(Anno_Drugs_GDSC$Synonyms))
Anno_Drugs_CCLE$GDSC = (!is.na(idx1) | !is.na(idx2))

Anno_Drugs_CCLE$GDSC %>% sum   # 17
Anno_Drugs_CCLE = Anno_Drugs_CCLE %>% 
  mutate(Drug_CID=ifelse(is.na(idx1), Drug_CID, Anno_Drugs_GDSC$Drug_CID[idx1])) %>% 
  mutate(Drug_CID=ifelse(is.na(idx2), Drug_CID, Anno_Drugs_GDSC$Drug_CID[idx2]))

# Manual check for drugs not in GDSC but CCLE
sum(Anno_Drugs_CCLE$Drug_CID %in% Anno_Drugs_GDSC$Drug_CID)   # 17
Anno_Drugs_CCLE %>% subset(!GDSC) %>% pull(Drug_Name)
# AEW541  Nutlin-3  L-685458  ZD-6474  LBW242  RAF265  TKI258
# Anno_Drugs_GDSC %>% View

# Nutlin-3 [CCLE, CID 216345] & Nutlin-3a (-) [GDSC, CID 11433190] both exists...
Anno_Drugs_CCLE$GDSC[Anno_Drugs_CCLE$Drug_Name=="Nutlin-3"] = T
Anno_Drugs_CCLE$Drug_CID[Anno_Drugs_CCLE$Drug_Name=="Nutlin-3"] = 11433190

Anno_Drugs_CCLE$GDSC %>% sum   # 18
sum(Anno_Drugs_CCLE$Drug_CID %in% Anno_Drugs_GDSC$Drug_CID)   # 18

# Add Taget information from IC50_CCLE
col = c("Drug", "Target")
DTI_CCLE = IC50_CCLE[, col] %>% distinct
DTI_CCLE = DTI_CCLE %>% group_by(Drug) %>% 
  summarise(Target=paste0(Target, collapse="|")) %>% as.data.frame

idx = match(Anno_Drugs_CCLE$Drug_Name, DTI_CCLE$Drug)
Anno_Drugs_CCLE = Anno_Drugs_CCLE %>% 
  mutate(Target=DTI_CCLE$Target[idx]) %>% 
  relocate(Target, .before=GDSC) %>% as.data.frame


# Get SMILES using CIDs from PubChem Homepage
# https://pubchem.ncbi.nlm.nih.gov
dir = mkdir("../../processed_data/drug_data/CCLE")
file = sprintf("%s/CID_List.txt", dir)
cid_list = Anno_Drugs_CCLE$Drug_CID   # 24
cat(cid_list, file=file, sep="\n")

# Get SMILES of PubChem CIDs
file = sprintf("%s/PubChem_SMILES.csv", dir)
SMILES_CCLE = read.csv(file)

col = c("cid", "inchikey", "canonicalsmiles", "isosmiles")
SMILES_CCLE = SMILES_CCLE[, col]   # 417
colnames(SMILES_CCLE) = c("Drug_CID", "INCHI_Key", "SMILES_CAN", "SMILES_ISO")

file = sprintf("%s/SMILES_CCLE.csv", dir)
write.csv(SMILES_CCLE, file=file, row.names=F)



##### 4. Comparison of IC50 Data between GDSC and CCLE

# Complete BROAD_ID and SANGER_ID in CCLE IC50 Data
idx1 = match(IC50_CCLE$Cell, Anno_Cells_CCLE$CCLE_ID)
idx2 = match(IC50_CCLE$Cell_Primary, Anno_Cells_CCLE$MODEL_NAME)
broad_ccle1 = Anno_Cells_CCLE$BROAD_ID[coalesce(idx1, idx2)]
sanger_ccle1 = Anno_Cells_CCLE$SANGER_MODEL_ID[coalesce(idx1, idx2)]

idx1 = match(IC50_CCLE$Cell, Anno_Cells_SANGER$CCLE_ID)
idx2 = match(IC50_CCLE$Cell_Primary, Anno_Cells_SANGER$MODEL_NAME)
broad_ccle2 = Anno_Cells_SANGER$BROAD_ID[coalesce(idx1, idx2)]
sanger_ccle2 = Anno_Cells_SANGER$SANGER_MODEL_ID[coalesce(idx1, idx2)]

IC50_CCLE = IC50_CCLE %>% 
  mutate(Cell_BROAD_ID=coalesce(broad_ccle1, broad_ccle2), 
         Cell_SANGER_ID=coalesce(sanger_ccle1, sanger_ccle2)) %>% 
  relocate(Cell_BROAD_ID, Cell_SANGER_ID, .after=Cell_Primary) %>% as.data.frame

# # Find missing BROAD_ID and SANGER_ID by manual googling... [X]
IC50_CCLE$Cell_BROAD_ID %>% is.na %>% sum    # 48
IC50_CCLE$Cell_SANGER_ID %>% is.na %>% sum   # 94
all(which(is.na(IC50_CCLE$Cell_BROAD_ID)) %in% which(is.na(IC50_CCLE$Cell_SANGER_ID)))   # T

# IC50_CCLE %>% subset(is.na(Cell_BROAD_ID)) %>% pull(Cell) %>% unique
# # COLO677_LUNG, MB157_BREAST
# IC50_CCLE %>% subset(is.na(Cell_BROAD_ID)) %>% pull(Cell_Primary) %>% unique
# # COLO-677, MB 157

# Complete PubChem CID in CCLE IC50 Data
idx = match(IC50_CCLE$Drug, Anno_Drugs_CCLE$Drug_Name)
IC50_CCLE = IC50_CCLE %>% 
  mutate(Drug_CID=as.character(Anno_Drugs_CCLE$Drug_CID[idx])) %>% 
  relocate(Drug_CID, .after=Drug) %>% as.data.frame


# Comparison of LN_IC50
# This plot will be drawn again with prediction errors of GCNPath
# Just save the object [IC50_GDSC_CCLE] for later use

col1 = c("SANGER_MODEL_ID", "BROAD_ID", "CELL_LINE_NAME", "DRUG_NAME", "DRUG_CID", "LN_IC50")
col2 = c("Cell_BROAD_ID", "Drug_CID", "LN_IC50", "Capping")
by = c("BROAD_ID"="Cell_BROAD_ID", "DRUG_CID"="Drug_CID")

IC50_GDSC$DRUG_CID = IC50_GDSC$DRUG_CID %>% as.character
IC50_CCLE$Drug_CID = IC50_CCLE$Drug_CID %>% as.character
IC50_GDSC_CCLE = inner_join(IC50_GDSC[, col1], IC50_CCLE[, col2], by=by)   # 5553

IC50_GDSC_CCLE = IC50_GDSC_CCLE %>% 
  dplyr::rename(LN_IC50_GSDC=LN_IC50.x, LN_IC50_CCLE=LN_IC50.y) %>% as.data.frame
IC50_GDSC_CCLE %>% plot_def(LN_IC50_GSDC, LN_IC50_CCLE, xy_line=T, force_bold=T, save=F)
IC50_GDSC_CCLE[idx, ] %>% plot_def(LN_IC50_GSDC, LN_IC50_CCLE, xy_line=T, force_bold=T, save=F)

idx = !(IC50_GDSC_CCLE$Capping)   # 2555
IC50_GDSC_CCLE %>% with(cor(LN_IC50_GSDC, LN_IC50_CCLE)) %>% round(3)           # 0.788
IC50_GDSC_CCLE %>% with(RMSE(LN_IC50_GSDC, LN_IC50_CCLE)) %>% round(3)          # 2.04
IC50_GDSC_CCLE[idx, ] %>% with(cor(LN_IC50_GSDC, LN_IC50_CCLE)) %>% round(3)    # 0.726
IC50_GDSC_CCLE[idx, ] %>% with(RMSE(LN_IC50_GSDC, LN_IC50_CCLE)) %>% round(3)   # 2.215



##### 5. Save files...

dir = mkdir("../../processed_data/ic50_data/CCLE")
file = sprintf("%s/IC50_CCLE.csv", dir)
write.csv(IC50_CCLE, file=file, row.names=F)

col = c("Cell_Primary", "Cell_SANGER_ID", "Cell_BROAD_ID", "Drug_CID", "Capping", "LN_IC50")
file = sprintf("%s/IC50_CCLE.txt", dir)
write.table(IC50_CCLE[, col], file=file, row.names=F, sep="\t")

dir = "../../processed_data/drug_data/CCLE"
file = sprintf("%s/Anno_Drugs.csv", dir)
write.csv(Anno_Drugs_CCLE, file=file, row.names=F)

dir = "../../processed_data/ic50_data/GDSC"
file = sprintf("%s/IC50_GDSC_CCLE.csv", dir)
write.csv(IC50_GDSC_CCLE, file=file, row.names=F)

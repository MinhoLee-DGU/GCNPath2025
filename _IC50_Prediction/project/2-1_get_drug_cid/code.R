#!/usr/bin/env Rscript

##### 1. Packages and Data

source("../functions.R")
loadings()

dir = "../../raw_data/GDSC"
file = sprintf("%s/Drug_listThu Jun  8 07_40_56 2023.csv", dir)
Anno_Drugs = read.csv(file, header=T, na.strings="")
colnames(Anno_Drugs) = colnames(Anno_Drugs) %>% standard_cols(specific="ID")

Anno_Drugs = Anno_Drugs %>% 
  mutate(Drug_CID=PubCHEM) %>% 
  relocate(Drug_CID, .after=Synonyms) %>% 
  subset(select=-c(PubCHEM, Number_Of_Cell_Lines))

Anno_Drugs$Drug_CID %>% is.na %>% sum    # 240 [from 700]
# Anno_Drugs$Drug_CID %>% table
# none, None, several, multiple CIDs with comma

idx = grep(",", Anno_Drugs$Drug_CID)
Anno_Drugs[idx, c("Name", "Drug_CID")]

# All Have almost the same structures in Drug_CID
# Name          Drug_CID
# AZD5363       25227436, 42602260
# Buparlisib    66577015, 16654980
# Ulixertinib   11719003, 58641927
# Ulixertinib   11719003, 58641927

cids_na = c("none", "None", "several")
Anno_Drugs$Drug_CID[idx] = gsub(",.*", "", Anno_Drugs$Drug_CID[idx])
Anno_Drugs$Drug_CID[Anno_Drugs$Drug_CID %in% cids_na] = NA
Anno_Drugs$Drug_CID = Anno_Drugs$Drug_CID %>% as.numeric

Anno_Drugs$Name %>% length             # 700
Anno_Drugs$Name %>% unique %>% length  # 544



# 2. Get CIDs manually from drug names and synonyms by the below
# https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

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

col = c("Drug_ID", "Name", "Synonyms", "Targets", "Drug_CID")
Drugs_No_CIDs = Anno_Drugs[, col] %>% subset(is.na(Drug_CID)) %>% 
  separate_rows(Name, sep=", |; ") %>% separate_rows(Synonyms, sep=", |; ") %>% as.data.frame

dir = mkdir("../../processed_data/drug_data/GDSC")
file1 = sprintf("%s/Drugs_No_CIDs_Names.txt", dir)
file2 = sprintf("%s/Drugs_No_CIDs_Synonyms.txt", dir)

no_cids_names = Drugs_No_CIDs$Name %>% na.omit %>% as.character
no_cids_syns = Drugs_No_CIDs$Synonyms %>% na.omit %>% as.character

cat(no_cids_names, file=file1, sep="\n")
cat(no_cids_syns, file=file2, sep="\n")

# Search by manual...
file1 = sprintf("%s/PubChem_CIDs_Names.txt", dir)
file2 = sprintf("%s/PubChem_CIDs_Synonyms.txt", dir)
PubChem_Names = read.csv(file1, sep="\t", header=F)   # 388
PubChem_Syns = read.csv(file2, sep="\t", header=F)    # 166

col = c("Query", "Drug_CID")
colnames(PubChem_Names) = col
colnames(PubChem_Syns) = col

PubChem_Names_ = PubChem_Names %>% cids_first_match   # 138
PubChem_Syns_ = PubChem_Syns %>% cids_first_match     # 51

# Names are epistasis than synonyms
# Mismatch in Drugs [Synonyms]
# alpha-lipoic acid   aLA   5950
# aLA is mistaken as Alanine...

idx1 = match(Drugs_No_CIDs$Name, PubChem_Names_$Query)
idx2 = match(Drugs_No_CIDs$Synonyms, PubChem_Syns_$Query)

cids1 = PubChem_Names_$Drug_CID[idx1]
cids2 = PubChem_Syns_$Drug_CID[idx2]
Drugs_No_CIDs = Drugs_No_CIDs %>% mutate(Drug_CID = cids1) %>% 
  mutate(Drug_CID = ifelse(!is.na(Drug_CID), Drug_CID, cids2))
# alpha-lipoic acid is corrected as 864 [not 5950]

Drugs_No_CIDs$Source_CID = NA
source = "PubChem [ID Exchange]"
idx = !is.na(Drugs_No_CIDs$Drug_CID)
Drugs_No_CIDs$Source_CID[idx] = source



##### 4. Search by WebChem

suppressMessages(library(webchem))
idx = Drugs_No_CIDs$Drug_CID %>% is.na
no_cids_names = Drugs_No_CIDs$Name[idx] %>% na.omit %>% as.character
no_cids_syns = Drugs_No_CIDs$Synonyms[idx] %>% na.omit %>% as.character

# PubChem [4>3, 0]
WebChem_PC_Names = get_cid(no_cids_names, from="name", domain="compound", match="all", verbose=T) 
WebChem_PC_Syns = get_cid(no_cids_syns, from="name", domain="compound", match="all", verbose=T) 

# Found Just 3...
# Mirin         [CID] 137319711
# CAP-232       [CID] 6918265
# ICL-SIRT078   [CID] 45168044

WebChem_PC_Names = WebChem_PC_Names %>% 
  setNames(c("Query", "Drug_CID")) %>% as.data.frame

WebChem_CIDs = WebChem_PC_Names %>% 
  subset(!is.na(Drug_CID)) %>% distinct %>% as.data.frame

source = "PubChem [WebChem]"
idx = match(Drugs_No_CIDs$Name, WebChem_CIDs$Name)
Drugs_No_CIDs$Drug_CID[!is.na(idx)] = WebChem_CIDs$Drug_CID[na.omit(idx)]
Drugs_No_CIDs$Source_CID[!is.na(idx)] = source



##### 5. Search by Files

match_to_db = function(query, DataBase, key, sep=NULL) {
  
  key_lower = sprintf("%s_", key)
  DataBase[[key_lower]] = DataBase[[key]] %>% tolower
  Query_List = data.frame(Query=query, Query_Lower=tolower(query))
  
  if (!is.null(sep)) {
    DataBase = DataBase %>% separate_rows(all_of(key_lower), sep=sep)
  }
  
  DataBase = DataBase[DataBase[[key_lower]] %in% Query_List$Query_Lower, ]
  idx = match(DataBase[[key_lower]], Query_List$Query_Lower)
  
  DataBase = DataBase[!is.na(idx), ] %>% 
    mutate(Query=Query_List$Query[na.omit(idx)]) %>% 
    relocate(Query, .before=everything()) %>% as.data.frame
  
  DataBase[[key_lower]] = NULL
  sprintf("Query Found : %s", nrow(DataBase)) %>% print
  return(DataBase)
}

dir = "../../raw_data/Drug_List"
file1 = sprintf("%s/drugbank vocabulary.csv", dir)
file2 = sprintf("%s/small_molecule_20230609151933.csv", dir)
file3 = sprintf("%s/DOWNLOAD-bGxOa7WMXUE0NN6ZKzAhPRjpjNgO-KtkZjIlyPHO01c=.csv", dir)

DrugBank = read.csv(file1, header=T, na.strings="")
HMS_LINCS = read.csv(file2, header=T, na.strings="")
ChEMBl = read.csv(file3, sep=";", header=T, na.strings="")

idx = Drugs_No_CIDs$Drug_CID %>% is.na
no_cids_names = Drugs_No_CIDs$Name[idx] %>% na.omit %>% as.character
no_cids_syns = Drugs_No_CIDs$Synonyms[idx] %>% na.omit %>% as.character

LINCS_Names1 = match_to_db(no_cids_names, HMS_LINCS, key="Name")
LINCS_Names2 = match_to_db(no_cids_names, HMS_LINCS, key="Alternative.Names", sep="; ")
LINCS_Syns1 = match_to_db(no_cids_syns, HMS_LINCS, key="Name")
LINCS_Syns2 = match_to_db(no_cids_syns, HMS_LINCS, key="Alternative.Names", sep="; ")
# 16, 0, 0, 0

ChEMBl_Names1 = match_to_db(no_cids_names, ChEMBl, key="Name")
ChEMBl_Names2 = match_to_db(no_cids_names, ChEMBl, key="Synonyms", sep="|")
ChEMBl_Syns1 = match_to_db(no_cids_syns, ChEMBl, key="Name")
ChEMBl_Syns2 = match_to_db(no_cids_syns, ChEMBl, key="Synonyms", sep="|")
# 0, 0, 0, 0

DrugBank_Names1 = match_to_db(no_cids_names, DrugBank, key="Common.name")
DrugBank_Names2 = match_to_db(no_cids_names, DrugBank, key="Synonyms", sep=" | ")
DrugBank_Syns1 = match_to_db(no_cids_syns, DrugBank, key="Common.name")
DrugBank_Syns2 = match_to_db(no_cids_syns, DrugBank, key="Synonyms", sep=" | ")
# 0, 0, 0, 0

##### Results
# ChEMBl      [0]
# HMS-LINCS   [16]
# DrugBank    [0]

# Drugs without CIDs in HMS-LINCS
cids = c(58525121, 73265214, 71576693, 42628535, 68234129, 73265211)
names = c("QL-VIII-58", "QL-XI-92", "QL-XII-61", "XMD11-85h", "XMD13-2", "XMD14-99")
df_temp = data.frame(Name=names, PubChem.CID=cids)

idx = match(LINCS_Names1$Name, df_temp$Name)
LINCS_Names1$PubChem.CID[!is.na(idx)] = df_temp$PubChem.CID[na.omit(idx)]
LINCS_Names1_ = LINCS_Names1 %>% subset(!is.na(PubChem.CID))   # 16 > 14
LINCS_Names1_$Query %>% unique %>% length   # 16 [from 16, no duplicate]

source = "HMS-LINCS"
idx = match(Drugs_No_CIDs$Name, LINCS_Names1$Name)
Drugs_No_CIDs$Drug_CID[!is.na(idx)] = LINCS_Names1$PubChem.CID[na.omit(idx)]
Drugs_No_CIDs$Source_CID[!is.na(idx)] = source



##### Last Manual Search 
##### Googling with or without Targets and/or Synonyms

# ascorbate (vitamin C) [CID 54670067]
# 
# BDP-00009066 [CID 132275018] 
# [BDP9066, Target MRCKB]
# www.rcsb.org/structure/5OTF
# doi: 10.1158/0008-5472.CAN-17-2870
# 
# CRT0105446 [CID 117918919]
# [CRT-0105446]
# https://medkoo.com/products/17685
# doi: 10.18632/oncotarget.6288
# 
# GSK-LSD1-2HCl [CID 71522234]
# Same as GSK-LSD1 [CID 71522234]
# 
# HG-6-71-01 [CID 71550931] 
# [HG-6-71-1, DDR1-IN-2, Target DDR1 & DDR2 & SRC]
# doi.org/10.1016/j.chembiol.2020.07.014
# doi.org/10.1021/cb400430t
# 
# JQ12 [CID 6918878]
# [BRD-6929, Target HDAC1 & HDAC2]
# www.medkoo.com/products/48225
# www.caymanchem.com/product/34149
# 
# LIMK1 inhibitor BMS4 [CID 46192505]
# [LIMK1 inhibitor BMS 4, Target LIMK1]
# www.axonmedchem.com/product/1949
# www.chemexpress.cn/905298-84-2.htm
# www.medchemexpress.com/limk1-inhibitor-bms-4.html?locale=ko-KR
# www.bocsci.com/limk1-inhibitor-bms-4-cas-905298-84-2-item-55022.html
#
# HKMTI-1-005 [CID] 53315436
# www.chemspider.com/Chemical-Structure.26645415.html

##### Synonyms or Mis-Written Names
# [Bleomycin → Bleomycin (50 uM), CID 5460769]
# [Bleomycin (10 uM) → Bleomycin (50 uM), CID 5460769]
# [Veneclexta → Venotoclax, Target BCL-2, CID 49846579]
# [Dimethyloxalylglcine → Dimethyloxalylglycine, CID 560326]
# [LIMK1 inhibitor BMS4 → LIMK1 inhibitor BMS 4, CID 46192505]

names = c("ascorbate (vitamin C)", "BDP-00009066", "Bleomycin", 
          "Bleomycin (10 uM)", "CRT0105446", "GSK-LSD1-2HCl ", 
          "HG-6-71-01", "JQ12", "LIMK1 inhibitor BMS4", "Venotoclax", "HKMTI-1-005")
cids = c(54670067, 132275018, 5460769, 5460769, 117918919, 
         71522234, 71550931, 6918878, 46192505, 49846579, 53315436)
Manual_CIDs = data.frame(Name=names, Drug_CID=cids)

source = "Manual Search"
idx = match(Drugs_No_CIDs$Name, Manual_CIDs$Name)
Drugs_No_CIDs$Drug_CID[!is.na(idx)] = Manual_CIDs$Drug_CID[na.omit(idx)]
Drugs_No_CIDs$Source_CID[!is.na(idx)] = source
# Drugs_No_CIDs %>% subset(Name %in% Manual_CIDs$Name)


# Paste all CIDs into Anno_Drugs
Drugs_No_CIDs_ = Drugs_No_CIDs %>% subset(!is.na(Drug_CID))
idx = match(Anno_Drugs$Drug_ID, Drugs_No_CIDs_$Drug_ID)

cids = Anno_Drugs$Drug_CID
Anno_Drugs$Drug_CID = ifelse(!is.na(cids), cids, Drugs_No_CIDs_$Drug_CID[idx])
Anno_Drugs$Source_CID = ifelse(!is.na(cids), "GDSC", Drugs_No_CIDs_$Source_CID[idx])

# Anno_Drugs %>% subset(Name %in% Manual_CIDs$Name) %>% View
# Anno_Drugs %>% subset(Name %in% Drugs_No_CIDs_$Name) %>% View



##### 6. Relation [DRUG-NAME & CID]

# Catch Drugs [1-DRUG-NAME ~ multi-CID]
names_cids_list = split(Anno_Drugs$Drug_CID, Anno_Drugs$Name) %>% lapply(unique)
check_multi_cids = names_cids_list %>% sapply(length)
drugs_multi_cids = check_multi_cids[check_multi_cids==2] %>% names

names_cids_list[drugs_multi_cids]
col = c("Name", "Drug_CID", "Source_CID")
Anno_Drugs[Anno_Drugs$Name %in% drugs_multi_cids, col]

# 5 drugs have multi-cids
# BMS-345541    [GDSC] 9813758    [GDSC] 9926054        Similar structure
# Cisplatin     [GDSC] 84691      [PubChem] 5702198     Similar structure
# JNK-9L        [GDSC] 25222038   [PubChem] 146156476   Different structure
# MIM1          [GDSC] 16241412   [PubChem] 135691163   Different structure
# Oxaliplatin   [GDSC] 5310940    [PubChem] 9887053     Similar structure

# Rules of CID choice
# 1. For no other reasons, Follow GDSC CIDs
# 2. If GDSC CIDs becomes legacy, Follow Manual CIDs (MIM1)
# 3. Cisplatin CID should be set as 5702198, not 84691

# GCN cannot process CID 84691 as all the atoms are separated as ions
# In General, hydrogens [H] are not explicitly processed
# RDKit and SMILES cannot recognize which atom pairs form ionic bonds
# [CID 5702198] N.N.Cl[Pt]Cl
# [CID 84691] N.N.[Cl-].[Cl-].[Pt+2]

Anno_Drugs$Drug_CID[Anno_Drugs$Name=="BMS-345541"] = 9813758
Anno_Drugs$Drug_CID[Anno_Drugs$Name=="Cisplatin"] = 5702198
Anno_Drugs$Drug_CID[Anno_Drugs$Name=="JNK-9L"] = 25222038
Anno_Drugs$Drug_CID[Anno_Drugs$Name=="MIM1"] = 135691163
Anno_Drugs$Drug_CID[Anno_Drugs$Name=="Oxaliplatin"] = 5310940

Anno_Drugs$Source_CID[Anno_Drugs$Name=="BMS-345541"] = "GDSC"
Anno_Drugs$Source_CID[Anno_Drugs$Name=="Cisplatin"] = "PubChem [ID Exchange]"
Anno_Drugs$Source_CID[Anno_Drugs$Name=="JNK-9L"] = "GDSC"
Anno_Drugs$Source_CID[Anno_Drugs$Name=="MIM1"] = "PubChem [ID Exchange]"
Anno_Drugs$Source_CID[Anno_Drugs$Name=="Oxaliplatin"] = "GDSC"


# Get SMILES using CIDs from PubChem Homepage
# https://pubchem.ncbi.nlm.nih.gov

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/CID_List.txt", dir)
cid_list = Anno_Drugs$Drug_CID %>% na.omit %>% unique   # 432
cat(cid_list, file=file, sep="\n")

# Get SMILES of PubChem CIDs
file = sprintf("%s/PubChem_SMILES.csv", dir)
SMILES_GDSC = read.csv(file)

col = c("cid", "inchikey", "canonicalsmiles", "isosmiles")
SMILES_GDSC = SMILES_GDSC[, col]   # 417
cids_legacy = setdiff(cid_list, SMILES_GDSC$cid)   # 17
names_legacy = Anno_Drugs$Name[Anno_Drugs$Drug_CID %in% cids_legacy] %>% unique

# Get CIDs of Legacy Drugs
PubChem_Legacy = get_cid(names_legacy, from="name", domain="compound", match="all", verbose=T)
PubChem_Legacy$query %>% unique %>% length   # 17 [from 17]
PubChem_Legacy$cid[PubChem_Legacy$query=="KIN001-260"] = 135465539

# KIN001-260 [CID 135465539]
# https://lincs.hms.harvard.edu/db/sm/10197-102-1/

colnames(PubChem_Legacy) = c("Query", "Drug_CID_New")
PubChem_Legacy = PubChem_Legacy %>% as.data.frame

idx = match(PubChem_Legacy$Query, Anno_Drugs$Name)
PubChem_Legacy$Drug_CID_Legacy = Anno_Drugs$Drug_CID[idx]
all(Anno_Drugs$Source_CID[idx]=="GDSC")   # T

row = c("MIM1", 135691163, 16241412)
PubChem_Legacy = PubChem_Legacy %>% rbind(row)

# Get SMILES of PubChem CIDs [Legacy Drugs]
file = sprintf("%s/Drug_List_Legacy.csv", dir)
write.csv(PubChem_Legacy, file=file, row.names=F)

file = sprintf("%s/PubChem_SMILES_Legacy.csv", dir)
SMILES_GDSC_Legacy = read.csv(file)
all(SMILES_GDSC_Legacy$cid %in% PubChem_Legacy$Drug_CID_New)   # T

col = c("cid", "inchikey", "canonicalsmiles", "isosmiles")
SMILES_GDSC_Legacy = SMILES_GDSC_Legacy[, col]
SMILES_GDSC = SMILES_GDSC %>% rbind(SMILES_GDSC_Legacy) %>% distinct
SMILES_GDSC$cid %>% unique %>% length   # 432
colnames(SMILES_GDSC) = c("Drug_CID", "INCHI_Key", "SMILES_CAN", "SMILES_ISO")

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/SMILES_GDSC.csv", dir)
write.csv(SMILES_GDSC, file=file, row.names=F)


# Correct legacy CIDs from Anno_Drugs
idx = match(Anno_Drugs$Drug_CID, PubChem_Legacy$Drug_CID_Legacy)
Anno_Drugs$Drug_CID[!is.na(idx)] = PubChem_Legacy$Drug_CID_New[na.omit(idx)]
Anno_Drugs$Drug_CID %>% na.omit %>% unique %>% length         # 432
all(na.omit(Anno_Drugs$Drug_CID) %in% SMILES_GDSC$Drug_CID)   # T

source_legacy = sprintf("GDSC [Legacy CID %s]", PubChem_Legacy$Drug_CID_Legacy)
Anno_Drugs$Source_CID[!is.na(idx)] = source_legacy[na.omit(idx)]

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/Anno_Drugs.csv", dir)
write.csv(Anno_Drugs, file=file, row.names=F)


supplementary = T
if (supplementary) {
  suppressMessages(library(openxlsx))
  
  ### [Supplementary Data] Supplementary Data 3
  Anno_Drugs_ = Anno_Drugs %>% 
    subset(!is.na(Drug_CID)) %>% 
    arrange(Drug_ID) %>% as.data.frame   # 597
  
  col = c("Drug_CID", "InCHI_Key", "Canonical_SMILES", "Isomeric_SMILES")
  SMILES_GDSC_ = SMILES_GDSC %>% setNames(col)
  
  Anno_Drugs_$Drug_CID = Anno_Drugs_$Drug_CID %>% as.character
  SMILES_GDSC_$Drug_CID = SMILES_GDSC_$Drug_CID %>% as.character
  Anno_Drugs_ = Anno_Drugs_ %>% left_join(SMILES_GDSC_, by="Drug_CID")   # 597 x 12
  
  
  # Manual Search
  drug_names = c("ascorbate (vitamin C)", "BDP-00009066", "Bleomycin", 
                 "Bleomycin (10 uM)", "CRT0105446", "GSK-LSD1-2HCl ", 
                 "HG-6-71-01", "JQ12", "LIMK1 inhibitor BMS4", "Venotoclax", "HKMTI-1-005")
  
  sources = c(
    "Equal to Vitamin C", 
    "doi.org/10.1158/0008-5472.CAN-17-2870|www.rcsb.org/structure/5OTF", 
    "Equal to Bleomycin (50 uM)", 
    "Equal to Bleomycin (50 uM)",
    "doi.org/10.18632/oncotarget.6288|www.medkoo.com/products/17685",
    "Equal to GSK-LSD1",
    "doi.org/10.1021/cb400430t|doi.org/10.1016/j.chembiol.2020.07.014", 
    "www.caymanchem.com/product/34149|www.medkoo.com/products/48225",
    "www.axonmedchem.com/product/1949|www.chemexpress.cn/905298-84-2.htm|www.medchemexpress.com/limk1-inhibitor-bms-4.html?locale=ko-KR|www.bocsci.com/limk1-inhibitor-bms-4-cas-905298-84-2-item-55022.html", 
    "Equal to Venetoclax", 
    "www.chemspider.com/Chemical-Structure.26645415.html"
  )
  
  Anno_Drugs_Manual = data.frame(Name=drug_names, Source=sources)
  idx = match(Anno_Drugs_$Name, Anno_Drugs_Manual$Name)
  Anno_Drugs_$Source_CID_Manual = Anno_Drugs_Manual$Source[idx]   # 597 x 13
  
  dir = "../../processed_data/drug_data/GDSC"
  sheets = "Supplementary Data 3"
  file = sprintf("%s/%s.xlsx", dir, sheets)
  write.xlsx(Anno_Drugs_, file=file, sheetName=sheets, rowNames=F)
}

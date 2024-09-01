#!/usr/bin/env Rscript

##### 1. Packages and Data

source("../functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

# ### Test cells included in not IC50 but RNA-Seq [Cell-Blind like]
# 
# file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
# SANGER_RNA_GSVA = read.csv(file, row.names=1)
# 
# dir = "../../processed_data/ic50_data/GDSC"
# file = sprintf("%s/IC50_GDSC.txt", dir)
# IC50_GDSC = read.csv(file, sep="\t")
# 
# IC50_GDSC_Num = IC50_GDSC %>% group_by(Cell) %>% 
#   summarise(Num_IC50=n()) %>% as.data.frame
# idx = match(Anno_Cells$SANGER_MODEL_ID, IC50_GDSC_Num$Cell)
# 
# Anno_Cells = Anno_Cells %>% 
#   mutate(RNA=SANGER_MODEL_ID %in% rownames(SANGER_RNA_GSVA), 
#          IC50=SANGER_MODEL_ID %in% unique(IC50_GDSC$Cell), 
#          Num_IC50=IC50_GDSC_Num$Num_IC50[idx])
# 
# Anno_Cells %>% group_by(TCGA_CODE) %>%
#   summarise(Num_Screen_O=sum(RNA & IC50),
#             Num_Screen_X=sum(RNA & !IC50)) %>%
#   arrange(desc(Num_Screen_X)) %>% as.data.frame %>% View

standard_cols = function(x, symbols=c("\\.", "-", " ", "\\(", "\\)", "#")) {
  
  if (!is.null(symbols)) {
    symbols = paste0(symbols, collapse="|")
    x = gsub(symbols, "_", x)
  }
  
  toupper_first = function(a) sub("(.)", "\\U\\1", a, perl=T)
  x = gsub("[_]{2, }", "_", x) %>% strsplit("_") %>% sapply(toupper_first)
  x = x %>% sapply(function(a) paste0(a[a!=""], collapse="_"))
  return(x)
}

unit_into_uM = function(ChEMBL) {
  
  values = ChEMBL$Standard_Value
  units = ChEMBL$Standard_Units
  mw = ChEMBL$Molecular_Weight
  
  values = ifelse(units=="nM", values/1000, values)
  values = ifelse(units=="ug.mL-1", values*1000/mw, values)
  units = ifelse(units %in% c("nM", "ug.mL-1"), "µM", units)
  
  ChEMBL$Standard_Value = values
  ChEMBL$Standard_Units = units
  info_units = ChEMBL$Standard_Units %>% table
  
  for (i in 1:length(info_units)) {
    sprintf("# ChEMBL Units : %s [%s]", names(info_units[i]), info_units[i]) %>% print
  }
  
  return(ChEMBL)
}

read_chembl = function(dir, Anno_Cells=NULL, ln_ic50=T) {
  
  ChEMBL = data.frame()
  files = list.files(dir, full.names=T)
  cells = files %>% strsplit("/") %>% sapply(function(x) x[length(x)])
  cells = gsub(".csv", "", cells)
  
  for (i in 1:length(cells)) {
    Temp = read.csv(files[i], sep=";", na.strings="")
    Temp$SANGER_MODEL_ID = cells[i]
    ChEMBL = ChEMBL %>% rbind(Temp)
  }
  
  colnames(ChEMBL) = colnames(ChEMBL) %>% standard_cols
  col = c("Cell_ChEMBL_ID", "SANGER_MODEL_ID", "Assay_Cell_Type", 
          "Molecule_ChEMBL_ID", "Molecule_Name", 
          "Standard_Type", "Standard_Relation", "Standard_Value", 
          "Standard_Units", "Smiles")
  
  col_ = col[!(col %in% "Molecule_Name")]
  idx = complete.cases(ChEMBL[, col_])
  ChEMBL = ChEMBL[idx, ] %>% 
    mutate(Standard_Relation=gsub("\'", "", Standard_Relation)) %>% 
    relocate(all_of(col), .before=everything()) %>% as.data.frame
  
  # Units of IC50 values are set into uM
  ChEMBL = ChEMBL %>% subset(Standard_Type=="IC50")
  ChEMBL = ChEMBL %>% unit_into_uM
  
  if (ln_ic50) {
    # Log-normalize IC50 values
    ChEMBL = ChEMBL %>% mutate(LN_IC50=log(Standard_Value)) %>% 
      relocate(LN_IC50, .after=Standard_Units) %>% as.data.frame
  }
  
  if (!is.null(Anno_Cells)) {
    idx = match(ChEMBL$SANGER_MODEL_ID, Anno_Cells$SANGER_MODEL_ID)
    ChEMBL = ChEMBL %>% 
      mutate(BROAD_ID=Anno_Cells$BROAD_ID[idx], 
             COSMIC_ID=Anno_Cells$COSMIC_ID[idx]) %>% 
      relocate(BROAD_ID, COSMIC_ID, .after=SANGER_MODEL_ID) %>% as.data.frame
  }
  
  return(ChEMBL)
}

dir = "../../raw_data/ChEMBL/COAD"
ChEMBL_COAD = read_chembl(dir, Anno_Cells=Anno_Cells)   # 4362

# dir = "../../raw_data/ChEMBL/HCC"
# ChEMBL_HCC = read_chembl(dir)   # 187
# 
# dir = "../../raw_data/ChEMBL/LUAD"
# ChEMBL_LUAD = read_chembl(dir)   # 317
# 
# dir = "../../raw_data/ChEMBL/LUSC"
# ChEMBL_LUSC = read_chembl(dir)   # 553
# 
# dir = "../../raw_data/ChEMBL/MEL"
# ChEMBL_MEL = read_chembl(dir)   # 946
# 
# dir = "../../raw_data/ChEMBL/PCM"
# ChEMBL_PCM = read_chembl(dir)   # 98
# 
# dir = "../../raw_data/ChEMBL/SCLC"
# ChEMBL_SCLC = read_chembl(dir)   # 289
# 
# dir = "../../raw_data/ChEMBL/UCEC"
# ChEMBL_UCEC = read_chembl(dir)   # 88


ChEMBL_COAD$SANGER_MODEL_ID %>% table %>% sort(decreasing=T)
# SIDM00840  SIDM00071  SIDM00823  SIDM00500  SIDM00780 
# 4029       287        35         7          4

dir = mkdir("../../processed_data/ic50_data/ChEMBL")
col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", 
        "Molecule_ChEMBL_ID", "Standard_Relation", "LN_IC50")

file = sprintf("%s/IC50_ChEMBL_COAD.csv", dir)
write.csv(ChEMBL_COAD, file=file, row.names=F, fileEncoding="UTF-8")

file = sprintf("%s/IC50_ChEMBL_COAD.txt", dir)
write.table(ChEMBL_COAD[, col], file=file, sep="\t", row.names=F)



# Get PubChem CIDs from the drugs from Inchi-Key
# Remember deep learning models are week in Drug-Blind and Strict-Blind Test
# https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

dir = "../../raw_data/Drug_List"
file = sprintf("%s/ChEMBL_2023_08_25.csv", dir)
Drug_ChEMBL = fread(file, na.strings="")
colnames(Drug_ChEMBL) = colnames(Drug_ChEMBL) %>% standard_cols

drug_coad = unique(ChEMBL_COAD$Molecule_ChEMBL_ID)   # 4007
Drug_COAD = Drug_ChEMBL %>% subset(ChEMBL_ID %in% drug_coad)
drug_coad[!(drug_coad %in% Drug_ChEMBL$ChEMBL_ID)]   # 0

dir = mkdir("../../processed_data/drug_data/ChEMBL")

file = sprintf("%s/SMILES_ChEMBL_COAD.csv", dir)
write.csv(Drug_COAD, file=file, row.names=F)


# Open the results from PubChem Exchange
# https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
file = sprintf("%s/CID_ChEMBL_COAD1.txt", dir)
CID_COAD1 = read.csv(file, sep="\t", header=F)
colnames(CID_COAD1) = c("Inchi_Key", "Drug_CID")

file = sprintf("%s/CID_ChEMBL_COAD2.txt", dir)
CID_COAD2 = read.csv(file, sep="\t", header=F)
colnames(CID_COAD2) = c("SMILES", "Drug_CID")

idx1 = match(Drug_COAD$Inchi_Key, CID_COAD1$Inchi_Key)
idx2 = match(Drug_COAD$Smiles, CID_COAD2$SMILES)

Drug_COAD = Drug_COAD %>% 
  mutate(Drug_CID=CID_COAD1$Drug_CID[idx1]) %>% 
  mutate(Drug_CID=ifelse(!is.na(Drug_CID), Drug_CID, CID_COAD2$Drug_CID[idx2])) %>% 
  relocate(Drug_CID, .after=ChEMBL_ID) %>% as.data.frame
Drug_COAD$Drug_CID %>% is.na %>% sum   # 86


dir = mkdir("../../processed_data/drug_data/ChEMBL")

file = sprintf("%s/SMILES_ChEMBL_COAD.csv", dir)
write.csv(Drug_COAD, file=file, row.names=F)

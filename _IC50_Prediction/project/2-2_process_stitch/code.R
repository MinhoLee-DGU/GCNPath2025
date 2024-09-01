#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(reshape2))
source("../functions.R")
loadings()

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/Anno_Drugs.csv", dir)
Anno_Drugs = read.csv(file, header=T)
Anno_Drugs = Anno_Drugs %>% subset(!is.na(Drug_CID))

file = sprintf("%s/SMILES_GDSC.csv", dir)
SMILES_GDSC = read.csv(file)



##### 2. Process STRING & STiTCH Database

dir = "../../processed_data/net_data/STRING"
file = sprintf("%s/Anno_Genes.csv", dir)
STRING_Target = fread(file)   # 19566
STRING_Target = STRING_Target[, 1:2]
STRING_Target = STRING_Target %>% mutate(ENSEMBL_ID = gsub("^9606.", "", ENSEMBL_ID))

dir = "../../raw_data/STiTCH"
file1 = sprintf("%s/chemicals.v5.0.tsv", dir)
file2 = sprintf("%s/9606.protein_chemical.links.detailed.v5.0.tsv", dir)

STiTCH = fread(file2, sep="\t")                   # 15473939
STiTCH = STiTCH %>% subset(combined_score>=700)   # 466669
col = colnames(STiTCH) %>% standard_cols
colnames(STiTCH) = col
stitch_cid = STiTCH$chemical %>% unique           # 156510


# DON'T RUN THIS CODE WHEN YOU ARE NOT USING A SERVER [~25.6GB...]
STiTCH_Drug = fread(file1, sep="\t")
STiTCH_Drug = STiTCH_Drug %>% subset(chemical %in% stitch_cid)   # 156488
gc()

col = colnames(STiTCH_Drug) %>% standard_cols
colnames(STiTCH_Drug) = col

dir = mkdir("../../processed_data/drug_data/STiTCH")
file1 = sprintf("%s/STiTCH_DTI.csv", dir)
file2 = sprintf("%s/STiTCH_Drug.csv", dir)

fwrite(STiTCH, file=file1, row.names=F)
fwrite(STiTCH_Drug, file=file2, row.names=F)



##### 3-1. Summarize STRING & STiTCH dataset

# Some STRING targets have ensembl gene names
STRING_Target$SYMBOL %>% is.na %>% sum         # 0
grepl("^ENSP", STRING_Target$SYMBOL) %>% sum   # 0
grepl("^ENSG", STRING_Target$SYMBOL) %>% sum   # 393

suppressMessages(library(AnnotationDbi))
suppressMessages(library(EnsDb.Hsapiens.v86))

# keytypes(EnsDb.Hsapiens.v86)
sym_ensg = grep("^ENSG", STRING_Target$SYMBOL, value=T) %>% unique   # 390
Anno_Sym_ENSG = select(EnsDb.Hsapiens.v86, keys=sym_ensg, keytype="GENEID", columns="SYMBOL")

idx = match(STRING_Target$SYMBOL, Anno_Sym_ENSG$GENEID)
STRING_Target$SYMBOL = ifelse(is.na(idx), STRING_Target$SYMBOL, Anno_Sym_ENSG$SYMBOL[idx])

cond = grepl("^ENSG", STRING_Target$SYMBOL)       # 22
STRING_Target = STRING_Target %>% subset(!cond)   # 19566 > 19544


# Match Ensembl ~ Symbol of target in STiTCH DTI
STiTCH = STiTCH %>% mutate(Protein = gsub("^9606.", "", Protein))
idx = match(STiTCH$Protein, STRING_Target$ENSEMBL_ID)
STiTCH = STiTCH %>% mutate(Symbol = STRING_Target$SYMBOL[idx]) %>% 
  relocate(Symbol, .after=Protein)
STiTCH$Symbol %>% is.na %>% sum   # 84226

# keytypes(EnsDb.Hsapiens.v86)
sym_na = STiTCH$Protein[is.na(STiTCH$Symbol)] %>% unique   # 2171
Anno_Sym = select(EnsDb.Hsapiens.v86, keys=sym_na, keytype="PROTEINID", columns="SYMBOL")   # 1728

Anno_Sym$SYMBOL %>% is.na %>% sum   # 0
idx = match(STiTCH$Protein, Anno_Sym$PROTEINID)
STiTCH$Symbol = ifelse(!is.na(STiTCH$Symbol), STiTCH$Symbol, Anno_Sym$SYMBOL[idx])

STiTCH$Symbol %>% is.na %>% sum              # 84226 > 14114
STiTCH = STiTCH %>% subset(!is.na(Symbol))   # 466669 > 452555 [-14114]

colnames(Anno_Sym) = colnames(STRING_Target)
STRING_Target = STRING_Target %>% rbind(Anno_Sym)   # 21272



##### 3-2. Mining of GDSC Target from STiTCH

smiles_can_rcdk = function(smiles, cores=20) {
  
  suppressMessages(library(rcdk))
  
  if (cores<2) {
    smiles = smiles %>% parse.smiles
    smiles_can = lapply(smiles, function(x) get.smiles(x, smiles.flavors("Canonical")))
  } else {
    # suppressMessages(library(doFuture))
    # plan(multisession, workers=cores)
    cluster = clust_future(cores=cores)
    on.exit(stopCluster(cluster))
    smiles_can = foreach(i=1:length(smiles), .combine=c) %dofuture% {
      smiles_ = smiles[i] %>% parse.smiles
      smiles_can_ = smiles_[[1]] %>% get.smiles(smiles.flavors("Canonical"))
      return(smiles_can_)
    }
  }
  return(smiles_can)
}

smiles_can_gdsc = SMILES_GDSC$SMILES_CAN %>% smiles_can_rcdk
smiles_can_stitch = STiTCH_Drug$SMILES_String %>% smiles_can_rcdk
identical(length(smiles_can_gdsc), nrow(SMILES_GDSC))     # T
identical(length(smiles_can_stitch), nrow(STiTCH_Drug))   # T
sum(smiles_can_gdsc %in% smiles_can_stitch)   # 281

SMILES_GDSC$SMILES_CAN_RCDK = smiles_can_gdsc
idx = match(Anno_Drugs$Drug_CID, SMILES_GDSC$Drug_CID)
Anno_Drugs$SMILES_CAN_RCDK = SMILES_GDSC$SMILES_CAN_RCDK[idx]
STiTCH_Drug$SMILES_CAN_RCDK = smiles_can_stitch

# Cf. Differences between CIDs and CIDm in PubChem CID? [X]
STiTCH_Drug$Drug_CID = gsub("CIDs|CIDm", "", STiTCH_Drug$Chemical) %>% as.numeric   # 156488
STiTCH_Drug_Dup = STiTCH_Drug %>% filter(n_distinct(Drug_CID)>1)   # 156488

STiTCH_Drug_Dup %>% group_by(Drug_CID) %>% 
  summarise(N_SMILES=n_distinct(SMILES_String)) %>% 
  pull(N_SMILES) %>% table   # 1 [92432]

rm(STiTCH_Drug_Dup)
STiTCH_Drug = STiTCH_Drug %>% relocate(Drug_CID, .after=Chemical)

# Match GDSC-STiTCH by Name, Synonyms, CID and Canonical SMILES
col = c("Name", "Synonyms", "Drug_CID", "SMILES_CAN_RCDK")
GDSC_Drug_Temp = Anno_Drugs[, col] %>% 
  separate_rows(Name, Synonyms, sep=", |,") %>% 
  mutate(Name=trimws(Name), Synonyms=trimws(Synonyms)) %>% as.data.frame

sum(STiTCH_Drug$Drug_CID %in% GDSC_Drug_Temp$Drug_CID)                 # 469
sum(tolower(STiTCH_Drug$Name) %in% tolower(GDSC_Drug_Temp$Name))       # 417
sum(tolower(STiTCH_Drug$Name) %in% tolower(GDSC_Drug_Temp$Synonyms))   # 100
sum(STiTCH_Drug$SMILES_CAN_RCDK %in% GDSC_Drug_Temp$SMILES_CAN_RCDK)   # 629

STiTCH_Drug$Name %>% is.na %>% sum          # 0
GDSC_Drug_Temp$Name %>% is.na %>% sum       # 0
GDSC_Drug_Temp$Synonyms %>% is.na %>% sum   # 154

idx1 = match(tolower(GDSC_Drug_Temp$Name), tolower(STiTCH_Drug$Name), incomparables=NA)
idx2 = match(tolower(GDSC_Drug_Temp$Synonyms), tolower(STiTCH_Drug$Name), incomparables=NA)
idx3 = match(tolower(GDSC_Drug_Temp$Drug_CID), tolower(STiTCH_Drug$Drug_CID), incomparables=NA)
idx4 = match(tolower(GDSC_Drug_Temp$SMILES_CAN_RCDK), tolower(STiTCH_Drug$SMILES_CAN_RCDK), incomparables=NA)

cid_stitch1 = STiTCH_Drug$Drug_CID[idx1]
cid_stitch2 = STiTCH_Drug$Drug_CID[idx2]
cid_stitch3 = STiTCH_Drug$Drug_CID[idx3]
cid_stitch4 = STiTCH_Drug$Drug_CID[idx4]

GDSC_Drug_Temp = GDSC_Drug_Temp %>% 
  mutate(Drug_CID_STiTCH=cid_stitch1) %>% 
  mutate(Drug_CID_STiTCH=ifelse(!is.na(Drug_CID_STiTCH), Drug_CID_STiTCH, cid_stitch2)) %>% 
  mutate(Drug_CID_STiTCH=ifelse(!is.na(Drug_CID_STiTCH), Drug_CID_STiTCH, cid_stitch3)) %>% 
  mutate(Drug_CID_STiTCH=ifelse(!is.na(Drug_CID_STiTCH), Drug_CID_STiTCH, cid_stitch4))

GDSC_Drug_Temp$Drug_CID_STiTCH %>% is.na %>% sum   # 247
STiTCH$Drug_CID = gsub("CIDs|CIDm", "", STiTCH$Chemical) %>% as.numeric
STiTCH = STiTCH %>% relocate(Drug_CID, .after=Chemical)

cid_gdsc_stitch = GDSC_Drug_Temp$Drug_CID_STiTCH %>% na.omit %>% unique   # 307
GDSC_DTI_STiTCH = STiTCH %>% subset(Drug_CID %in% cid_gdsc_stitch) %>% as.data.frame   # 9986

idx = match(GDSC_DTI_STiTCH$Drug_CID, GDSC_Drug_Temp$Drug_CID_STiTCH)
GDSC_DTI_STiTCH = GDSC_DTI_STiTCH %>% 
  mutate(Drug_Name=GDSC_Drug_Temp$Name[idx], 
         Drug_CID=GDSC_Drug_Temp$Drug_CID[idx])

GDSC_DTI_STiTCH = GDSC_DTI_STiTCH %>% 
  subset(!is.na(Drug_CID)) %>% 
  rename(Target_Symbol=Symbol, Target_Ensembl=Protein) %>% 
  subset(select=-c(Chemical, Experimental, Prediction, Database, Textmining)) %>% 
  relocate(Drug_Name, Drug_CID, Target_Symbol, Target_Ensembl, Combined_Score) %>% 
  distinct %>% as.data.frame   # 6989
all(GDSC_DTI_STiTCH$Drug_CID %in% GDSC_Drug_Temp$Drug_CID)   # T


col = c("Drug_CID", "Target_Symbol")
GDSC_DTI = GDSC_DTI_STiTCH[, col] %>% distinct   # 9986 > 6777
colnames(GDSC_DTI) = c("Drug_CID", "Target")

any_ = function(x) ifelse(length(x)>=1, 1, 0)
GDSC_DTI_Wide = GDSC_DTI %>% reshape2::acast(Drug_CID~Target, value.var="Target", fun.aggregate=any_)   # 302 x 2188
sum(GDSC_DTI_Wide==1)   # 6777
sum(GDSC_DTI_Wide>1)    # 0


# Save files
dir = "../../processed_data/drug_data/STiTCH"

file = sprintf("%s/STiTCH.csv", dir)
fwrite(STiTCH, file=file, row.names=F)

file = sprintf("%s/STRING_Target.csv", dir)
write.csv(STRING_Target, file=file, row.names=F)

file = sprintf("%s/GDSC_DTI.csv", dir)
write.csv(GDSC_DTI_STiTCH, file=file, row.names=F)

file = sprintf("%s/GDSC_DTI_Wide.csv", dir)
write.csv(GDSC_DTI_Wide, file=file, row.names=T)

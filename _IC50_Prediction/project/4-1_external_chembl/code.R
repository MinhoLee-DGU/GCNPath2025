#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
source("../functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

dir = "../../processed_data/ic50_data/GDSC"
file = sprintf("%s/IC50_GDSC.csv", dir)
IC50_GDSC = fread(file)

col = c("SANGER_MODEL_ID", "MODEL_NAME", "SYNONYMS", "BROAD_ID", "COSMIC_ID", "RRID")
Anno_Cells_SANGER = Anno_Cells[, col] %>% 
  separate_rows(SYNONYMS, sep=";") %>% as.data.frame %>% 
  reshape2::melt(measure.vars=c("MODEL_NAME", "SYNONYMS")) %>% 
  distinct %>% as.data.frame   # 4348

colnames(Anno_Cells_SANGER)[5:6] = c("Name_Category", "Cell_Line_Name")



##### 2-1. Process cell annotation
# Run first "sbatch process_assay.sh"

dir = "../../processed_data/ic50_data/ChEMBL"
file = sprintf("%s/Anno_Cells_ChEMBL_Temp.csv", dir)
Anno_Cells_ChEMBL = read.csv(file, na.strings="")   # 1483 x 11

specific = c("ChEMBL", "LINCS", "CLO", "EFO", "CL", "ID")
colnames(Anno_Cells_ChEMBL) = colnames(Anno_Cells_ChEMBL) %>% standard_cols(specific=specific)

Anno_Cells_ChEMBL$Cell_Source_Organism %>% unique   # Homo sapiens
Anno_Cells_ChEMBL$Cell_Name %>% unique %>% length   # 1483

# Match cells between GDSC-ChEMBL by cell-line names and Cellosaurus RRID
cells_sanger = Anno_Cells_SANGER$Cell_Line_Name %>% toupper %>% gsub("-", "", .)
cells_chembl = Anno_Cells_ChEMBL$Cell_Name %>% toupper %>% gsub("-", "", .)

cells_sanger %>% is.na %>% sum   # 1614
cells_chembl %>% is.na %>% sum   # 0
Anno_Cells_SANGER$RRID %>% is.na %>% sum   # 1028
Anno_Cells_ChEMBL$Cellosaurus_ID %>% is.na %>% sum   # 148

idx1 = match(cells_chembl, cells_sanger, incomparables=NA)
idx2 = match(Anno_Cells_ChEMBL$Cellosaurus_ID, Anno_Cells_SANGER$RRID, incomparables=NA)
idx = coalesce(idx1, idx2)

# Cf. SIDM00400 has multiple RRIDs
# SIDM00400, CVCL_2717;CVCL_1888, ACH-002392;ACH-002303
Anno_Cells_SANGER %>% subset(grepl(";", RRID))

# # Already matched by cell-line name
# which(Anno_Cells_ChEMBL$Cellosaurus_ID=="CVCL_2717")   # NA
# which(Anno_Cells_ChEMBL$Cellosaurus_ID=="CVCL_1888")   # 1398
# Anno_Cells_ChEMBL$Cell_Name[1398]             # SC1
# Anno_Cells_SANGER$Cell_Line_Name[idx[1398]]   # SC-1

Anno_Cells_ChEMBL = Anno_Cells_ChEMBL %>% cbind(Anno_Cells_SANGER[idx, ])
rownames(Anno_Cells_ChEMBL) = NULL

# Filter out cell-lines not included in GDSC
Anno_Cells_ChEMBL$Cell_GDSC = Anno_Cells_ChEMBL$SANGER_MODEL_ID %in% unique(IC50_GDSC$SANGER_MODEL_ID)
Anno_Cells_ChEMBL$Cell_GDSC %>% sum   # 765

# Filter out cell-lines without names and ids
col = c("SANGER_MODEL_ID", "Cell_Line_Name", "BROAD_ID", "COSMIC_ID")
Anno_Cells_ChEMBL[, col] %>% apply(2, function(x) sum(!is.na(x)))
# SANGER_MODEL_ID  Cell_Line_Name  BROAD_ID  COSMIC_ID 
# 915              915             891       828

idx = Anno_Cells_ChEMBL[, col] %>% is.na %>% apply(1, all)
Anno_Cells_ChEMBL = Anno_Cells_ChEMBL[!idx, ]   # 1483 > 915

# Fortunately, all cells have SANGER IDs...
Anno_Cells_ChEMBL$SANGER_MODEL_ID %>% is.na %>% sum   # 0

dir = mkdir("../../processed_data/ic50_data/ChEMBL")
file = sprintf("%s/Anno_Cells_ChEMBL.csv", dir)
write.csv(Anno_Cells_ChEMBL, file=file, row.names=F)



##### 2-2. Process assay annotation

dir = "../../processed_data/ic50_data/ChEMBL"
file = sprintf("%s/Assay_ChEMBL.csv", dir)
Assay_ChEMBL = fread(file, na.strings="")   # 280888 x 28

specific = c("ChEMBL", "BAO", "ID")
names(Assay_ChEMBL) = names(Assay_ChEMBL) %>% standard_cols(specific=specific)
all(Assay_ChEMBL$Cell_ChEMBL_ID %in% Anno_Cells_ChEMBL$Cell_ChEMBL_ID)   # T

# Exclude Drug-Target Interaction data
# https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
Assay_ChEMBL = Assay_ChEMBL[Confidence_Score==1]           # 280888 > 242053
Assay_ChEMBL = Assay_ChEMBL[Assay_Type %in% c("F", "T")]   # 242053 > 230282

# Exclude Non-Human data
Assay_ChEMBL = Assay_ChEMBL[Assay_Organism=="Homo sapiens"]   # 230282 > 221909

# Select Cell-based Format
Assay_ChEMBL = Assay_ChEMBL[BAO_Label=="cell-based format"]   # 221909 > 213881

# Exclude high-throughput response data [SANGER, GDSC, CCLE, CTRP, NCI-60, PRISM]
# https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/document-and-data-source-questions
Assay_ChEMBL = Assay_ChEMBL[Src_ID!=5]   # 213881 > 213185
grep("SANGER|GDSC|CCLE|CTRP|NCI-60|PRISM|DepMap", Assay_ChEMBL$Description, value=T, ignore.case=T)   # 57 [NCI-60]
grep("Cancer Cell Line Encyclopedia", Assay_ChEMBL$Description, value=T, ignore.case=T)               # 0 [CCLE]
grep("Cancer Therapeutics Response Portal", Assay_ChEMBL$Description, value=T, ignore.case=T)         # 0 [CTRP]
grep("Genomics of Drug Sensitivity in Cancer", Assay_ChEMBL$Description, value=T, ignore.case=T)      # 0 [GDSC]
Assay_ChEMBL = Assay_ChEMBL %>% subset(!grepl("NCI-60", Description))   # 213185 > 213128

dir  = "../../processed_data/ic50_data/ChEMBL"
file = sprintf("%s/Assay_ChEMBL.csv", dir)
fwrite(Assay_ChEMBL, file=file)



##### 3. Process IC50 Annotation
# Run first "sbatch process_ic50.sh"

pattern = "IC50_ChEMBL_Temp_[0-9]+.csv"
dir = "../../processed_data/ic50_data/ChEMBL/Temp"
file_list = list.files(dir, pattern=pattern, full.names=T)

IC50_ChEMBL_List = list()
for (i in 1:length(file_list)) {
  IC50_ChEMBL_List[[i]] = fread(file_list[i], na.strings="")
}

IC50_ChEMBL = rbindlist(IC50_ChEMBL_List)   # 434542 x 46
specific = c("ChEMBL", "pChEMBL", "SMILES", "QUDT", "TOID", "BAO", "ID", "UO")
names(IC50_ChEMBL) = names(IC50_ChEMBL) %>% standard_cols(specific=specific)

all(IC50_ChEMBL$Standard_Type=="IC50")   # T
all(IC50_ChEMBL$Assay_ChEMBL_ID %in% Assay_ChEMBL$Assay_ChEMBL_ID)   # T
idx = match(IC50_ChEMBL$Assay_ChEMBL_ID, Assay_ChEMBL$Assay_ChEMBL_ID)

IC50_ChEMBL = IC50_ChEMBL %>% 
  mutate(Cell_ChEMBL_ID = Assay_ChEMBL$Cell_ChEMBL_ID[idx]) %>% 
  relocate(Cell_ChEMBL_ID, .before=Molecule_ChEMBL_ID)

all(IC50_ChEMBL$Cell_ChEMBL_ID %in% Anno_Cells_ChEMBL$Cell_ChEMBL_ID)   # T
col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "RRID", "Cell_Line_Name", "Cell_GDSC")
idx = match(IC50_ChEMBL$Cell_ChEMBL_ID, Anno_Cells_ChEMBL$Cell_ChEMBL_ID)

IC50_ChEMBL = IC50_ChEMBL %>% cbind(Anno_Cells_ChEMBL[idx, col]) %>% 
  relocate(SANGER_MODEL_ID, BROAD_ID, COSMIC_ID, 
           RRID, Cell_Line_Name, Cell_GDSC, .after=Cell_ChEMBL_ID)


### 3-1. Filter out NA in IC50 values 
# [Standard_Value, Standard_Units, Standard_Relation]

IC50_ChEMBL$Standard_Value %>% is.na %>% sum      # 0
IC50_ChEMBL$Standard_Units %>% is.na %>% sum      # 0
IC50_ChEMBL$Standard_Relation %>% is.na %>% sum   # 0

### 3-2. Filter out units of IC50 values with trivial sample numbers
# The main unit is nM, consistent with the paper...
# [2014] The ChEMBL bioactivity database: an update
# [2014] 10.1093/nar/gkt1031

IC50_ChEMBL$Standard_Units %>% table %>% sort(decreasing=T)
# nM      ug.mL-1  ug  (The Rest...)
# 416761  17614    88  

IC50_ChEMBL = IC50_ChEMBL %>% 
  subset(Standard_Units %in% c("nM", "ug.mL-1"))   # 434542 > 434375


### 3-3. Filter out NA in drugs

IC50_ChEMBL$Canonical_SMILES %>% is.na %>% sum     # 1499
IC50_ChEMBL$Molecule_ChEMBL_ID %>% is.na %>% sum   # 0
IC50_ChEMBL = IC50_ChEMBL %>% subset(!is.na(Canonical_SMILES))   # 434375 > 432876
drug_chembl = IC50_ChEMBL$Molecule_ChEMBL_ID %>% unique   # 148455

dir = mkdir("../../processed_data/drug_data/ChEMBL")
file = sprintf("%s/Drug_ChEMBL.txt", dir)
cat(drug_chembl, file=file, sep="\n")

# Find CIDs from PubChem
# https://pubchem.ncbi.nlm.nih.gov/idexchange

file = sprintf("%s/CIDs_ChEMBL.txt", dir)
CIDs_ChEMBL = read.csv(file, sep="\t", header=F)   # 148455 x 2
colnames(CIDs_ChEMBL) = c("Drug_ChEMBL_ID", "Drug_CID")

CIDs_ChEMBL$Drug_CID %>% is.na %>% sum   # 8
idx = match(IC50_ChEMBL$Molecule_ChEMBL_ID, CIDs_ChEMBL$Drug_ChEMBL_ID)

IC50_ChEMBL = IC50_ChEMBL %>% 
  mutate(Drug_CID = CIDs_ChEMBL$Drug_CID[idx]) %>% 
  mutate(Drug_GDSC = Drug_CID %in% unique(IC50_GDSC$DRUG_CID)) %>% 
  relocate(Drug_CID, Drug_GDSC, .after=Molecule_ChEMBL_ID)



### 3-4. Transform IC50 values to ln(uM) scales like GDSC...

IC50_ChEMBL = IC50_ChEMBL %>% 
  mutate(LN_IC50 = log(Standard_Value/1000)) %>% 
  relocate(LN_IC50, .after=Standard_Value)

# Filter out infinite or minus values in ln(IC50)
IC50_ChEMBL$LN_IC50 %>% range %>% round(3)               # [NaN, NaN]
IC50_ChEMBL$LN_IC50 %>% na.omit %>% range %>% round(3)   # [-Inf, 14.286]

IC50_ChEMBL %>% subset(is.na(LN_IC50)) %>% pull(Standard_Value) %>% range   # [-100000, -563]
IC50_ChEMBL = IC50_ChEMBL %>% subset(!is.na(LN_IC50))   # 432876 > 432859

sum(IC50_ChEMBL$Standard_Value==0)   # 44
idx = (IC50_ChEMBL$Standard_Relation=="=" & IC50_ChEMBL$Standard_Value==0)   # 44
IC50_ChEMBL = IC50_ChEMBL[!idx, ]    # 432859 > 432815
sum(IC50_ChEMBL$Standard_Value==0)   # 0

IC50_ChEMBL$LN_IC50 %>% range %>% round(3)   # [-25.328, 14.286]
IC50_ChEMBL$LN_IC50 %>% hist


### 3-5. Average ln(IC50) of duplicated cell x drug combination

average_dup = function(IC50) {
  
  n_before = IC50 %>% nrow
  num_tags = IC50$TAG %>% table
  dup_tags = names(num_tags)[num_tags>1]
  
  col = c("Cell_ChEMBL_ID", "Molecule_ChEMBL_ID", "LN_IC50", "TAG")
  IC50_Dup_Raw = IC50[, col] %>% subset(TAG %in% dup_tags)
  IC50_Dup = IC50_Dup_Raw %>% group_by(TAG) %>% 
    summarise(LN_IC50_Mean=mean(LN_IC50), N_Dup=n()) %>% as.data.frame
  
  idx = match(IC50_Dup_Raw$TAG, IC50_Dup$TAG)
  IC50_Dup_Raw = IC50_Dup_Raw %>% 
    mutate(LN_IC50_Mean=IC50_Dup$LN_IC50_Mean[idx]) %>% 
    relocate(TAG, .after=everything())
  
  idx = match(IC50$TAG, IC50_Dup$TAG)
  IC50$LN_IC50 = ifelse(is.na(idx), IC50$LN_IC50, IC50_Dup$LN_IC50_Mean[idx])
  IC50 = IC50 %>% subset(!duplicated(TAG))
  
  n_after = IC50 %>% nrow
  sprintf("# Remove Duplicates : %s > %s", n_before, n_after) %>% print
  IC50_List = list(IC50=IC50, IC50_Dup=IC50_Dup, IC50_Dup_Raw=IC50_Dup_Raw)
  return(IC50_List)
}

# Drug Annotation
col = c("Molecule_ChEMBL_ID", "Molecule_Pref_Name", "Drug_CID", "Drug_GDSC", "Canonical_SMILES")
Anno_Drugs_ChEMBL = IC50_ChEMBL[, col, with=F] %>% distinct %>% as.data.frame   # 148415
Anno_Drugs_ChEMBL$Canonical_SMILES %>% unique %>% length     # 148410
Anno_Drugs_ChEMBL$Molecule_ChEMBL_ID %>% unique %>% length   # 148415

# Unify ChEMBL IDs sharing the same SMILES
Anno_Drugs_ChEMBL = Anno_Drugs_ChEMBL %>% subset(!duplicated(Canonical_SMILES))   # 148410
idx = match(IC50_ChEMBL$Canonical_SMILES, Anno_Drugs_ChEMBL$Canonical_SMILES)
IC50_ChEMBL$Molecule_ChEMBL_ID = Anno_Drugs_ChEMBL$Molecule_ChEMBL_ID[idx]

IC50_ChEMBL$Molecule_ChEMBL_ID %>% unique %>% length   # 148410
Anno_Drugs_ChEMBL$Drug_GDSC %>% sum   # 233

IC50_ChEMBL_EQ = IC50_ChEMBL %>% subset(Standard_Relation=="=")   # 341436
IC50_ChEMBL_NE = IC50_ChEMBL %>% subset(Standard_Relation!="=")   # 91379

col = c("Assay_ChEMBL_ID", "Activity_ID", "Cell_ChEMBL_ID", "Molecule_ChEMBL_ID", "LN_IC50")
IC50_ChEMBL_EQ_List = IC50_ChEMBL_EQ[, col] %>% 
  mutate(TAG = paste0(Cell_ChEMBL_ID, "@", Molecule_ChEMBL_ID))
IC50_ChEMBL_EQ_List = IC50_ChEMBL_EQ_List %>% as.data.frame %>% average_dup

rmse = IC50_ChEMBL_EQ_List$IC50_Dup_Raw %>% with(RMSE(LN_IC50, LN_IC50_Mean))   # 1.64929
corr = IC50_ChEMBL_EQ_List$IC50_Dup_Raw %>% with(cor(LN_IC50, LN_IC50_Mean))    # 0.8735085

dir = "../../processed_data/ic50_data/ChEMBL"
main = sprintf("%s/LN_IC50 & LN_IC50_Mean [ChEMBL]", dir)

xlab = bquote(ln(IC[50]))
ylab = bquote(Mean~ln(IC[50]))

IC50_ChEMBL_EQ_List$IC50_Dup_Raw %>% 
  plot_def(LN_IC50, LN_IC50_Mean, main=main, alpha=0.5, xlab=xlab, ylab=ylab, 
           axis_tl=25, axis_tx=20, xy_line=T, force_bold=F, save=T, dpi=1500)

IC50_ChEMBL_EQ_ = IC50_ChEMBL_EQ_List$IC50 %>% subset(select=-c(TAG))
IC50_ChEMBL = rbind(IC50_ChEMBL_EQ_, IC50_ChEMBL_NE)

IC50_ChEMBL$Standard_Relation %>% table
# <     <=   =       >      >=   >>  ~
# 3293  276  302648  86814  974  7   15



##### 3. Save files...

IC50_ChEMBL_ = IC50_ChEMBL %>% 
  subset(select=c(Assay_ChEMBL_ID, Cell_ChEMBL_ID, SANGER_MODEL_ID, 
                  BROAD_ID, COSMIC_ID, RRID, Cell_Line_Name, Cell_GDSC,
                  Molecule_ChEMBL_ID, Molecule_Pref_Name, Drug_CID, Drug_GDSC, 
                  Standard_Relation, Standard_Value, LN_IC50))

col = c("SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "Cell_Line_Name")
IC50_ChEMBL_[, col] %>% sapply(function(x) x %>% na.omit %>% unique %>% length)
# SANGER_MODEL_ID   Cell_Line_Name   BROAD_ID   COSMIC_ID 
# 606               586              540        607 

IC50_ChEMBL_$Standard_Relation %>% table
# <    <=  =      >      >= 
# 243  28  38886  12827  116


# IC50 ChEMBL
dir = "../../processed_data/drug_data/ChEMBL"
file = sprintf("%s/SMILES_ChEMBL.csv", dir)
fwrite(Anno_Drugs_ChEMBL, file=file)

dir = "../../processed_data/ic50_data/ChEMBL"
file = sprintf("%s/IC50_ChEMBL.csv", dir)
fwrite(IC50_ChEMBL, file=file)

file = sprintf("%s/IC50_ChEMBL.txt", dir)
fwrite(IC50_ChEMBL_, file=file, sep="\t")

file = sprintf("%s/IC50_ChEMBL.xlsx", dir)
write.xlsx(IC50_ChEMBL_EQ_List, file=file, rowNames=F)


# IC50_GDSC_ChEMBL = data.frame(LN_IC50=c(IC50_GDSC$LN_IC50, IC50_ChEMBL_EQ_$LN_IC50), 
#                               Dataset=c(rep("GDSC", nrow(IC50_GDSC)), rep("ChEMBL", nrow(IC50_ChEMBL_EQ_))))
# 
# levels = c("GDSC", "ChEMBL")
# IC50_GDSC$LN_IC50 %>% round(3) %>% range   # -9.800 13.847
# IC50_GDSC_ChEMBL$Dataset = IC50_GDSC_ChEMBL$Dataset %>% factor(levels=levels)
# 
# suppressMessages(library(ggpubr))
# margin = margin(10, 10, 10, 10, unit="pt")
# 
# font_label = font("xylab", size=20)
# font_text = font("xy.text", size=16, color="grey30", margin=margin)
# font_lgt = font("legend.title", size=14.4)
# font_lgx = font("legend.text", size=14.4)
# 
# font = font_label + font_text + font_lgt + font_lgx
# color = scale_color_manual(values=c("firebrick2", "royalblue2"), labels=levels)
# fill = scale_fill_manual(values=c("firebrick1", "royalblue1"), labels=levels)
# 
# pl = IC50_GDSC_ChEMBL %>% 
#   ggdensity(x="LN_IC50", color="Dataset", fill="Dataset") %>% 
#   ggpar(legend="right") + labs(x=xlab) + font + color + fill
# 
# main = sprintf("%s/Histogram of LN_IC50 [GDSC & ChEMBL]", dir)
# pl %>% save_fig_ggpubr(main=main, width=16, height=12, svg=T)


### Cf. Re-Visualization of IC50 Distribution

prepared_tag_ref = T
if (prepared_tag_ref) {
  dir = "../../do_better/_performance_chembl"
  file = sprintf("%s/Cell_Drug_Pair.csv", dir)
  Cell_Drug_Pair = fread(file)   # 237580
  
  IC50_ChEMBL_List_Final = list()
  tag_ref = Cell_Drug_Pair %>% with(paste0(Cell_ChEMBL_ID, "@", Drug_ChEMBL_ID))
  by = c("Cell_ChEMBL_ID"="Cell_ChEMBL_ID", "Molecule_ChEMBL_ID"="Drug_ChEMBL_ID")
  
  IC50_ChEMBL_List_Final$IC50 = right_join(IC50_ChEMBL_EQ_List$IC50, Cell_Drug_Pair, by=by)
  IC50_ChEMBL_List_Final$IC50_Dup = IC50_ChEMBL_EQ_List$IC50_Dup %>% subset(TAG %in% tag_ref)
  IC50_ChEMBL_List_Final$IC50_Dup_Raw = IC50_ChEMBL_EQ_List$IC50_Dup_Raw %>% subset(TAG %in% tag_ref)
  IC50_ChEMBL_Final = IC50_ChEMBL_List_Final$IC50
  
  IC50_ChEMBL_EQ_List %>% sapply(nrow)
  # [IC50] 302648, [IC50_Dup] 17704, [IC50_Dup_Raw] 56492
  IC50_ChEMBL_List_Final %>% sapply(nrow)
  # [IC50] 237580, [IC50_Dup] 12976, [IC50_Dup_Raw] 33213 
  
  
  # 1. Duplicated IC50
  rmse = IC50_ChEMBL_List_Final$IC50_Dup_Raw %>% with(RMSE(LN_IC50, LN_IC50_Mean))   # 1.250103
  corr = IC50_ChEMBL_List_Final$IC50_Dup_Raw %>% with(cor(LN_IC50, LN_IC50_Mean))    # 0.9302273
  
  dir = "../../processed_data/ic50_data/ChEMBL"
  main = sprintf("%s/LN_IC50 & LN_IC50_Mean [ChEMBL, Final]", dir)
  
  xlab = bquote(ln(IC[50]))
  ylab = bquote(Mean~ln(IC[50]))
  
  IC50_ChEMBL_List_Final$IC50_Dup_Raw %>% 
    plot_def(LN_IC50, LN_IC50_Mean, main=main, xlab=xlab, ylab=ylab, 
             size=1.5, alpha=0.5, axis_tl=30, axis_tx=24, dpi=1200,
             xy_line=T, raster=T, save=T, save_svg=T)
  
  
  # 2. ChEMBL vs GDSC
  col = c("Cell", "Drug", "LN_IC50")
  col_gdsc = c("SANGER_MODEL_ID", "DRUG_CID", "LN_IC50")
  col_chembl = c("Cell_ChEMBL_ID", "Molecule_ChEMBL_ID", "LN_IC50")
  dataset = c(rep("GDSC", nrow(IC50_GDSC)), rep("ChEMBL", nrow(IC50_ChEMBL_Final)))
  
  # idx1 = match(IC50_ChEMBL_Final$Cell_ChEMBL_ID, Anno_Cells_ChEMBL$Cell_ChEMBL_ID)
  # idx2 = match(IC50_ChEMBL_Final$Molecule_ChEMBL_ID, Anno_Drugs_ChEMBL$Molecule_ChEMBL_ID)
  # 
  # IC50_ChEMBL_Temp = IC50_ChEMBL_Final %>% 
  #   mutate(Cell=Anno_Cells_ChEMBL$SANGER_MODEL_ID[idx1], 
  #          Drug=Anno_Drugs_ChEMBL$Drug_CID[idx2]) %>% 
  #   subset(select=c(Cell, Drug, LN_IC50)) %>% as.data.frame
  # 
  # IC50_ChEMBL_Temp$Cell %>% is.na %>% sum   # 0
  # IC50_ChEMBL_Temp$Drug %>% is.na %>% sum   # 0
  # 
  # IC50_GDSC_ChEMBL = rbind(
  #   IC50_GDSC[, col_gdsc, with=F] %>% setNames(col) %>% as.data.frame, 
  #   IC50_ChEMBL_Temp
  # ) %>% mutate(Dataset=dataset)
  
  IC50_GDSC_ChEMBL = rbind(
    IC50_GDSC[, col_gdsc, with=F] %>% setNames(col) %>% as.data.frame, 
    IC50_ChEMBL_Final[, col_chembl] %>% setNames(col) %>% as.data.frame
  ) %>% mutate(Dataset=dataset)
  
  levels = c("GDSC", "ChEMBL")
  IC50_GDSC$LN_IC50 %>% round(3) %>% range   # -9.800 13.847
  IC50_GDSC_ChEMBL$Dataset = IC50_GDSC_ChEMBL$Dataset %>% factor(levels=levels)
  IC50_GDSC_ChEMBL %>% group_by(Dataset) %>% summarise(LN_IC50_Median=median(LN_IC50))
  # [GDSC] 2.85, [ChEMBL] 1.46
  
  
  suppressMessages(library(ggpubr))
  margin1 = margin(0.2, 0.2, 0.2, 0.2, unit="cm")
  margin2 = margin(b=0.5, unit="cm")
  margin3 = margin(b=0.1, l=0.25, t=0.1, unit="cm")
  theme_lg = theme(legend.key.width=unit(0.8, "cm"), 
                   legend.key.height=unit(0.6, "cm"))
  
  font_label = font("xylab", size=25, margin=margin1)
  font_text = font("xy.text", size=20, color="grey30", margin=margin1)
  font_lgt = font("legend.title", size=20, margin=margin2)
  font_lgx = font("legend.text", size=20, margin=margin3)
  
  ylab = "Density"
  font = font_label + font_text + font_lgt + font_lgx
  color = scale_color_manual(values=c("firebrick2", "royalblue2"), labels=levels)
  fill = scale_fill_manual(values=c("firebrick1", "royalblue1"), labels=levels)
  
  pl = IC50_GDSC_ChEMBL %>% 
    ggdensity(x="LN_IC50", color="Dataset", fill="Dataset") %>% 
    ggpar(legend="right") + labs(x=xlab, y=ylab) + font + color + fill + theme_lg
  
  main = sprintf("%s/Histogram of LN_IC50 [GDSC & ChEMBL, Final]", dir)
  pl %>% save_fig_ggpubr(main=main, width=18, height=12, svg=T)
}


supplementary = T
if (supplementary) {
  suppressMessages(library(openxlsx))
  
  ### [Supplementary Data] Supplementary Data 7
  Anno_Drugs_ChEMBL_ = Anno_Drugs_ChEMBL %>% 
    subset(Molecule_ChEMBL_ID %in% unique(Cell_Drug_Pair$Drug_ChEMBL_ID)) %>% 
    rename(GDSC_Drug=Drug_GDSC) %>% arrange() %>% as.data.frame   # 112946
  all(Anno_Drugs_ChEMBL_$GDSC_Drug)   # F
  
  sheets = "Supplementary Data 7"
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(Anno_Drugs_ChEMBL_, file=file, sheetName=sheets, rowNames=F)
  
  
  ### [Supplementary Data] Supplementary Data 8
  by = c("Molecule_ChEMBL_ID"="Drug_ChEMBL_ID", "Cell_ChEMBL_ID"="Cell_ChEMBL_ID")
  IC50_ChEMBL_EQ_1 = IC50_ChEMBL_ %>% 
    subset(Standard_Relation=="=") %>% 
    right_join(Cell_Drug_Pair, by=by) %>% 
    rename(GDSC_Cell=Cell_GDSC, GDSC_Drug=Drug_GDSC) %>% as.data.frame   # 237580
  
  sheets = "Supplementary Data 8"
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(IC50_ChEMBL_EQ_1, file=file, sheetName=sheets, rowNames=F)
  
  
  ### [Source Data] Supplementary Fig. 34
  IC50_Dup_ = IC50_ChEMBL_List_Final$IC50_Dup_Raw %>% subset(select=-TAG)
  IC50_GDSC_ChEMBL_ = IC50_GDSC_ChEMBL %>% 
    mutate(Dataset=recode(Dataset, "GDSC"="GDSC1+2"))
  df_list = list(IC50_Dup_, IC50_GDSC_ChEMBL_)
  
  num_fig = letters[1:2]
  sheets = sprintf("Supplementary Fig. 34%s", num_fig)
  file = "SourceData_SupplementaryFig34.xlsx"
  write.xlsx(df_list, file=file, sheetName=sheets, rowNames=F)
}

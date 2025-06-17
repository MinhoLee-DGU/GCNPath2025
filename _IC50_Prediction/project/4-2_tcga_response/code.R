#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(cogena))
suppressMessages(library(TCGAbiolinks))
suppressMessages(library(SummarizedExperiment))

source("../functions.R")
loadings()

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.entrez.gmt", dir)
Path_List = gmt2list(file)
geneset = Path_List %>% unlist %>% unique   # 1509

dir = "../../processed_data/cell_data/SANGER_Passports"

file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

file = sprintf("%s/TPM_Ent.csv", dir)
SANGER_RNA_Ent = fread_def(file, check_names=F)

dir = "../../processed_data/cell_data/BIOCARTA"
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
SANGER_RNA_GSVA = read.csv(file, row.names=1, check.names=F)



##### 2-1. Parse TCGA Data [RNA-Seq]

gsva_def = function(Omics, Path_List, method="gsva", 
                    filt_genes=T, do_qnorm=F, cores=10, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  suppressMessages(library(stats))
  gene_list = Path_List %>% unlist %>% unique
  
  n_before = length(gene_list)
  n_after = sum(colnames(Omics) %in% gene_list)
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    n_before = ncol(Omics)
    Omics = Omics[, colnames(Omics) %in% gene_list]
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  Omics = Omics %>% t %>% as.matrix
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = Omics %>% normalize.quantiles(keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = Omics_Path %>% t %>% as.data.frame
  return(Omics_Path)
}

dir = mkdir("../../raw_data/TCGA")
gdc_code = getGDCprojects()$project_id
tcga_code = grep("TCGA", gdc_code, value=T)   # 33

query_tcga = GDCquery(project=tcga_code, access="open",
                      data.category="Transcriptome Profiling", 
                      data.type="Gene Expression Quantification")

# Download from GDC repository
GDCdownload(query_tcga, directory=dir)
# Make R object from the downloaded data
data = GDCprepare(query_tcga, directory=dir)
# Extract Gene expression matrix
TCGA_RNA = assay(data, "tpm_unstrand")
# Genes x Samples [60660 x 11274]
TCGA_RNA = TCGA_RNA %>% t %>% as.data.frame
# Samples x Genes [11274 x 60660]
gc()

# # TCGA_RNA is truly the TPM data
# tcga_rowsum = TCGA_RNA %>% apply(1, sum)
# tcga_rowsum %>% range   # 1e+06 1e+06

TCGA_Info = query_tcga$results[[1]]
TCGA_Info$project %>% table %>% sort(decreasing=T)
identical(rownames(TCGA_RNA), TCGA_Info$cases)   # T

dir = mkdir("../../processed_data/cell_data/TCGA")
file = sprintf("%s/TCGA_RNA.csv", dir)
fwrite(TCGA_RNA, file=file, row.names=T)

file = sprintf("%s/TCGA_Info.csv", dir)
fwrite(TCGA_Info, file=file)

file = sprintf("%s/TCGA_RNA.RData", dir)
save(TCGA_RNA, TCGA_Info, query_tcga, file=file)


# [TCGA_RNA] Ensembl ID > Entrez ID
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))

keytypes(org.Hs.eg.db)
columns = c("ENTREZID", "SYMBOL")
key = gsub("\\.[0-9]+", "", colnames(TCGA_RNA))
Anno_Genes_TCGA = AnnotationDbi::select(org.Hs.eg.db, key, keytype="ENSEMBL", columns=columns)   # 60895

col = c("ENSEMBL_GENE_ID", "ENTREZ_ID", "HGNC_SYMBOL")
Anno_Genes_TCGA_ = Anno_Genes[, col] %>% subset(!is.na(ENSEMBL_GENE_ID))
colnames(Anno_Genes_TCGA_) = colnames(Anno_Genes_TCGA)

Anno_Genes_TCGA = Anno_Genes_TCGA %>% rbind(Anno_Genes_TCGA_)   # 98347
Anno_Genes_TCGA$ENSEMBL %>% is.na %>% sum   # 0

Anno_Genes_TCGA = Anno_Genes_TCGA %>% 
  group_by(ENSEMBL) %>% tidyr::fill(ENTREZID, SYMBOL, .direction="downup")

Anno_Genes_TCGA$ENTREZID %>% is.na %>% sum   # 19882
Anno_Genes_TCGA = Anno_Genes_TCGA %>% distinct(.keep_all=T)   # 69135
sum(geneset %in% Anno_Genes_TCGA$ENTREZID)   # 1509 [1509]
# All BIOCARTA genes exist in TCGA genes

key %>% unique %>% length
# 60660 [from 60660, all unique]
all(key %in% Anno_Genes_TCGA$ENSEMBL)
# T [All Ensembl IDs match to Entrez IDs]
Anno_Genes_TCGA$ENTREZID %>% is.na %>% sum
# 19877 [Some genes do not have the corresponding Entrez IDs]

# Ensembl ID > Entrez ID
TCGA_RNA_Ent = TCGA_RNA
idx = match(key, Anno_Genes_TCGA$ENSEMBL)
col = Anno_Genes_TCGA$ENTREZID[idx]

# Duplicated genes exist, drop them
col %>% na.omit %>% duplicated %>% sum   # 156
idx_no = is.na(col) | duplicated(col)
TCGA_RNA_Ent = TCGA_RNA_Ent[, !idx_no]   # 60660 > 40629
colnames(TCGA_RNA_Ent) = col[!idx_no]

# Log2 transformation with pseudo count 1
TCGA_RNA_Ent %>% rowSums %>% range   # 866270.5 998629.4
TCGA_RNA_Ent = log2(TCGA_RNA_Ent+1) %>% as.data.frame

# Calculate GSVA pathway scores
cores = 48
sum(colnames(TCGA_RNA_Ent) %in% geneset)   # 1505/1509 [99.73%]
TCGA_RNA_GSVA = TCGA_RNA_Ent %>% gsva_def(Path_List, method="gsva", cores=cores)


# Ensembl ID > Symbol
TCGA_RNA_Sym = TCGA_RNA
idx = match(key, Anno_Genes_TCGA$ENSEMBL)
col = Anno_Genes_TCGA$SYMBOL[idx]

# Duplicated genes exist, drop them
col %>% na.omit %>% duplicated %>% sum   # 156
idx_no = is.na(col) | duplicated(col)
TCGA_RNA_Sym = TCGA_RNA_Sym[, !idx_no]   # 60660 > 40634
colnames(TCGA_RNA_Sym) = col[!idx_no]

# Log2 transformation with pseudo count 1
TCGA_RNA_Sym %>% rowSums %>% range   # 866282.4 998629.4
TCGA_RNA_Sym = log2(TCGA_RNA_Sym+1) %>% as.data.frame


# ### Cf. How many tumor-normal samples are from the same tissues?
# # Those samples will be also utilize to compare predictions [Tumor vs Normal]
# col = c("cases.submitter_id", "cases", "sample_type", "project")
# TCGA_Info_TN = TCGA_Info[, col] %>% group_by(cases.submitter_id) %>% 
#   filter(n_distinct(cases)>1) %>% as.data.frame %>% 
#   mutate(tumor_normal=ifelse(sample_type=="Solid Tissue Normal", "Normal", "Tumor"))   # 1714
# 
# TCGA_Info_TN = TCGA_Info_TN %>% group_by(cases.submitter_id) %>% 
#   filter(n_distinct(tumor_normal)>1 & n_distinct(project)==1) %>% as.data.frame   # 1714 > 1462
# 
# TCGA_Info_TN$project %>% unique %>% length              # 23
# TCGA_Info_TN$cases.submitter_id %>% unique %>% length   # 716
# TCGA_Info_TN = TCGA_Info_TN %>% relocate(project, .after=everything()) %>% as.data.frame
# colnames(TCGA_Info_TN) = c("Patient", "Sample", "Sample_Type", "Tumor_Normal", "Project")


dir = mkdir("../../processed_data/cell_data/TCGA")

# file = sprintf("%s/TCGA_Info_TN.csv", dir)
# fwrite(TCGA_Info_TN, file=file)

file = sprintf("%s/TPM_Ent.csv", dir)
fwrite(TCGA_RNA_Ent, row.names=T, file=file)

file = sprintf("%s/Anno_Genes.csv", dir)
fwrite(Anno_Genes_TCGA, file=file)

dir = "../../processed_data/cell_data/BIOCARTA"
file = sprintf("%s/TCGA_RNA_GSVA.csv", dir)
fwrite(TCGA_RNA_GSVA, row.names=T, file=file)



##### 2-1. Parse TCGA Data [RNA-Seq, PCA Analysis]
# Does GSVA virtually do batch correction?

# # Cf. Is the result of parallelPCA equal to prcomp?
# suppressMessages(library(PCAtools))
# data(iris)
# pca_x1 = iris %>% prcomp(center=T, scale=T)
# pca_x2 = iris %>% parallelPCA(center=T, scale=T)

pca_exp = function(EXP, db=NULL, tcga_code=NULL, db_lvl=NULL, rank=10) {
  
  sd_zero = (sapply(EXP, sd)==0)
  if (any(sd_zero)) sprintf("# Num of SD=0 : %s", sum(sd_zero)) %>% print
  
  EXP_PCA = EXP[, !sd_zero] %>% prcomp(center=T, scale=T, rank.=rank)
  EXP_PCA = EXP_PCA$x %>% as.data.frame
  
  if (!is.null(db) & !is.null(tcga_code)) {
    EXP_PCA = EXP_PCA %>% mutate(Database=db, TCGA_Code=tcga_code) %>% 
      relocate(Database, TCGA_Code, .before=everything())
  }
  
  if (!is.null(db_lvl)) {
    EXP_PCA$Database = EXP_PCA$Database %>% factor(levels=db_lvl)
  }
  
  return(EXP_PCA)
}

identical(rownames(SANGER_RNA_Ent), rownames(SANGER_RNA_GSVA))   # T
identical(rownames(TCGA_RNA_Ent), rownames(TCGA_RNA_GSVA))       # T

db1 = rep("SANGER", nrow(SANGER_RNA_GSVA))
db2 = rep("TCGA", nrow(TCGA_RNA_GSVA))
db = c(db1, db2)

idx = match(rownames(SANGER_RNA_GSVA), Anno_Cells$SANGER_MODEL_ID)
tcga_code1 = Anno_Cells$TCGA_CODE[idx]
tcga_code1 = ifelse(!is.na(tcga_code1), tcga_code1, "UNCLASSIFIED")

all(rownames(TCGA_RNA_GSVA) %in% TCGA_Info$cases)   # T
idx = match(rownames(TCGA_RNA_GSVA), TCGA_Info$cases)
tcga_code2 = gsub("^TCGA-", "", TCGA_Info$project[idx])

setdiff(unique(tcga_code1), unique(tcga_code2)) %>% length   # 82
setdiff(unique(tcga_code2), unique(tcga_code1)) %>% length   # 8
tcga_code = c(tcga_code1, tcga_code2)

db_lvl = c("TCGA", "SANGER")


### 1. PCA Plot [GSVA]

identical(colnames(SANGER_RNA_GSVA), colnames(TCGA_RNA_GSVA))   # T
SANGER_TCGA_GSVA_PCA = rbind(SANGER_RNA_GSVA, TCGA_RNA_GSVA) %>% 
  pca_exp(db=db, tcga_code=tcga_code)

tcga_main = SANGER_TCGA_GSVA_PCA$TCGA_Code %>% 
  table %>% sort(decreasing=T) %>% head(5) %>% names

xlim = SANGER_TCGA_GSVA_PCA$PC1 %>% range
ylim = SANGER_TCGA_GSVA_PCA$PC2 %>% range
xlim = c(xlim[1]-0.1, xlim[2]+0.1)
ylim = c(ylim[1]-0.1, ylim[2]+0.1)

dir_ = "GSVA_Analysis/Batch Correction [TCGA]"
dir = mkdir(sprintf("../../processed_data/cell_data/BIOCARTA/%s", dir_))

file = sprintf("%s/PC2 Plot [GSVA]", dir)
SANGER_TCGA_GSVA_PCA %>% arrange(Database) %>% 
  plot_def(PC1, PC2, color=Database, main=file, 
           xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
           alpha=0.8, size=0.8, width=20, height=16, force_bold=T, save=T)

file = sprintf("%s/PC2 Plot [GSVA, SANGER, Top 5]", dir)
SANGER_TCGA_GSVA_PCA %>% 
  subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>% 
  plot_def(PC1, PC2, color=TCGA_Code, main=file, 
           xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
           alpha=0.8, size=2, width=20, height=16, force_bold=T, save=T)

file = sprintf("%s/PC2 Plot [GSVA, TCGA, Top 5]", dir)
SANGER_TCGA_GSVA_PCA %>% 
  subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>% 
  plot_def(PC1, PC2, color=TCGA_Code, main=file, 
           xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
           alpha=0.8, size=1.2, width=20, height=16, force_bold=T, save=T)


file = sprintf("%s/PCA_GSVA.csv", dir)
write.csv(SANGER_TCGA_GSVA_PCA, file=file, row.names=T)




##### 2-2. Parse TCGA Data [Clinical]

clin_query = list()
clin_drug = list()
clin_radia = list()
clin_patient = list()

project = TCGA_Info$project %>% unique   # 33
dir = "../../raw_data/TCGA"

for (proj in project) {
  # Make R object from the downloaded data
  query = GDCquery(project=proj, access="open",
                   data.category="Clinical", file.type="xml",
                   data.type="Clinical Supplement")
  GDCdownload(query, directory=dir)
  
  clin_query[[proj]] = query
  clin_drug[[proj]] = GDCprepare_clinic(query, clinical.info="drug", directory=dir)
  clin_radia[[proj]] = GDCprepare_clinic(query, clinical.info="radiation", directory=dir)
  clin_patient[[proj]] = GDCprepare_clinic(query, clinical.info="patient", directory=dir)
}

# If down for some reasons, replace for loop like below
# project_ = setdiff(project, names(clin_query))
# for (proj in project_)

col_common = function(df_list) {
  col_list = df_list %>% lapply(colnames)
  col = Reduce(intersect, col_list)
  col_num = col_list %>% sapply(length) %>% range
  sprintf("# Column range : [%s, %s]", col_num[1], col_num[2]) %>% print
  return(col)
}

col_drug = col_common(clin_drug)         # 24, 29
col_radia = col_common(clin_radia)       # 20, 20
col_patient = col_common(clin_patient)   # 44, 114

rbind_clin = function(df1, df2, col) {
  # All factor variables are converted into character ones
  df1[sapply(df1, is.factor)] = sapply(df1[sapply(df1, is.factor)], as.character)
  df2[sapply(df2, is.factor)] = sapply(df2[sapply(df2, is.factor)], as.character)
  df = rbind(df1[, col], df2[, col])
  df[df==""] = NA
  return(df)
}

rbind_drug = function(df1, df2) rbind_clin(df1, df2, col_drug)
rbind_radia = function(df1, df2) rbind_clin(df1, df2, col_radia)
rbind_patient = function(df1, df2) rbind_clin(df1, df2, col_patient)

TCGA_Drug = Reduce(rbind_drug, clin_drug)            # 13371 x 27
TCGA_Radia = Reduce(rbind_radia, clin_radia)         # 4230 x 27
TCGA_Patient = Reduce(rbind_patient, clin_patient)   # 12218 x 27

TCGA_Drug = TCGA_Drug %>% distinct %>% as.data.frame         # 13371 > 12571
TCGA_Radia = TCGA_Radia %>% distinct %>% as.data.frame       # 4230 > 3976
TCGA_Patient = TCGA_Patient %>% distinct %>% as.data.frame   # 12218 > 11167

dir = "../../processed_data/cell_data/TCGA"
file = sprintf("%s/TCGA_Clinical_Ori.RData", dir)
save(clin_query, clin_drug, clin_radia, clin_patient, 
     TCGA_Drug, TCGA_Radia, TCGA_Patient, file=file)



##### 2-2. Process TCGA Data [Clinical, Search PubChem CIDs]
# Process TCGA Response data personally
# Standardized drug names were obtained from GDISC [DrugCorrection1.csv]
# https://gdisc.bme.gatech.edu/cgi-bin/gdisc/tap5.cgi

# Principle of processing clinical data
# 1. Patients treated with only one drug
# 2. Patients with explicit and consistent responses

suppressMessages(library(webchem))

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/Anno_Drugs.csv", dir)
Anno_Drugs = read.csv(file)

dir = "../../raw_data/TCGA"

file = sprintf("%s/DrugCorrection1.csv", dir)
Anno_Drugs_TCGA1 = read.csv(file)
Anno_Drugs_TCGA1 = Anno_Drugs_TCGA1 %>% 
  mutate(OldName=trimws(OldName), 
         Correction=trimws(Correction)) %>% as.data.frame

file = sprintf("%s/DrugCorrection2.csv", dir)
Anno_Drugs_TCGA2 = read.csv(file, row.names=1)
Anno_Drugs_TCGA2 = Anno_Drugs_TCGA2 %>% 
  separate_rows(OldName, sep=",|, ") %>% 
  mutate(OldName=trimws(OldName), 
         Correction=trimws(Correction)) %>% as.data.frame

Anno_Drugs_TCGA = rbind(Anno_Drugs_TCGA1, Anno_Drugs_TCGA2) %>% distinct   # 940
drugs_tcga = Anno_Drugs_TCGA$Correction %>% unique   # 294
CIDs_TCGA = get_cid(drugs_tcga, from="name", domain="compound", match="all", verbose=T) 
CIDs_TCGA = CIDs_TCGA %>% as.data.frame
colnames(CIDs_TCGA) = c("Name", "CID")

# Manual Search
# Googling with drug names
# These drugs were not included in GDSC

# 'Tioguanine' [Thioguanine, CID 2723601]
# https://pubchem.ncbi.nlm.nih.gov/compound/2723601#section=DSSTox-Substance-ID
# 
# CEP -11981 [CEP-11981, CID 11751922]
# https://pubchem.ncbi.nlm.nih.gov/compound/11751922
# 
# LY2228820 [CID 11570805]
# https://pubchem.ncbi.nlm.nih.gov/compound/11570805
# 
# Porfimer [CID 3086257]
# https://pubchem.ncbi.nlm.nih.gov/compound/3086257
# 
# Tiomolibdate [CID 5245480]
# https://pubchem.ncbi.nlm.nih.gov/compound/5245480
#
# Vamydex [CID 135430970]
# https://ncithesaurus.nci.nih.gov/ncitbrowser/pages/concept_details.jsf
# https://pubchem.ncbi.nlm.nih.gov/compound/135430970

names_temp = c("Thioguanine", "CEP-11981", "LY2228820", "Porfimer", "Tiomolibdate", "Vamydex")
cids_temp = c(2723601, 11751922, 11570805, 3086257, 5245480, 135430970)
CIDs_TCGA_Temp = data.frame(Name=names_temp, CID=cids_temp)
CIDs_TCGA = CIDs_TCGA %>% rbind(CIDs_TCGA_Temp)
CIDs_TCGA$CID[CIDs_TCGA$Name=="NOS"] = NA
# NOS = "Not Otherwise Specified"

Anno_Drugs_TCGA$Correction[Anno_Drugs_TCGA$Correction=="'Tioguanine'"] = "Thioguanine"
Anno_Drugs_TCGA$Correction[Anno_Drugs_TCGA$Correction=="CEP -11981"] = "CEP-11981"

idx = match(Anno_Drugs_TCGA$Correction, CIDs_TCGA$Name)
Anno_Drugs_TCGA$CID = CIDs_TCGA$CID[idx]

# Allocate CIDs from each drug in TCGA Response Data
idx = match(TCGA_Drug$drug_name, Anno_Drugs_TCGA$OldName)
col = c("bcr_patient_barcode", "drug_name", "measure_of_response", "project")

TCGA_Resp = TCGA_Drug[, col] %>% 
  mutate(drug_cid=Anno_Drugs_TCGA$CID[idx]) %>% 
  relocate(drug_cid, .after=drug_name) %>% as.data.frame

# Patients with unknown or no drug treatments were excluded
TCGA_Resp = TCGA_Resp %>% subset(!is.na(drug_cid))   # 12571 > 11203

# Patients with multiple drug treatments were excluded
TCGA_Resp = TCGA_Resp %>% group_by(bcr_patient_barcode) %>% 
  filter(n_distinct(drug_cid)==1) %>% as.data.frame   # 11203 > 1598

# Patients with unknown or no responses were excluded
TCGA_Resp = TCGA_Resp %>% subset(!is.na(measure_of_response))   # 1598 > 515

# Patients with multiple responses were excluded
TCGA_Resp = TCGA_Resp %>% group_by(bcr_patient_barcode) %>% 
  filter(n_distinct(measure_of_response)==1) %>% as.data.frame    # 515 > 494

# Duplicated information is removed
TCGA_Resp = TCGA_Resp %>% 
  distinct(bcr_patient_barcode, drug_cid, measure_of_response, .keep_all=T)   # 494 > 482

# Drugs at least 10 patients were left
TCGA_Resp$drug_cid %>% unique %>% length   # 46
TCGA_Resp = TCGA_Resp %>% group_by(drug_cid) %>% 
  filter(length(unique(bcr_patient_barcode))>=10) %>% as.data.frame   # 482 > 414

TCGA_Resp$drug_cid %>% unique %>% length   # 11
colnames(TCGA_Resp) = c("Patient", "Drug_Name", "Drug_CID", "Response", "TCGA_Code")
TCGA_Resp$TCGA_Code = gsub("^TCGA-", "", TCGA_Resp$TCGA_Code) 

# Patients not included in RNA-Seq were excluded
# Ex. TCGA[Project]-02[Tissue]-0001[Patient]-01C[Sample]-01D[Portion]-0182[Plate]-01[Center]
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode

pat_rna = rownames(TCGA_RNA_Ent) %>% strsplit("-")
pat_rna = pat_rna %>% sapply(function(x) paste(x[1:3], collapse="-"))
pat_rna %>% unique %>% length   # 10385 [from 11274]

sum(TCGA_Resp$Patient %in% pat_rna)   # 399 [from 414]
TCGA_Resp = TCGA_Resp %>% subset(Patient %in% pat_rna) %>% 
  group_by(Drug_CID) %>% filter(n_distinct(Patient)>=10) %>% as.data.frame   # 414 > 399

# Match RNA-Seq Samples-Patients
col = c("cases", "sample_type", "cases.submitter_id")
TCGA_Resp = left_join(TCGA_Resp, TCGA_Info[, col], by=c("Patient"="cases.submitter_id"))

TCGA_Resp = TCGA_Resp %>% 
  mutate(Sample=cases, Sample_Type=sample_type) %>% 
  relocate(Patient, Sample, Sample_Type, TCGA_Code, .before=everything()) %>% 
  select(subset=-c(cases, sample_type)) %>% as.data.frame   # 414

class_resp = c("Complete Response", "Partial Response")
TCGA_Resp_Sum = TCGA_Resp %>% group_by(Drug_CID) %>% 
  summarise(Resp_C1=sum((Response %in% class_resp)), 
            Non_C1=sum(!(Response %in% class_resp)), 
            Resp_C2=sum((Response %in% class_resp[1])), 
            Non_C2=sum(!(Response %in% class_resp[1])), 
            Total=Resp_C1+Non_C1) %>% as.data.frame

TCGA_Resp_Sum$Total %>% sum   # 399
cids_tcga = TCGA_Resp$Drug_CID %>% unique                # 11
cids_gdsc = Anno_Drugs$Drug_CID %>% na.omit %>% unique   # 432
cids_tcga[!(cids_tcga %in% cids_gdsc)]  

# Manual searching from GDSC Drug Annotation
# Could not found all those 3 drugs from GDSC
# 426756  Carboplatin
# 60953   Capecitabine, Xeloda
# 657181  Eligard, Lupron, Leuprorelin


# Get SMILES of TCGA Drugs from PubChem
# https://pubchem.ncbi.nlm.nih.gov

col = c("Correction", "CID")
TCGA_Drug_Info = Anno_Drugs_TCGA[, col] %>% 
  subset(CID %in% unique(TCGA_Resp$Drug_CID)) %>% distinct

colnames(TCGA_Drug_Info) = c("Drug_Name", "Drug_CID")
TCGA_Drug_Info = TCGA_Drug_Info %>% arrange(Drug_Name)
TCGA_Drug_Info = left_join(TCGA_Drug_Info, TCGA_Resp_Sum, by="Drug_CID")


dir = "../../processed_data/cell_data/TCGA"
file = sprintf("%s/TCGA_Drug_Info.csv", dir)
write.csv(TCGA_Drug_Info, row.names=F, file=file)

file = sprintf("%s/TCGA_Drug_SMILES.csv", dir)
TCGA_Drug_SMILES = read.csv(file)

idx = match(TCGA_Drug_Info$Drug_CID, TCGA_Drug_SMILES$cid)
TCGA_Drug_Info = TCGA_Drug_Info %>% 
  mutate(Inchi_Key = TCGA_Drug_SMILES$inchikey[idx],
         SMILES_Can = TCGA_Drug_SMILES$canonicalsmiles[idx])

TCGA_Drug_Info = TCGA_Drug_Info %>% 
  mutate(GDSC=(Drug_CID %in% Anno_Drugs$Drug_CID)) %>% 
  relocate(GDSC, .after=Drug_CID) %>% as.data.frame


# # Average row-wise using a certain column 
# a = 1:24 %>% matrix(ncol=4) %>% as.data.frame
# a$X = sprintf("X%s", c(1, 1, 2, 3, 3, 3))
# a %>% group_by(X) %>% summarise_if(is.numeric, ~ mean(.x, na.rm=T)) %>% as.data.frame
# a %>% group_by(X) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm=T))) %>% as.data.frame


dir = "../../processed_data/cell_data/TCGA"

file = sprintf("%s/TCGA_Response.csv", dir)
write.csv(TCGA_Resp, row.names=F, file=file)

file = sprintf("%s/TCGA_Drug_Info.csv", dir)
write.csv(TCGA_Drug_Info, row.names=F, file=file)

file = sprintf("%s/Anno_Drugs_TCGA.csv", dir)
write.csv(Anno_Drugs_TCGA, row.names=F, file=file)




##### cf. Prepare benchmark tests

unify_genes = function(Omics, genes, value=0, cell_row=T) {
  if (!is.data.frame(Omics)) Omics = as.data.frame(Omics)
  genes_omics = if (cell_row) colnames(Omics) else rownames(Omics)
  genes_na = genes[!(genes %in% genes_omics)]
  
  if (length(genes_na)==length(genes)) {
    print("No gene matching, check the categories of gene names")
  } else {
    print(sprintf("The number of genes missing : %s", len(genes_na)))
  }
  
  if (cell_row) {
    Omics[, genes_na] = value
    Omics = Omics[, genes]
  } else {
    Omics[genes_na, ] = value
    Omics = Omics[genes, ]
  }
  return(Omics)
}

combat_def = function(Omics_List, tcga_list, db_list, 
                      ref_batch=1, par_prior=T, cell_row=T, cores=T, ...) {
  
  # Omics_List [List of data.frame] : Omics_Source, Omics_Target1, Omics_Target2, ...
  # tcga_list [List of character] : tcga_source, tcga_target1, tcga_target2, ...
  # db_list [Character] : db_source, db_target1, db_target2, ... 
  
  suppressMessages(library(sva))
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  BPPARAM = BiocParallel::MulticoreParam(workers=cores)
  
  t_def = function(df) as.data.frame(t(df))
  if (cell_row) Omics_List = lapply(Omics_List, t_def)
  cell_list = lapply(Omics_List, colnames)
  
  n_data = sapply(Omics_List, ncol)
  Omics_List = Reduce(cbind, Omics_List)
  
  if (is.null(tcga_list)) {
    Mod = NULL
  } else {
    Pheno = data.frame(Tissue=unlist(tcga_list))
    rownames(Pheno) = colnames(Omics_List)
    Mod = model.matrix(~as.factor(Tissue), data=Pheno)
  }
  
  # Reference-batch ComBat with covariates
  batch = rep(1:length(db_list), n_data)
  Omics_List = ComBat(dat=Omics_List, batch=batch, mod=Mod, 
                      par.prior=par_prior, ref.batch=ref_batch, BPPARAM=BPPARAM, ...)
  
  Omics_List = t_def(Omics_List)
  Omics_List = lapply(cell_list, function(x) Omics_List[x, ])
  names(Omics_List) = db_list
  return(Omics_List)
}

pca_exp = function(EXP, db=NULL, tcga_code=NULL, db_lvl=NULL, rank=10) {
  
  sd_zero = (sapply(EXP, sd)==0)
  if (any(sd_zero)) sprintf("# Num of SD=0 : %s", sum(sd_zero)) %>% print
  
  EXP_PCA = EXP[, !sd_zero] %>% prcomp(center=T, scale=T, rank.=rank)
  EXP_PCA = EXP_PCA$x %>% as.data.frame
  
  if (!is.null(tcga_code)) {
    EXP_PCA = EXP_PCA %>% mutate(TCGA_Code=tcga_code)
  }
  
  if (!is.null(db)) {
    EXP_PCA = EXP_PCA %>% mutate(Database=db) %>% 
      relocate(Database, .before=everything())
  }
  
  if (!is.null(db_lvl)) {
    EXP_PCA$Database = EXP_PCA$Database %>% factor(levels=db_lvl)
  }
  
  return(EXP_PCA)
}

plot_pca_batch = function(Omics_List, db_list=NULL, tcga_list=NULL, tcga_color=NULL, 
                          main=NULL, add_db=NULL, add_tcga=NULL, legend="TCGA Code", 
                          size=0.5, width=16.5, width2=16.5, height=13.5, height2=13.5, 
                          save=T, save_svg=F, return_pca=T, ...) {
  
  # Omics_List [List of data.frame] : Omics_Source, Omics_Target1, Omics_Target2, ...
  # tcga_list [List of character] : tcga_source, tcga_target1, tcga_target2, ...
  # db_list [Character] : db_source, db_target1, db_target2, ... 
  
  # Annotation Label [TCGA Code]
  if (!is.null(tcga_list)) {
    tcga_code = tcga_list %>% unlist
    cond = identical(unname(sapply(tcga_list, length)), unname(sapply(Omics_List, nrow)))
    if (!cond) stop("Check the TCGA Codes...")
  } else {
    tcga_code = NULL
  }
  
  # Annotation Label [DB]
  n_samples = sapply(Omics_List, nrow)
  db = rep(db_list, n_samples)
  
  # Compress Omics Data by PCA
  col = Reduce(intersect, lapply(Omics_List, colnames))
  rbind_ = function(df1, df2) rbind(df1[, col], df2[, col])
  
  Omics_PCA = Reduce(rbind_, Omics_List) %>% 
    pca_exp(db=db, tcga_code=tcga_code, db_lvl=db_list)
  
  # Set X and Y Ranges
  xlim = Omics_PCA$PC1 %>% range
  ylim = Omics_PCA$PC2 %>% range
  xlim = c(xlim[1]-0.1, xlim[2]+0.1)
  ylim = c(ylim[1]-0.1, ylim[2]+0.1)
  
  # Plot PCA [Database]
  Omics_PCA %>% arrange(Database) %>% 
    plot_def(PC1, PC2, color=Database, shape=Database, main=main, add=add_db,
             legend="Database", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
             size=size, width=width, height=height, save=save, save_svg=save_svg, ...)
  
  if (!is.null(tcga_color)) {
    for (i in 1:length(db_list)) {
      add_tcga = c(add_tcga, list(labs(shape=legend)))
      main_ = sprintf("%s (%s colored by %s)", main, db_list[i], legend)
      Omics_PCA %>% subset(Database==db_list[i] & TCGA_Code %in% tcga_color) %>% 
        mutate(TCGA_Code = factor(TCGA_Code, levels=tcga_color)) %>% 
        plot_def(PC1, PC2, color=TCGA_Code, shape=TCGA_Code, main=main_, add=add_tcga,
                 legend=legend, xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
                 size=size, width=width2, height=height2, save=save, save_svg=save_svg, ...)
    }
  }
  
  if (return_pca) return(Omics_PCA)
}

plot_pca_tcga = function(Omics_List, Omics_List_CB, db_list=NULL, tcga_list=NULL, 
                         main=NULL, idx_tcga=NULL, train_db="SANGER", 
                         width=18.6, width2=18.6, height=13.5) {
  
  if (!is.null(idx_tcga)) {
    tcga_list[[2]] = tcga_list[[2]][idx_tcga]
    Omics_List[[2]] = Omics_List[[2]][idx_tcga, ]
    Omics_List_CB[[2]] = Omics_List_CB[[2]][idx_tcga, ]
  }
  
  tissue_top7 = tcga_list[[2]] %>% 
    table %>% sort(decreasing=T) %>% head(7) %>% names
  
  shape_db = c(21, 24)
  color_db = c("brown1", "royalblue1")
  names(color_db) = c(train_db, "TCGA")
  
  shape_tcga = 0:6
  color_tcga = c("brown1", "coral", "gold", "seagreen3", "royalblue1", "mediumorchid", "black")
  
  add_db = list(scale_color_manual(values=color_db), 
                scale_shape_manual(values=shape_db), 
                theme(legend.key.size=unit(1, "cm")), 
                guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
  
  add_tcga = list(scale_color_manual(values=color_tcga), 
                  scale_shape_manual(values=shape_tcga), 
                  theme(legend.key.size=unit(1, "cm")), 
                  guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
  
  PCA = Omics_List %>% plot_pca_batch(
    db_list, tcga_list, tcga_color=tissue_top7,
    main=main[1], size=1.5, alpha=0.25, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, width2=width2, height=height,
    add_db=add_db, add_tcga=add_tcga, save=T, save_svg=T)
  
  PCA_CB = Omics_List_CB %>% plot_pca_batch(
    db_list, tcga_list, tcga_color=tissue_top7, 
    main=main[2], size=1.5, alpha=0.25, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, width2=width2, height=height, 
    add_db=add_db, add_tcga=add_tcga, save=T, save_svg=T)
  
  PCA_List = list(PCA, PCA_CB)
  names(PCA_List) = c("ComBat_X", "ComBat_O")
  return(PCA_List)
}

# Set UNCLASSIFIED
Anno_Cells$TCGA_CODE %>% is.na %>% sum   # 487
Anno_Cells$TCGA_CODE[is.na(Anno_Cells$TCGA_CODE)] = "UNCLASSIFIED"

category = c("ComBat X", "ComBat O")
dir_pca = mkdir("Batch Effects with PCA Plot")


### GCNPath
idx1 = match(rownames(SANGER_RNA_GSVA), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(TCGA_RNA_GSVA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(SANGER_RNA_GSVA, TCGA_RNA_GSVA)
TCGA_RNA_GSVA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
file = sprintf("%s/TCGA_RNA_GSVA_ComBat.csv", dir)
fwrite(TCGA_RNA_GSVA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [GCNPath, %s]", dir_pca, category)
idx_tcga = rownames(TCGA_RNA_GSVA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(SANGER_RNA_GSVA, TCGA_RNA_GSVA)

PCA_GCNPath = plot_pca_tcga(Temp_List, TCGA_RNA_GSVA_, db_list, tcga_list, 
                            main=main, idx_tcga=idx_tcga, train_db="SANGER")


### Genes
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
SANGER_RNA_TPM = fread_def(file, col_numeric=T)
geneset_ = Reduce(intersect, list(geneset, colnames(SANGER_RNA_TPM), colnames(TCGA_RNA_Ent)))   # 1503

idx1 = match(rownames(SANGER_RNA_TPM), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(TCGA_RNA_Ent), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(SANGER_RNA_TPM[, geneset_], TCGA_RNA_Ent[, geneset_])
TCGA_RNA_ = combat_def(exp_list, tcga_list, db_list)

# dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
# file = sprintf("%s/TCGA_RNA_GSVA_ComBat.csv", dir)
# fwrite(TCGA_RNA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [BIOCARTA Genes, %s]", dir_pca, category)
idx_tcga = rownames(TCGA_RNA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(SANGER_RNA_TPM[, geneset_], TCGA_RNA_Ent[, geneset_])

PCA_Genes = plot_pca_tcga(Temp_List, TCGA_RNA_, db_list, tcga_list, 
                          main=main, idx_tcga=idx_tcga, train_db="SANGER")



### PaccMann
dir = "../../benchmark_test/PaccMann/data/gene_expression"
file = sprintf("%s/gdsc-rnaseq_gene-expression.csv", dir)
EXP_PaccMann = fread_def(file)

library(reticulate)
pickle <- import("pickle")
builtins <- import_builtins()

dir = "../../benchmark_test/PaccMann/data"
file = sprintf("%s/2128_genes.pkl", dir)
con = builtins$open(file, "rb")
genes_paccmann = pickle$load(con)
con$close()

TCGA_RNA_Sym_ = 2**TCGA_RNA_Sym - 1
TCGA_RNA_Sym_ = TCGA_RNA_Sym_ %>% asinh %>% as.data.frame
EXP_TCGA_PaccMann = TCGA_RNA_Sym_ %>% unify_genes(genes_paccmann)   # 2128 [NA 40]

dir = "../../benchmark_test/PaccMann/_data"
file = sprintf("%s/EXP_TCGA.csv", dir)
fwrite(EXP_TCGA_PaccMann, file=file, row.names=T, col.names=T)


### PaccMann_SANGER
dir = "../../benchmark_test/PaccMann_SANGER/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_PaccMann_SG = fread_def(file)

genes_paccmann_sg = EXP_PaccMann_SG %>% colnames
TCGA_RNA_Sym_ = 2**TCGA_RNA_Sym - 1
TCGA_RNA_Sym_ = TCGA_RNA_Sym_ %>% asinh %>% as.data.frame
EXP_PaccMann_SG_TCGA = TCGA_RNA_Sym_ %>% unify_genes(genes_paccmann_sg)   # 2093 [NA 5]

dir = "../../benchmark_test/PaccMann_SANGER/_data"
file = sprintf("%s/EXP_TCGA.csv", dir)
fwrite(EXP_PaccMann_SG_TCGA, file=file, row.names=T, col.names=T)


idx1 = match(rownames(EXP_PaccMann_SG), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(EXP_PaccMann_SG_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_PaccMann_SG, EXP_PaccMann_SG_TCGA)
EXP_PaccMann_SG_TCGA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/PaccMann_SANGER/_data"
file = sprintf("%s/EXP_TCGA_ComBat.csv", dir)
fwrite(EXP_PaccMann_SG_TCGA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [PaccMann_SANGER, %s]", dir_pca, category)
idx_tcga = rownames(EXP_PaccMann_SG_TCGA_[[2]]) %in% TCGA_Resp$Sample
Temp_List = list(EXP_PaccMann_SG, EXP_PaccMann_SG_TCGA)

PCA_PaccMann_SG = plot_pca_tcga(Temp_List, EXP_PaccMann_SG_TCGA_, 
                                db_list, tcga_list, main=main, 
                                idx_tcga=idx_tcga, train_db="SANGER")


### HiDRA
scale_ = function(df) {
  df = df %>% scale
  df[is.na(df)] = 0
  df[is.infinite(df)] = 0
  df = df %>% as.data.frame
  return(df)
}

dir = "../../benchmark_test/HiDRA/_data/ProcessedFile"
file_list = sprintf("%s/%s.csv", dir, 0:185)

EXP_HiDRA = list()
Gene_List = list()
for (i in 1:length(file_list)) {
  EXP_HiDRA[[i]] = fread_def(file_list[i])
  Gene_List[[i]] = EXP_HiDRA[[i]] %>% colnames
}

genes_hidra = Gene_List %>% unlist %>% unique   # 4590
EXP_HiDRA_TCGA = TCGA_RNA_Sym %>% t %>% scale %>% t %>% as.data.frame
EXP_HiDRA_TCGA = EXP_HiDRA_TCGA %>% unify_genes(genes_hidra)

dir = mkdir("../../benchmark_test/HiDRA/_data/ProcessedFile_TCGA")
for (i in 1:length(file_list)) {
  file = sprintf("%s/%s.csv", dir, i-1)
  Temp = EXP_HiDRA_TCGA[, Gene_List[[i]]]
  fwrite(Temp, file=file, row.names=T, col.names=T)
}


EXP_HiDRA = Reduce(cbind, EXP_HiDRA)[, genes_hidra]
idx1 = match(rownames(EXP_HiDRA), Anno_Cells$COSMIC_ID)
idx2 = match(rownames(EXP_HiDRA_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("GDSC", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_HiDRA, EXP_HiDRA_TCGA)
EXP_HiDRA_TCGA_ = combat_def(exp_list, tcga_list, db_list)

dir = mkdir("../../benchmark_test/HiDRA/_data/ProcessedFile_TCGA_ComBat")
for (i in 1:length(file_list)) {
  file = sprintf("%s/%s.csv", dir, i-1)
  Temp = EXP_HiDRA_TCGA_[[2]][, Gene_List[[i]]]
  fwrite(Temp, file=file, row.names=T, col.names=T)
}


main = sprintf("%s/PCA Plot [HiDRA, %s]", dir_pca, category)
idx_tcga = rownames(EXP_HiDRA_TCGA_[[2]]) %in% TCGA_Resp$Sample
Temp_List = list(EXP_HiDRA, EXP_HiDRA_TCGA)

PCA_HiDRA = plot_pca_tcga(Temp_List, EXP_HiDRA_TCGA_, 
                          db_list, tcga_list, main=main, 
                          idx_tcga=idx_tcga, train_db="GDSC")


### TGDRP
dir = "../../benchmark_test/TGSA/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_TGSA = read.csv(file, row.names=1, check.names=F)

genes_tgsa = EXP_TGSA %>% colnames
EXP_TGSA_TCGA = TCGA_RNA_Ent %>% unify_genes(genes_tgsa)   # 706 [NA 0]

# dir = mkdir("../../benchmark_test/TGSA/_data_tcga")
# file = sprintf("%s/EXP.csv", dir)
# fwrite(EXP_TGSA_TCGA, file=file, row.names=T, col.names=T)


# Get Mutation & CNV Data
gdc_code = getGDCprojects()$project_id
tcga_code = grep("TCGA", gdc_code, value=T)   # 33

TCGA_Resp_Ori = TCGA_Resp
patient = TCGA_Resp$Patient %>% unique        # 399

TCGA_Resp = TCGA_Resp_Ori %>% 
  dplyr::rename(Sample_RNA = Sample) %>% 
  mutate(Sample = gsub("^(([^-]+-){3}[^-]+).*", "\\1", Sample_RNA)) %>% 
  relocate(Sample, .before=Sample_RNA) %>% as.data.frame

query_tcga_mut = GDCquery(project=tcga_code, access="open", 
                          data.category="Simple Nucleotide Variation", 
                          data.type="Masked Somatic Mutation", barcode=patient)

query_tcga_cnv = GDCquery(project=tcga_code, access="open", 
                          data.category="Copy Number Variation", 
                          data.type="Gene Level Copy Number", barcode=patient)

# Download from GDC repository
dir = "../../raw_data/TCGA"
GDCdownload(query_tcga_mut, directory=dir)
# Make R object from the downloaded data
data_mut = GDCprepare(query_tcga_mut, directory=dir)   # 88264 x 140
data_mut$Tumor_Sample_Barcode %>% unique %>% length    # 380

MUT_TGSA_TCGA = data_mut %>% subset(Entrez_Gene_Id %in% genes_tgsa)
MUT_TGSA_TCGA = MUT_TGSA_TCGA %>% as.data.frame %>% 
  reshape2::acast(Tumor_Sample_Barcode~Entrez_Gene_Id, 
                  value.var="Entrez_Gene_Id", fun.aggregate=any, fill=0)

MUT_TGSA_TCGA = MUT_TGSA_TCGA %>% unify_genes(genes_tgsa)   # [NA] 23

# Download from GDC repository
# Unfortunately, TCGA cannot provide gene-level log2 CN-ratio...
# The gene-level CNV values are discrete

dir = "../../raw_data/TCGA"
GDCdownload(query_tcga_cnv, directory=dir)
# Make R object from the downloaded data

# data_cnv = GDCprepare(query_tcga_cnv, directory=dir)
# Error in GDCprepare(query_tcga_cnv, directory = dir) : 
# There are samples duplicated. We will not be able to prepare it

query_tcga_cnv$results[[1]]$cases %>% unique %>% length
query_tcga_cnv$results[[1]] = query_tcga_cnv$results[[1]] %>% distinct(cases, .keep_all=T)
query_tcga_cnv$results[[1]] = query_tcga_cnv$results[[1]] %>% distinct(sample.submitter_id, .keep_all=T)
query_tcga_cnv$results[[1]] = query_tcga_cnv$results[[1]] %>% subset(sample.submitter_id %in% TCGA_Resp$Sample)

data_cnv = GDCprepare(query_tcga_cnv, directory=dir)
TCGA_CNV = assay(data_cnv) %>% t %>% as.data.frame   # 386 x 60623
TCGA_CNV_Ent = TCGA_CNV

ens = gsub("\\.[0-9]+", "", colnames(TCGA_CNV_Ent))
idx = match(ens, Anno_Genes_TCGA$ENSEMBL)
col = Anno_Genes_TCGA$ENTREZID[idx]

col %>% na.omit %>% duplicated %>% sum   # 156
idx_no = is.na(col) | duplicated(col)
TCGA_CNV_Ent = TCGA_CNV_Ent[, !idx_no]   # 60623 > 40592
colnames(TCGA_CNV_Ent) = col[!idx_no]

# CNV in TCGA were categorized
# Center CNV values for deep learning
# Train TGDRP & TGSA with there categorized CNV with retrain_total.sh
# And test these models with TCGA dataset with test_tcga.sh

# Val = round( 2 * 2^log2(C/Ploidy) )
# Equal to 2 * CN-Ratio...
# if Val == 0: Category = 'Deletion'
# if Val == 1: Category = 'Loss'
# if Val == 2: Category = 'Neutral'
# if Val == 3: Category = 'Gain'
# if Val >= 4: Category = 'Amplification'

CNV_TGSA_TCGA = TCGA_CNV_Ent %>% unify_genes(genes_tgsa, value=2)   # [NA 0]
CNV_TGSA_TCGA = CNV_TGSA_TCGA - 2
CNV_TGSA_TCGA[CNV_TGSA_TCGA>=2] = 2
CNV_TGSA_TCGA[CNV_TGSA_TCGA<=-2] = -2

CNV_TGSA_TCGA %>% is.na %>% sum   # 16142
CNV_TGSA_TCGA[is.na(CNV_TGSA_TCGA)] = 0
CNV_TGSA_TCGA %>% unlist %>% table
#  -2     -1       0      1      2 
# 144  20022  164404  37784  50162


sample_tgsa_exp = gsub("^(([^-]+-){3}[^-]+).*", "\\1", rownames(EXP_TGSA_TCGA))
sample_tgsa_mut = gsub("^(([^-]+-){3}[^-]+).*", "\\1", rownames(MUT_TGSA_TCGA))
sample_tgsa_cnv = gsub("^(([^-]+-){3}[^-]+).*", "\\1", rownames(CNV_TGSA_TCGA))
sample_tgsa = Reduce(intersect, list(sample_tgsa_exp, sample_tgsa_mut, sample_tgsa_cnv))   # 363
sample_tgsa %>% unique %>% length   # 363

EXP_TGSA_TCGA = EXP_TGSA_TCGA[sample_tgsa_exp %in% sample_tgsa, ]   # 363 x 706
MUT_TGSA_TCGA = MUT_TGSA_TCGA[sample_tgsa_mut %in% sample_tgsa, ]   # 363 x 706
CNV_TGSA_TCGA = CNV_TGSA_TCGA[sample_tgsa_cnv %in% sample_tgsa, ]   # 363 x 706

idx1 = match(sample_tgsa, sample_tgsa_exp[sample_tgsa_exp %in% sample_tgsa])
idx2 = match(sample_tgsa, sample_tgsa_mut[sample_tgsa_mut %in% sample_tgsa])
idx3 = match(sample_tgsa, sample_tgsa_cnv[sample_tgsa_cnv %in% sample_tgsa])

EXP_TGSA_TCGA = EXP_TGSA_TCGA[idx1, ]
MUT_TGSA_TCGA = MUT_TGSA_TCGA[idx2, ]
CNV_TGSA_TCGA = CNV_TGSA_TCGA[idx3, ]


idx1 = match(rownames(EXP_TGSA), Anno_Cells$BROAD_ID)
idx2 = match(rownames(EXP_TGSA_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("CCLE", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_TGSA, EXP_TGSA_TCGA)
EXP_TGSA_TCGA_ = combat_def(exp_list, tcga_list, db_list)


main = sprintf("%s/PCA Plot [TGSA, %s]", dir_pca, category)
idx_tcga = rownames(EXP_TGSA_TCGA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(EXP_TGSA %>% scale_, EXP_TGSA_TCGA %>% scale_)

PCA_TGSA = plot_pca_tcga(Temp_List, EXP_TGSA_TCGA_, 
                         db_list, tcga_list, main=main, 
                         idx_tcga=idx_tcga, train_db="CCLE")

dir = mkdir("../../benchmark_test/TGSA/_data_tcga")
file1 = sprintf("%s/EXP.csv", dir)
file2 = sprintf("%s/MUT.csv", dir)
file3 = sprintf("%s/CNV_Category.csv", dir)
file4 = sprintf("%s/EXP_ComBat.csv", dir)

fwrite(EXP_TGSA_TCGA, file=file1, row.names=T, col.names=T)
fwrite(MUT_TGSA_TCGA, file=file2, row.names=T, col.names=T)
fwrite(CNV_TGSA_TCGA, file=file3, row.names=T, col.names=T)
fwrite(EXP_TGSA_TCGA_[[2]], file=file4, row.names=T, col.names=T)


# CCLE DepMap do not have categorized CNV data
# Instead, we reverse-transformed and categorized log2(CN-Ratio + 1) data
dir = "../../benchmark_test/TGSA/_data"
file = sprintf("%s/CNV.csv", dir)
CNV_TGSA_CCLE = read.csv(file, row.names=1, check.names=F)

CNV_TGSA_CCLE = 2**CNV_TGSA_CCLE - 1
CNV_TGSA_CCLE %>% unlist %>% range    # 1.926979e-10 5.958854e+01
CNV_TGSA_CCLE %>% unlist %>% mean     # 1.046076
CNV_TGSA_CCLE %>% unlist %>% median   # 1.007339

CNV_TGSA_CCLE = round(CNV_TGSA_CCLE * 2)
CNV_TGSA_CCLE %>% unlist %>% table
#    0       1       2       3      4  ...      
# 1017  100700  739051  126007  14796  ...

CNV_TGSA_CCLE = CNV_TGSA_CCLE - 2
CNV_TGSA_CCLE[CNV_TGSA_CCLE>=2] = 2
CNV_TGSA_CCLE %>% unlist %>% table
#   -2      -1       0       1      2 
# 1017  100700  739051  126007  20919 

dir = "../../benchmark_test/TGSA/_data"
file = sprintf("%s/CNV_Category.csv", dir)
write.csv(CNV_TGSA_CCLE, file=file, row.names=T)


### TGDRP_SANGER
dir = "../../benchmark_test/TGSA_SANGER/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_TGSA_SG = read.csv(file, row.names=1, check.names=F)
genes_tgsa_sg = EXP_TGSA_SG %>% colnames
identical(genes_tgsa, genes_tgsa_sg)   # T


idx1 = match(rownames(EXP_TGSA_SG), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(EXP_TGSA_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_TGSA_SG, EXP_TGSA_TCGA)
EXP_TGSA_SG_TCGA_ = combat_def(exp_list, tcga_list, db_list)


main = sprintf("%s/PCA Plot [TGSA_SANGER, %s]", dir_pca, category)
idx_tcga = rownames(EXP_TGSA_SG_TCGA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(EXP_TGSA_SG %>% scale_, EXP_TGSA_TCGA %>% scale_)

PCA_TGSA_SG = plot_pca_tcga(Temp_List, EXP_TGSA_SG_TCGA_, 
                            db_list, tcga_list, main=main, 
                            idx_tcga=idx_tcga, train_db="SANGER")


dir = mkdir("../../benchmark_test/TGSA_SANGER/_data_tcga")
file1 = sprintf("%s/EXP.csv", dir)
file2 = sprintf("%s/MUT.csv", dir)
file3 = sprintf("%s/CNV_Category.csv", dir)
file4 = sprintf("%s/EXP_ComBat.csv", dir)

fwrite(EXP_TGSA_TCGA, file=file1, row.names=T, col.names=T)
fwrite(MUT_TGSA_TCGA, file=file2, row.names=T, col.names=T)
fwrite(CNV_TGSA_TCGA, file=file3, row.names=T, col.names=T)
fwrite(EXP_TGSA_SG_TCGA_[[2]], file=file4, row.names=T, col.names=T)


# SANGER Passports also have categorized CNV data
dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/WES_pureCN_CNV_genes_cn_category_20221213.csv", dir)
CNV_TGSA_SG = fread(file, header=T, na.strings=c("NA", ""))

cells = CNV_TGSA_SG[1, 2:ncol(CNV_TGSA_SG)] %>% as.character              # 1357 [1253]
genes = CNV_TGSA_SG[4:nrow(CNV_TGSA_SG), 1] %>% unlist %>% as.character   # 18348 [18348]

CNV_TGSA_SG = CNV_TGSA_SG[4:nrow(CNV_TGSA_SG), 2:ncol(CNV_TGSA_SG)] %>% as.data.frame
colnames(CNV_TGSA_SG) = cells
rownames(CNV_TGSA_SG) = genes

CNV_TGSA_SG = CNV_TGSA_SG[, !duplicated(cells)] %>% t %>% as.data.frame

idx = match(colnames(CNV_TGSA_SG), Anno_Genes$HGNC_SYMBOL)
genes_ = Anno_Genes$ENTREZ_ID[idx]   # 18348
sum(genes_tgsa %in% genes_)          # 688 / 706 [NA 18]
CNV_TGSA_SG = CNV_TGSA_SG[, !is.na(genes_)] %>% setNames(na.omit(genes_))   # 1253 x 17975
CNV_TGSA_SG = CNV_TGSA_SG %>% unify_genes(genes_tgsa, "Neutral")            # 1253 x 706 [NA] 18
CNV_TGSA_SG[is.na(CNV_TGSA_SG)] = "Neutral"

CNV_TGSA_SG[CNV_TGSA_SG=="Deletion"] = -2
CNV_TGSA_SG[CNV_TGSA_SG=="Loss"] = -1
CNV_TGSA_SG[CNV_TGSA_SG=="Neutral"] = 0
CNV_TGSA_SG[CNV_TGSA_SG=="Gain"] = 1
CNV_TGSA_SG[CNV_TGSA_SG=="Amplification"] = 2
CNV_TGSA_SG[] = lapply(CNV_TGSA_SG, as.numeric)

all(rownames(EXP_TGSA_SG) %in% rownames(CNV_TGSA_SG))    # T
CNV_TGSA_SG = CNV_TGSA_SG[rownames(EXP_TGSA_SG), ]   # 1183 x 706

dir = "../../benchmark_test/TGSA_SANGER/_data"
file = sprintf("%s/CNV_Category.csv", dir)
write.csv(CNV_TGSA_SG, file=file, row.names=T)


### DRPreter
dir = "../../benchmark_test/DRPreter/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_DRPreter = fread_def(file)

genes_drpreter = EXP_DRPreter %>% colnames
EXP_DRPreter_TCGA = TCGA_RNA_Sym %>% unify_genes(genes_drpreter)   # 2369 [NA 3]

dir = "../../benchmark_test/DRPreter/_data"
file = sprintf("%s/EXP_TCGA.csv", dir)
fwrite(EXP_DRPreter_TCGA, file=file, row.names=T, col.names=T)


idx1 = match(rownames(EXP_DRPreter), Anno_Cells$BROAD_ID)
idx2 = match(rownames(EXP_DRPreter_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("CCLE", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_DRPreter, EXP_DRPreter_TCGA)
EXP_DRPreter_TCGA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/DRPreter/_data"
file = sprintf("%s/EXP_TCGA_ComBat.csv", dir)
fwrite(EXP_DRPreter_TCGA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [DRPreter, %s]", dir_pca, category)
idx_tcga = rownames(EXP_DRPreter_TCGA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(EXP_DRPreter %>% scale_, EXP_DRPreter_TCGA %>% scale_)

PCA_DRPreter = plot_pca_tcga(Temp_List, EXP_DRPreter_TCGA_, 
                             db_list, tcga_list, main=main, 
                             idx_tcga=idx_tcga, train_db="CCLE")


### DRPreter_SANGER
dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_DRPreter_SG = fread_def(file)

genes_drpreter_sg = EXP_DRPreter_SG %>% colnames
EXP_DRPreter_SG_TCGA = TCGA_RNA_Sym %>% unify_genes(genes_drpreter_sg)   # 2369 [NA 3]

dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP_TCGA.csv", dir)
fwrite(EXP_DRPreter_SG_TCGA, file=file, row.names=T, col.names=T)


idx1 = match(rownames(EXP_DRPreter_SG), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(EXP_DRPreter_SG_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_DRPreter_SG, EXP_DRPreter_SG_TCGA)
EXP_DRPreter_SG_TCGA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP_TCGA_ComBat.csv", dir)
fwrite(EXP_DRPreter_SG_TCGA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [DRPreter_SANGER, %s]", dir_pca, category)
idx_tcga = rownames(EXP_DRPreter_SG_TCGA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(EXP_DRPreter_SG %>% scale_, EXP_DRPreter_SG_TCGA %>% scale_)

PCA_DRPreter_SG = plot_pca_tcga(Temp_List, EXP_DRPreter_SG_TCGA_, 
                                db_list, tcga_list, main=main, 
                                idx_tcga=idx_tcga, train_db="SANGER")


### RF
# robust_scaler <- function(train_df) {
#   medians <- apply(train_df, 2, median)
#   iqrs <- apply(train_df, 2, IQR)
#   transform <- function(new_df) {
#     scaled <- sweep(new_df, 2, medians, "-")
#     scaled <- sweep(scaled, 2, iqrs, "/")
#     return(as.data.frame(scaled))
#   }
#   return(transform)
# }

robust_scaler <- function(df) {
  medians <- apply(df, 2, median)
  iqrs <- apply(df, 2, IQR)
  scaled <- sweep(df, 2, medians, "-")
  scaled <- sweep(scaled, 2, iqrs, "/")
  scaled <- as.matrix(scaled)
  scaled[is.na(scaled)] <- 0
  scaled[is.infinite(scaled)] <- 0
  return(as.data.frame(scaled))
}

dir = "../../benchmark_test/RF/_data"
file = sprintf("%s/SANGER_RNA_TPM.csv", dir)
EXP_RF = fread_def(file, col_numeric=T)

genes_rf = EXP_RF %>% colnames
EXP_RF_TCGA = TCGA_RNA_Ent %>% unify_genes(genes_rf)   # 1001 [NA 1]

dir = "../../benchmark_test/RF/_data"
file = sprintf("%s/EXP_TCGA.csv", dir)
fwrite(EXP_RF_TCGA, file=file, row.names=T, col.names=T)


idx1 = match(rownames(EXP_RF), Anno_Cells$SANGER_MODEL_ID)
idx2 = match(rownames(EXP_RF_TCGA), TCGA_Info$cases)
tissue_ori = Anno_Cells$TCGA_CODE[idx1]
tissue_tcga = gsub("TCGA-", "", TCGA_Info$project[idx2])
tissue_ori[is.na(tissue_ori)] = "UNCLASSIFIED"
tissue_tcga[is.na(tissue_tcga)] = "UNCLASSIFIED"

db_list = c("SANGER", "TCGA")
tcga_list = list(tissue_ori, tissue_tcga)
exp_list = list(EXP_RF, EXP_RF_TCGA)
EXP_RF_TCGA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/RF/_data"
file = sprintf("%s/EXP_TCGA_ComBat.csv", dir)
fwrite(EXP_RF_TCGA_[[2]], file=file, row.names=T, col.names=T)


main = sprintf("%s/PCA Plot [RF, %s]", dir_pca, category)
idx_tcga = rownames(EXP_RF_TCGA_[[2]]) %in% TCGA_Resp$Sample_RNA
Temp_List = list(EXP_RF %>% robust_scaler, EXP_RF_TCGA %>% robust_scaler)

PCA_RF = plot_pca_tcga(Temp_List, EXP_RF_TCGA_, 
                       db_list, tcga_list, main=main, 
                       idx_tcga=idx_tcga, train_db="SANGER")


# library(caret)
# library(dplyr)
# library(ggplot2)
# 
# set.seed(42)
# a <- data.frame(matrix(rnorm(100 * 5, mean = 5), ncol = 5))
# b <- data.frame(matrix(rnorm(50 * 5, mean = 5.5), ncol = 5))
# colnames(a) <- colnames(b) <- paste0("V", 1:5)
# 
# # 1. StandardScaler
# # preproc <- caret::preProcess(a, method = c("center", "scale"))
# # a_scaled <- predict(preproc, a)
# # b_scaled <- predict(preproc, b)
# 
# # 2. RobustScaler
# scaler <- robust_scaler(a)
# a_scaled <- scaler(a)
# b_scaled <- scaler(b)
# 
# ab <- rbind(a, b)
# ab_scaled <- rbind(a_scaled, b_scaled)
# pca_model1 <- prcomp(ab, center = T, scale. = T)
# pca_model2 <- prcomp(ab_scaled, center = T, scale. = T)
# 
# group = c(rep("A", nrow(a)), rep("B", nrow(b)))
# pca_model1 = pca_model1$x %>% as.data.frame %>% mutate(Group=group)
# pca_model2 = pca_model2$x %>% as.data.frame %>% mutate(Group=group)
# 
# pca_model1 %>% ggplot(aes(PC1, PC2, color=group)) + geom_point() + ggtitle("PC2 Plot-1")
# pca_model2 %>% ggplot(aes(PC1, PC2, color=group)) + geom_point() + ggtitle("PC2 Plot-2")


# ### RF
# dir = "../../benchmark_test/RF/_data"
# file = sprintf("%s/MUT.csv", dir)
# MUT_RF = read.csv(file, row.names=1, check.names=F)
# 
# genes_rf = MUT_RF %>% colnames
# MUT_RF_TCGA = data_mut %>% subset(Hugo_Symbol %in% genes_rf)   # 1443 x 140
# MUT_RF_TCGA = MUT_RF_TCGA %>% as.data.frame %>% 
#   reshape2::acast(Tumor_Sample_Barcode~Hugo_Symbol, 
#                   value.var="Hugo_Symbol", fun.aggregate=any, fill=0)   # 340 x 136
# 
# MUT_RF_TCGA = MUT_RF_TCGA %>% unify_genes(genes_rf)   # [NA] 9
# MUT_RF_TCGA[is.na(MUT_RF_TCGA)] = 0
# 
# all(rownames(MUT_RF_TCGA) %in% TCGA_Resp$Sample)       # F
# all(rownames(MUT_RF_TCGA) %in% TCGA_Resp_Ori$Sample)   # F
# samples_rf = gsub("^(([^-]+-){3}[^-]+).*", "\\1", rownames(MUT_RF_TCGA))
# 
# idx = samples_rf %in% TCGA_Resp$Sample
# MUT_RF_TCGA = MUT_RF_TCGA[idx, ]   # 338 x 145
# 
# samples_rf = gsub("^(([^-]+-){3}[^-]+).*", "\\1", rownames(MUT_RF_TCGA))
# samples_rf %>% unique %>% length   # 338 [338]
# 
# Anno_RF_TCGA = data.frame(Sample=samples_rf, Sample_MUT=rownames(MUT_RF_TCGA))
# rownames(MUT_RF_TCGA) = samples_rf
# 
# dir = "../../benchmark_test/RF/_data"
# file = sprintf("%s/MUT_TCGA.csv", dir)
# fwrite(MUT_RF_TCGA, file=file, row.names=T, col.names=T)
# 
# file = sprintf("%s/TCGA_Response_.csv", dir)
# fwrite(TCGA_Resp, file=file, col.names=T)


if (supplementary) {
  save_for_nc = function(df_list, dir=".", num=1, num_fig=NULL, rowNames=F, suppl=T) {
    suppressMessages(library(openxlsx))
    is_list = inherits(df_list, "list")
    if (is_list & is.null(num_fig)) num_fig = letters[1:length(df_list)]
    
    if (!suppl) {
      sheets = sprintf("Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Fig. %s", num)
      file = sprintf("%s/SourceData_Fig%s.xlsx", dir, num)
    } else {
      sheets = sprintf("Supplementary Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Supplementary Fig. %s", num)
      file = sprintf("%s/SourceData_SupplementaryFig%s.xlsx", dir, num)
    }
    
    write.xlsx(df_list, file=file, sheetName=sheets, rowNames=rowNames)
  }
  
  
  ### Supplementary Data 9
  col = c("Drug_Name", "Drug_CID", "GDSC", "Inchi_Key", "SMILES_Can")
  TCGA_Drug_Info_ = TCGA_Drug_Info[, col]
  colnames(TCGA_Drug_Info_)[3:5] = c("GDSC_Drug", "InCHI_Key", "Canonical_SIMLES")
  
  TCGA_Resp_Num = TCGA_Resp %>% 
    mutate(Drug_CID=as.character(Drug_CID)) %>% group_by(Drug_CID) %>% 
    summarise(Num_CR=sum(Response=="Complete Response"), 
              Num_PR=sum(Response=="Partial Response"), 
              Num_SD=sum(Response=="Stable Disease"), 
              Num_PD=sum(Response=="Clinical Progressive Disease")) %>% as.data.frame
  
  TCGA_Drug_Info_ = TCGA_Drug_Info_ %>% 
    left_join(TCGA_Resp_Num, by="Drug_CID") %>% 
    relocate(Num_CR, Num_PR, Num_SD, Num_PD, .after=GDSC_Drug) %>% as.data.frame
  
  sheet = "Supplementary Data 9"
  file = sprintf("%s.xlsx", sheet)
  write.xlsx(TCGA_Drug_Info_, file=file, rowNames=F)
  
  
  ### Supplementary Data 10
  TCGA_Resp_ = TCGA_Resp %>% 
    rename(Cancer_Type=TCGA_Code) %>% 
    relocate(Cancer_Type, .after=Sample_Type)
  
  sheet = "Supplementary Data 10"
  file = sprintf("%s.xlsx", sheet)
  write.xlsx(TCGA_Resp_, file=file, sheetName=sheet, rowNames=F)
  
  
  ### [Source Data] Supplementary Fig. 34
  PCA_ = list(PCA_Genes[[1]] %>% relocate(TCGA_Code, .after=Database), 
              PCA_GCNPath[[1]] %>% relocate(TCGA_Code, .after=Database), 
              PCA_GCNPath[[2]] %>% relocate(TCGA_Code, .after=Database))
  
  PCA_ %>% save_for_nc(num=34, suppl=T, rowNames=T)
}

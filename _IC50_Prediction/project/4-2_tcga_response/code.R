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
SANGER_TCGA_GSVA_PCA %>% subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>% 
  plot_def(PC1, PC2, color=TCGA_Code, main=file, 
           xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
           alpha=0.8, size=2, width=20, height=16, force_bold=T, save=T)

file = sprintf("%s/PC2 Plot [GSVA, TCGA, Top 5]", dir)
SANGER_TCGA_GSVA_PCA %>% subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>% 
  plot_def(PC1, PC2, color=TCGA_Code, main=file, 
           xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
           alpha=0.8, size=1.2, width=20, height=16, force_bold=T, save=T)


file = sprintf("%s/PCA_GSVA.csv", dir)
write.csv(SANGER_TCGA_GSVA_PCA, file=file, row.names=T)


# ### 2. PCA Plot [GSVA, Reference ComBat]
# 
# identical(colnames(SANGER_RNA_GSVA), colnames(TCGA_RNA_GSVA))   # T
# SANGER_TCGA_GSVA = rbind(SANGER_RNA_GSVA, TCGA_RNA_GSVA)
# SANGER_TCGA_GSVA_CB1 = ComBat(dat=t(SANGER_TCGA_GSVA), batch=db, mod=NULL, ref.batch="SANGER")
# 
# idx = which(db=="SANGER")
# SANGER_TCGA_GSVA_CB1 = SANGER_TCGA_GSVA_CB1 %>% t %>% as.data.frame
# identical(SANGER_RNA_GSVA, SANGER_TCGA_GSVA_CB1[idx, ])   # T
# 
# SANGER_TCGA_GSVA_PCA_CB1 = SANGER_TCGA_GSVA_CB1 %>%
#   pca_exp(db=db, tcga_code=tcga_code)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% arrange(Database) %>%
#   plot_def(PC1, PC2, color=Database, main=file,
#            xlab="PC1", ylab="PC2", alpha=0.8, size=0.8,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref, SANGER, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref, TCGA, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2,
#            width=20, height=16, force_bold=T, save=T)
# 
# 
# ### 3. PCA Plot [GSVA, Reference ComBat (TCGA tissue as Co-variate)]
# 
# mod = as.data.frame(tcga_code)
# mod = model.matrix(~as.factor(tcga_code), data=mod)
# SANGER_TCGA_GSVA_CB2 = ComBat(dat=t(SANGER_TCGA_GSVA), batch=db, mod=mod, ref.batch="SANGER")
# 
# idx = which(db=="SANGER")
# SANGER_TCGA_GSVA_CB2 = SANGER_TCGA_GSVA_CB2 %>% t %>% as.data.frame
# identical(SANGER_RNA_GSVA, SANGER_TCGA_GSVA_CB2[idx, ])   # T
# 
# SANGER_TCGA_GSVA_PCA_CB2 = SANGER_TCGA_GSVA_CB2 %>%
#   pca_exp(db=db, tcga_code=tcga_code)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% arrange(Database) %>%
#   plot_def(PC1, PC2, color=Database, main=file,
#            xlab="PC1", ylab="PC2", alpha=0.8, size=0.8,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA, SANGER, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA, TCGA, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2,
#            width=20, height=16, force_bold=T, save=T)



# ### 2. PCA Plot [GSVA, Reference ComBat]
# 
# fill_col_na = function(df, col, value=0) {
#   
#   if (!is.data.frame(df)) df = df %>% as.data.frame
#   col_na = col[!(col %in% colnames(df))]
#   
#   df_col_na = rep(value, nrow(df)*length(col_na)) %>% 
#     matrix(nrow=nrow(df), ncol=length(col_na)) %>% as.data.frame
#   rownames(df_col_na) = rownames(df)
#   colnames(df_col_na) = col_na
#   
#   df = df %>% cbind(df_col_na)
#   df = df[, col]
#   return(df)
# }
# 
# rbind_union = function(df1, df2, value=0) {
#   col = union(colnames(df1), colnames(df2))
#   df1 = df1 %>% fill_col_na(col, value)
#   df2 = df2 %>% fill_col_na(col, value)
#   df = rbind(df1, df2)
#   return(df)
# }
# 
# suppressMessages(library(sva))
# # col = intersect(colnames(SANGER_RNA_Ent), colnames(TCGA_RNA_Ent))     # 37293
# # SANGER_TCGA_RNA = rbind(SANGER_RNA_Ent[, col], TCGA_RNA_Ent[, col])   # 12705 x 37293
# SANGER_TCGA_RNA = rbind_union(SANGER_RNA_Ent, TCGA_RNA_Ent)             # 12705 x 40902
# 
# # Remove columns with SD=0 to implement Combat
# # https://stackoverflow.com/questions/21532998/error-when-using-combat
# col_sd = sapply(SANGER_TCGA_RNA, sd)
# SANGER_TCGA_RNA = SANGER_TCGA_RNA[, col_sd!=0]   # 40902 > 40293
# 
# SANGER_TCGA_RNA_CB1 = ComBat(dat=t(SANGER_TCGA_RNA), batch=db, mod=NULL, ref.batch="SANGER")
# # Found 4231 genes with uniform expression within a single batch (all zeros)
# # these will not be adjusted for batch
# 
# idx = which(db=="SANGER")
# SANGER_TCGA_RNA_CB1 = SANGER_TCGA_RNA_CB1 %>% t %>% as.data.frame
# 
# identical(rownames(SANGER_TCGA_RNA), rownames(SANGER_TCGA_RNA_CB1))
# identical(colnames(SANGER_TCGA_RNA), colnames(SANGER_TCGA_RNA_CB1))
# identical(SANGER_TCGA_RNA[idx, ], SANGER_TCGA_RNA_CB1[idx, ])   # T
# TCGA_RNA_GSVA_CB1 = SANGER_TCGA_RNA_CB1[-idx, ] %>% gsva_def(Path_List, method="gsva", cores=cores)
# 
# identical(colnames(SANGER_RNA_GSVA), colnames(TCGA_RNA_GSVA_CB1))   # T
# SANGER_TCGA_GSVA_PCA_CB1 = rbind(SANGER_RNA_GSVA, TCGA_RNA_GSVA_CB1) %>% 
#   pca_exp(db=db, tcga_code=tcga_code)
# 
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% arrange(Database) %>% 
#   plot_def(PC1, PC2, color=Database, main=file, 
#            xlab="PC1", ylab="PC2", alpha=0.8, size=0.8, 
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref, SANGER, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>% 
#   plot_def(PC1, PC2, color=TCGA_Code, main=file, 
#            xlim=c(-18, 18), ylim=c(-15, 20), 
#            xlab="PC1", ylab="PC2", alpha=0.8, size=2, 
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref, TCGA, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB1 %>% subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>% 
#   plot_def(PC1, PC2, color=TCGA_Code, main=file, 
#            xlim=c(-18, 18), ylim=c(-15, 20), 
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2, 
#            width=20, height=16, force_bold=T, save=T)
# 
# 
# ### 3. PCA Plot [GSVA, Reference ComBat (TCGA tissue as Co-variate)]
# 
# mod = as.data.frame(tcga_code)
# mod = model.matrix(~as.factor(tcga_code), data=mod)
# 
# SANGER_TCGA_RNA_CB2 = ComBat(dat=t(SANGER_TCGA_RNA), batch=db, mod=mod, ref.batch="SANGER")
# # Found 4231 genes with uniform expression within a single batch (all zeros)
# # these will not be adjusted for batch
# 
# idx = which(db=="SANGER")
# SANGER_TCGA_RNA_CB2 = SANGER_TCGA_RNA_CB2 %>% t %>% as.data.frame
# 
# identical(rownames(SANGER_TCGA_RNA), rownames(SANGER_TCGA_RNA_CB2))
# identical(colnames(SANGER_TCGA_RNA), colnames(SANGER_TCGA_RNA_CB2))
# identical(SANGER_TCGA_RNA[idx, ], SANGER_TCGA_RNA_CB2[idx, ])   # T
# TCGA_RNA_GSVA_CB2 = SANGER_TCGA_RNA_CB2[-idx, ] %>% gsva_def(Path_List, method="gsva", cores=cores)
# 
# identical(colnames(SANGER_RNA_GSVA), colnames(TCGA_RNA_GSVA_CB2))   # T
# SANGER_TCGA_GSVA_PCA_CB2 = rbind(SANGER_RNA_GSVA, TCGA_RNA_GSVA_CB2) %>% 
#   pca_exp(db=db, tcga_code=tcga_code)
# 
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% arrange(Database) %>%
#   plot_def(PC1, PC2, color=Database, main=file,
#            xlab="PC1", ylab="PC2", alpha=0.8, size=0.8,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA, SANGER, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% subset(Database=="SANGER" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=2,
#            width=20, height=16, force_bold=T, save=T)
# 
# file = sprintf("%s/PC2 Plot [GSVA + Combat_Ref_TCGA, TCGA, Top 5]", dir)
# SANGER_TCGA_GSVA_PCA_CB2 %>% subset(Database=="TCGA" & TCGA_Code %in% tcga_main) %>%
#   plot_def(PC1, PC2, color=TCGA_Code, main=file,
#            xlim=c(-18, 18), ylim=c(-15, 20),
#            xlab="PC1", ylab="PC2", alpha=0.8, size=1.2,
#            width=20, height=16, force_bold=T, save=T)


### 4. PCA Plot [GSVA]

# dir = "../../processed_data/BIOCARTA/GSVA_Analysis/SANGER & TCGA"

# file1 = sprintf("%s/SANGER_TCGA_RNA.csv", dir)
# file2 = sprintf("%s/SANGER_TCGA_RNA_CB1.csv", dir)
# file3 = sprintf("%s/SANGER_TCGA_RNA_CB2.csv", dir)
# file4 = sprintf("%s/SANGER_TCGA_GSVA.csv", dir)

# fwrite(SANGER_TCGA_RNA, file=file1, row.names=T)
# fwrite(SANGER_TCGA_RNA_CB1, file=file2, row.names=T)
# fwrite(SANGER_TCGA_RNA_CB2, file=file3, row.names=T)
# fwrite(SANGER_TCGA_GSVA, file=file4, row.names=T)



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


if (supplementary) {
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
}

#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))
suppressMessages(library(cogena))
suppressMessages(library(impute))

# source("../1-3_gsva_pcn/functions.R")
source("../calc_pathway_scores.R")
source("../functions.R")
source("functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

cancer_type = "Small Cell Lung Carcinoma"
cells_sclc = Anno_Cells %>% 
  subset(CANCER_TYPE==cancer_type) %>% 
  pull(SANGER_MODEL_ID) %>% unique   # 81

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.symbols.gmt", dir)
Path_List = gmt2list(file)
geneset = Path_List %>% unlist %>% unique   # 1509

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Sym.csv", dir)
SANGER_RNA_TPM = fread_def(file, check_names=T)

dir = "../../processed_data/cell_data/BIOCARTA"
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
SANGER_RNA_GSVA = read.csv(file, row.names=1, check.names=F)


# file = "Table S1A_Clinical.csv"
# Clinical = fread(file)

# file = "1-s2.0-S0092867423013351-mmc6.xlsx"
# Clinical_ = read.xlsx(file, sheet=2)

file = "Table S1B_MUT.csv"
MUT = fread(file)
samples_t = MUT$Tumor_Sample_Barcode %>% unique          # 112
samples_n = MUT$Matched_Norm_Sample_Barcode %>% unique   # 112

file = "Table S1C_CNV.csv"
CNV = fread_def(file)
CNV %>% is.na %>% sum   # 0
CNV = CNV %>% t %>% as.data.frame
# CNV %>% unlist %>% hist   # Median 0.002

file = "Table S1D_RNA.csv"
RNA = fread_def(file)
RNA[is.na(RNA)] = 0
RNA = RNA %>% t %>% as.data.frame

file = "Table S1E_Prot.csv"
Prot = fread_def(file)
Prot = Prot %>% as.matrix %>% impute.knn
Prot = Prot$data %>% t %>% as.data.frame

category = c("ComBat X", "ComBat O")
dir_pca = mkdir("Batch Effects with PCA Plot")





##### 2. Process Data [Omics]

# Raw data from BIOCARTA genes [TPM]
cells = rownames(SANGER_RNA_TPM) %in% cells_sclc
geneset_ = Reduce(intersect, list(geneset, colnames(RNA), colnames(SANGER_RNA_TPM)))
SANGER_RNA_ = SANGER_RNA_TPM[cells, geneset_]   # 72 x 1443
RNA_ = RNA[, geneset_]   # 214 x 1443

anno_src = rep("Tumor", nrow(SANGER_RNA_))
anno_tgt = ifelse(grepl("T", rownames(RNA_)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(SANGER_RNA_, RNA_)
RNA_ = combat_def(exp_list, tcga_list, db_list)

main = sprintf("%s/PCA Plot [BIOCARTA, TPM & TPM]", dir_pca)
PCA_Genes = plot_pca_tcga(exp_list, RNA_, db_list, main=main, width=21)


# Raw data from BIOCARTA genes [Proteome]
cells = rownames(SANGER_RNA_TPM) %in% cells_sclc
geneset_ = Reduce(intersect, list(geneset, colnames(Prot), colnames(SANGER_RNA_TPM)))
SANGER_RNA_ = SANGER_RNA_TPM[cells, geneset_]   # 72 x 1014
Prot_ = Prot[, geneset_]   # 214 x 1014

anno_src = rep("Tumor", nrow(SANGER_RNA_))
anno_tgt = ifelse(grepl("T", rownames(Prot_)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(SANGER_RNA_, Prot_)
Prot_ = combat_def(exp_list, tcga_list, db_list)

main = sprintf("%s/PCA Plot [BIOCARTA, TPM & Proteom]", dir_pca)
PCA_Genes_Prot = plot_pca_tcga(exp_list, Prot_, db_list, main=main, width=21)


### GCNPath
# Calculate GSVA pathway scores
cores = 24
sum(colnames(RNA) %in% geneset)    # 1452/1509 [96.22%]
RNA_GSVA = RNA %>% calc_path_def(Path_List, method="gsva", cores=cores)

dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
file = sprintf("%s/Liu24_invivo_RNA_GSVA.csv", dir)
fwrite(RNA_GSVA, file=file, row.names=T, col.names=T)

sum(colnames(Prot) %in% geneset)   # 1019/1509 [67.52%]
Prot_GSVA = Prot %>% calc_path_def(Path_List, method="gsva", cores=cores)

dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
file = sprintf("%s/Liu24_invivo_Prot_GSVA.csv", dir)
fwrite(Prot_GSVA, file=file, row.names=T, col.names=T)


cells = rownames(SANGER_RNA_GSVA) %in% cells_sclc
SANGER_RNA_GSVA = SANGER_RNA_GSVA[cells, ]   # 72 x 292

anno_src = rep("Tumor", nrow(SANGER_RNA_GSVA))
anno_tgt = ifelse(grepl("T", rownames(RNA_GSVA)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(SANGER_RNA_GSVA, RNA_GSVA)
RNA_GSVA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
file = sprintf("%s/Liu24_invivo_RNA_GSVA_ComBat.csv", dir)
fwrite(RNA_GSVA_[[2]], file=file, row.names=T, col.names=T)

main = sprintf("%s/PCA Plot [GCNPath, %s]", dir_pca, category)
PCA_GCNPath = plot_pca_tcga(exp_list, RNA_GSVA_, db_list, main=main, width=21)


anno_src = rep("Tumor", nrow(SANGER_RNA_GSVA))
anno_tgt = ifelse(grepl("T", rownames(Prot_GSVA)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(SANGER_RNA_GSVA, Prot_GSVA)
Prot_GSVA_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/GCNPath/processed/cell_data_biocarta"
file = sprintf("%s/Liu24_invivo_Prot_GSVA_ComBat.csv", dir)
fwrite(Prot_GSVA_[[2]], file=file, row.names=T, col.names=T)

main = sprintf("%s/PCA Plot [GCNPath (Test with Proteome), %s]", dir_pca, category)
PCA_GCNPath_Prot = plot_pca_tcga(exp_list, Prot_GSVA_, db_list, main=main, width=21)


### DRPreter_SANGER
dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_DRPreter_SG = fread_def(file)

genes_drpreter_sg = EXP_DRPreter_SG %>% colnames
EXP_DRPreter_SG_Ext = RNA %>% unify_genes(genes_drpreter_sg)   # 2369 [NA 138]

dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP_Liu_2024.csv", dir)
fwrite(EXP_DRPreter_SG_Ext, file=file, row.names=T, col.names=T)


cells = rownames(EXP_DRPreter_SG) %in% cells_sclc
EXP_DRPreter_SG = EXP_DRPreter_SG[cells, ]   # 72 x 2369

anno_src = rep("Tumor", nrow(EXP_DRPreter_SG))
anno_tgt = ifelse(grepl("T", rownames(EXP_DRPreter_SG_Ext)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(EXP_DRPreter_SG, EXP_DRPreter_SG_Ext)
EXP_DRPreter_SG_Ext_ = combat_def(exp_list, tcga_list, db_list)

dir = "../../benchmark_test/DRPreter_SANGER/_data"
file = sprintf("%s/EXP_Liu_2024_ComBat.csv", dir)
fwrite(EXP_DRPreter_SG_Ext_[[2]], file=file, row.names=T, col.names=T)

main = sprintf("%s/PCA Plot [DRPreter_SANGER, %s]", dir_pca, category)
PCA_DRPreter_SG = plot_pca_tcga(exp_list, EXP_DRPreter_SG_Ext_, db_list, main=main, width=20)


### TGDRP_SANGER
dir = "../../benchmark_test/TGSA_SANGER/_data"
file = sprintf("%s/EXP.csv", dir)
EXP_TGSA_SG = read.csv(file, row.names=1, check.names=F)
genes_tgsa_sg = EXP_TGSA_SG %>% colnames

idx = match(genes_tgsa_sg, Anno_Genes$ENTREZ_ID)
genes_tgsa_sg_ = Anno_Genes$HGNC_SYMBOL[idx]
genes_tgsa_sg_ %>% is.na %>% sum   # 0

genes_tgsa_sg = EXP_TGSA_SG %>% colnames
EXP_TGSA_SG_Ext = RNA %>% unify_genes(genes_tgsa_sg_)   # 2369 [NA 27]
colnames(EXP_TGSA_SG_Ext) = genes_tgsa_sg

dir = mkdir("../../benchmark_test/TGSA_SANGER/_data_liu24_invivo")
file = sprintf("%s/EXP.csv", dir)
fwrite(EXP_TGSA_SG_Ext, file=file, row.names=T, col.names=T)


cells = rownames(EXP_TGSA_SG) %in% cells_sclc
EXP_TGSA_SG = EXP_TGSA_SG[cells, ]   # 72 x 706

anno_src = rep("Tumor", nrow(EXP_TGSA_SG))
anno_tgt = ifelse(grepl("T", rownames(EXP_TGSA_SG_Ext)), "Tumor", "NAT")
tcga_list = list(anno_src, anno_tgt)

db_list = c("SANGER", "Liu et al. (2024)")
exp_list = list(EXP_TGSA_SG, EXP_TGSA_SG_Ext)
EXP_TGSA_SG_Ext_ = combat_def(exp_list, tcga_list, db_list)

dir = mkdir("../../benchmark_test/TGSA_SANGER/_data_liu24_invivo")
file = sprintf("%s/EXP_ComBat.csv", dir)
fwrite(EXP_TGSA_SG_Ext_[[2]], file=file, row.names=T, col.names=T)

main = sprintf("%s/PCA Plot [TGSA_SANGER, %s]", dir_pca, category)
PCA_TGSA_SG = plot_pca_tcga(exp_list, EXP_TGSA_SG_Ext_, db_list, main=main, width=20)


# Mutation
MUT_TGSA_SG_Ext = MUT %>% subset(Hugo_Symbol %in% genes_tgsa_sg_)   # 182869 > 8621
MUT_TGSA_SG_Ext = MUT_TGSA_SG_Ext %>% 
  reshape2::acast(Tumor_Sample_Barcode~Hugo_Symbol, value.var="Hugo_Symbol",
                  fun.aggregate = function(x) as.integer(length(x)>0), fill=0)

MUT_TGSA_SG_Ext = MUT_TGSA_SG_Ext %>% unify_genes(genes_tgsa_sg_)   # [NA] 23
colnames(MUT_TGSA_SG_Ext) = genes_tgsa_sg

MUT_TGSA_SG_Ext_N = matrix(rep(0, length(samples_n)*ncol(MUT_TGSA_SG_Ext)), nrow=length(samples_n))
MUT_TGSA_SG_Ext_N = MUT_TGSA_SG_Ext_N %>% data.frame(row.names=samples_n) %>% setNames(colnames(MUT_TGSA_SG_Ext))
MUT_TGSA_SG_Ext = MUT_TGSA_SG_Ext %>% rbind(MUT_TGSA_SG_Ext_N)   # 224 x 706

all(rownames(EXP_TGSA_SG_Ext) %in% rownames(MUT_TGSA_SG_Ext))    # T
MUT_TGSA_SG_Ext = MUT_TGSA_SG_Ext[rownames(EXP_TGSA_SG_Ext), ]   # 214 x 706

dir = mkdir("../../benchmark_test/TGSA_SANGER/_data_liu24_invivo")
file = sprintf("%s/MUT.csv", dir)
fwrite(MUT_TGSA_SG_Ext, file=file, row.names=T, col.names=T)


# CNV
CNV_TGSA_SG_Ext = CNV %>% unify_genes(genes_tgsa_sg_)   # [NA] 6
colnames(CNV_TGSA_SG_Ext) = genes_tgsa_sg

CNV_TGSA_SG_Ext_N = matrix(rep(0, length(samples_n)*ncol(CNV_TGSA_SG_Ext)), nrow=length(samples_n))
CNV_TGSA_SG_Ext_N = CNV_TGSA_SG_Ext_N %>% data.frame(row.names=samples_n) %>% setNames(colnames(CNV_TGSA_SG_Ext))
CNV_TGSA_SG_Ext = CNV_TGSA_SG_Ext %>% rbind(CNV_TGSA_SG_Ext_N)   # 224 x 706

all(rownames(EXP_TGSA_SG_Ext) %in% rownames(CNV_TGSA_SG_Ext))    # T
CNV_TGSA_SG_Ext = CNV_TGSA_SG_Ext[rownames(EXP_TGSA_SG_Ext), ]   # 214 x 706

# Rescale log2 CN-Ratio into log2(x+1)
CNV_TGSA_SG_Ext %>% unlist %>% hist
CNV_TGSA_SG_Ext %>% unlist %>% mean %>% round(3)   # 0.022

CNV_TGSA_SG_Ext = log2(2**CNV_TGSA_SG_Ext+1)
CNV_TGSA_SG_Ext %>% unlist %>% hist
CNV_TGSA_SG_Ext %>% unlist %>% mean %>% round(3)   # 1.022

dir = mkdir("../../benchmark_test/TGSA_SANGER/_data_liu24_invivo")
file = sprintf("%s/CNV.csv", dir)
fwrite(CNV_TGSA_SG_Ext, file=file, row.names=T, col.names=T)




##### 3. Process Data [Drug Response]

samples = rownames(RNA_GSVA)
drugs = c("Etoposide", "Cisplatin", "Anlotinib", "Alisertib", "Barasertib", "AMG-900")
cids = c(36462, 5702198, 25017411, 24771867, 11497983, 24856041)

Pred = expand.grid(Sample=samples, Drug=drugs)
idx = match(Pred$Drug, drugs)
Pred$Drug_CID = cids[idx]

Pred$Status = ifelse(grepl("^T", Pred$Sample), "Tumor", "Normal")
Pred = Pred %>% relocate(Status, .after=Sample)

file = "Liu_2024_in_vivo.txt"
write.table(Pred, file=file, sep="\t", row.names=F, col.names=T, quote=F)

# Get SMILES using CIDs from PubChem Homepage
# https://pubchem.ncbi.nlm.nih.gov
# SMILES_Liu_2024.csv




supplementary = T
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
  

  ### [Source Data] Supplementary Fig. 36
  PCA_Genes_ = rbind(PCA_Genes[[1]] %>% mutate(Omics_Type="Log2(TPM+1)"), 
                     PCA_Genes_Prot[[1]] %>% mutate(Omics_Type="Proteome"))
  
  PCA_GSVA_ = rbind(PCA_GCNPath[[1]] %>% mutate(Omics_Type="Log2(TPM+1)"), 
                    PCA_GCNPath_Prot[[1]] %>% mutate(Omics_Type="Proteome"))
  
  PCA_GSVA_CB_ = rbind(PCA_GCNPath[[2]] %>% mutate(Omics_Type="Log2(TPM+1)"), 
                       PCA_GCNPath_Prot[[2]] %>% mutate(Omics_Type="Proteome"))
  
  Temp = list(PCA_Genes_ %>% relocate(Omics_Type, .after=Database), 
              PCA_GSVA_ %>% relocate(Omics_Type, .after=Database), 
              PCA_GSVA_CB_ %>% relocate(Omics_Type, .after=Database))
  
  Temp %>% save_for_nc(num=36, suppl=T, rowNames=T)
}

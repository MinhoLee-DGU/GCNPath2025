#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(cogena))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))

source("../utils/functions.R")
source("functions.R")
loadings()


dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC.txt", dir)
IC50_GDSC = read.csv(file, sep="\t")

dir = "../processed/cell_data_biocarta"
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
SANGER_RNA_GSVA = read.csv(file, row.names=1)

dir = "../data/path_data"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.symbols.gmt", dir)
Path_List = gmt2list(file)
geneset = unique(unlist(Path_List))

file = "TPM_Sym.csv"
SANGER_RNA_TPM = fread_def(file)
# SANGER_RNA_TPM = SANGER_RNA_TPM[, colnames(SANGER_RNA_TPM) %in% geneset]

file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)
Anno_Cells$TCGA_CODE[is.na(Anno_Cells$TCGA_CODE)] = "UNCLASSIFIED"

file = "Anno_Drugs.csv"
Anno_Drugs = read.csv(file)
Anno_Drugs$Drug_CID = Anno_Drugs$Drug_CID %>% as.character
Anno_Drugs$Target_Pathway[is.na(Anno_Drugs$Target_Pathway)] = "Unclassified"


IC50_GDSC_Miss = IC50_GDSC %>% complete(Cell, Drug) %>% 
  subset(is.na(LN_IC50)) %>% as.data.frame

idx = match(IC50_GDSC_Miss$Cell, IC50_GDSC$Cell)

IC50_GDSC_Miss = IC50_GDSC_Miss %>% 
  subset(select=-LN_IC50) %>% 
  mutate(Cell_BROAD=IC50_GDSC$Cell_BROAD[idx], 
         Cell_COSMIC=IC50_GDSC$Cell_COSMIC[idx])

dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC_Miss.txt", dir)
write.table(IC50_GDSC_Miss, file=file, row.names=F, sep="\t")



##### 2-1. Prediction of missing IC50 values (Interporation)

pattern = "pred_total_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Test = read_pred(dir, pattern=pattern)   # 3717510
Pred_Test$Drug = Pred_Test$Drug %>% as.character

pattern = "pred_miss_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Miss = read_pred(dir, pattern=pattern)   # 481530
Pred_Miss$Drug = Pred_Miss$Drug %>% as.character
Pred_Miss = Pred_Miss %>% mutate(LN_IC50=NA) %>% relocate(LN_IC50, .before=Prediction)

# GCNPath shows outstanding and stable performances...
Perf_Test = Pred_Test %>% group_by(Seed) %>% 
  summarise(RMSE=RMSE(LN_IC50, Prediction), 
            PCC=cor(LN_IC50, Prediction), 
            SCC=cor(LN_IC50, Prediction, method="spearman")) %>% as.data.frame

Perf_Test$RMSE %>% mean %>% round(3)   # 0.691
Perf_Test$PCC %>% mean %>% round(3)    # 0.967
Perf_Test$SCC %>% mean %>% round(3)    # 0.952

col = c("Cell", "Drug", "LN_IC50", "Prediction")
Pred_Test = Pred_Test[, col, with=F] %>% as.data.frame %>% mutate(Missing=F)
Pred_Miss = Pred_Miss[, col, with=F] %>% as.data.frame %>% mutate(Missing=T)

Pred_Test = Pred_Test %>% group_by(Cell, Drug) %>% 
  mutate(Prediction=mean(Prediction)) %>% 
  subset(select=c(Cell, Drug, LN_IC50, Prediction, Missing)) %>% distinct %>% as.data.frame

Pred_Miss = Pred_Miss %>% group_by(Cell, Drug) %>% 
  summarise(LN_IC50=NA, Prediction=mean(Prediction), Missing=T) %>% as.data.frame

Pred_Test = rbind(Pred_Test, Pred_Miss)
Pred_Test$Drug = Pred_Test$Drug %>% as.character

Pred_Test %>% nrow                     # 419904 [972*432]
Pred_Test$Cell %>% unique %>% length   # 972
Pred_Test$Drug %>% unique %>% length   # 432

Pred_Test %>% subset(Missing) %>% nrow    # 48153
Pred_Test %>% subset(!Missing) %>% nrow   # 371751
all(Pred_Test$Missing[is.na(Pred_Test$LN_IC50)])     # T
all(!Pred_Test$Missing[!is.na(Pred_Test$LN_IC50)])   # T



##### 2-2. Prediction of cell-lines without IC50 values

# Run the code "test_rest.sh" to predict the Pred_Rest
cells = rownames(SANGER_RNA_GSVA) %>% setdiff(unique(Pred_Test$Cell))   # 459
drugs = unique(Pred_Test$Drug)   # 432

dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC_Rest.txt", dir)
Combn_Rest = expand.grid(Cell=cells, Drug=drugs) %>% as.data.frame
write.table(Combn_Rest, file=file, sep="\t", quote=F, row.names=F)

pattern = "pred_rest_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Rest = read_pred(dir, pattern=pattern)   # 1982880 [198288*10]
Pred_Rest$Drug = Pred_Rest$Drug %>% as.character

Pred_Rest = Pred_Rest %>% group_by(Cell, Drug) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 198288

Pred_Test$Rest = F
Pred_Rest = Pred_Rest %>% 
  mutate(LN_IC50=NA, Missing=F, Rest=T) %>% 
  relocate(LN_IC50, .before=Prediction) %>% as.data.frame

Pred_Test = Pred_Test %>% rbind(Pred_Rest)
Pred_Wide = Pred_Test %>% reshape2::acast(Cell~Drug, value.var="Prediction") %>% as.data.frame   # 1431 x 432


# Counting TCGA Code & Target Pathway [Cell & Drug]
Anno_Cells_SANGER = Anno_Cells %>% 
  subset(SANGER_MODEL_ID %in% rownames(SANGER_RNA_GSVA)) %>% 
  group_by(CANCER_TYPE_DETAIL) %>% summarise(Num=n()) %>% 
  arrange(desc(Num)) %>% as.data.frame

Anno_Drugs_SANGER = Anno_Drugs %>% 
  subset(Drug_CID %in% unique(IC50_GDSC$Drug)) %>% 
  group_by(Target_Pathway) %>% summarise(Num=n()) %>% 
  arrange(desc(Num)) %>% as.data.frame

# Complete cell annotation
# cells included in IC50_GDSC & SANGER_RNA?

IC50_GDSC_Num = IC50_GDSC %>% group_by(Cell) %>% 
  summarise(Num_IC50=n()) %>% as.data.frame
idx = match(Anno_Cells$SANGER_MODEL_ID, IC50_GDSC_Num$Cell)

Anno_Cells = Anno_Cells %>% 
  mutate(RNA=SANGER_MODEL_ID %in% rownames(SANGER_RNA_GSVA), 
         IC50=SANGER_MODEL_ID %in% unique(IC50_GDSC$Cell), 
         Num_IC50=IC50_GDSC_Num$Num_IC50[idx])

Anno_Cells$MODEL_NAME %>% is.na %>% sum   # 0



##### Examine well-characterized cancer types
# Anno_Cells$CANCER_TYPE %>% table %>% sort(decreasing=T) %>% head(10)

Pred_Test = get_anno_cd(Pred_Test, Anno_Cells, Anno_Drugs)
identical(rownames(SANGER_RNA_GSVA), rownames(SANGER_RNA_TPM))   # T


### COREAD
# https://www.cancer.gov/about-cancer/treatment/drugs/colorectal
# https://www.globalcca.org/learn/colorectal-cancer-targeted-therapy

tissue = "Large Intestine"
Pred_COREAD1 = list()
Pred_COREAD2 = list()

dir = mkdir("Case Study [Well-known, COREAD]")
drugs = c("5-Fluorouracil", "Oxaliplatin", "Lapatinib", "Trametinib")
targets = c("TP53", "TP53", "EGFR", "MAP2K1")
pathways = c("BIOCARTA_P53_PATHWAY", "BIOCARTA_P53_PATHWAY", "BIOCARTA_HER2_PATHWAY", "BIOCARTA_MAPK_PATHWAY")

for (i in 1:length(drugs)) {
  main1 = sprintf("%s/%s [%s]", dir, drugs[i], targets[i])
  main2 = sprintf("%s/%s [%s]", dir, drugs[i], pathways[i])
  
  sprintf("\n# %s [%s]", drugs[i], targets[i]) %>% cat
  Pred_COREAD1[[drugs[i]]] = Pred_Test %>% 
    subset(Tissue==tissue & Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_TPM, gene=targets[i], main=main1, save=T)
  
  sprintf("# %s [%s]", drugs[i], pathways[i]) %>% cat
  Pred_COREAD2[[drugs[i]]] = Pred_Test %>% 
    subset(Tissue==tissue & Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_GSVA, gene=pathways[i], main=main2, save=T)
}



### Breast
# https://www.cancer.gov/about-cancer/treatment/drugs/breast
# https://www.cancer.org/cancer/types/breast-cancer/treatment/targeted-therapy-for-breast-cancer.html

tissue = "Breast"
Pred_Breast1 = list()
Pred_Breast2 = list()

dir = mkdir("Case Study [Well-known, Breast]")
drugs = c("Alpelisib", "AZD5363", "Docetaxel", "Palbociclib")
targets = c("PIK3CA", "AKT1", "TUBA1C", "CDK6")
pathways = c("BIOCARTA_AKT_PATHWAY", "BIOCARTA_AKT_PATHWAY", "BIOCARTA_DEATH_PATHWAY", "BIOCARTA_RB_PATHWAY")

for (i in 1:length(drugs)) {
  main1 = sprintf("%s/%s [%s]", dir, drugs[i], targets[i])
  main2 = sprintf("%s/%s [%s]", dir, drugs[i], pathways[i])
  
  sprintf("\n# %s [%s]", drugs[i], targets[i]) %>% cat
  Pred_Breast1[[drugs[i]]] = Pred_Test %>% 
    subset(Tissue==tissue & Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_TPM, gene=targets[i], main=main1, save=T)
  
  sprintf("\n# %s [%s]", drugs[i], pathways[i]) %>% cat
  Pred_Breast2[[drugs[i]]] = Pred_Test %>% 
    subset(Tissue==tissue & Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_GSVA, gene=pathways[i], main=main2, save=T)
}


##### Examine well-known therapies in SCLC
# Etoposide [CID 36462]
# Cisplatin [CID 5702198]
# Lurbinectedin [Pol2, CID 57327016]

##### Examine promising therapies in SCLC
# Rovalpituzumab tesirine, Rova-T (target DLL3, pyrrolobenzodiazepine, CID 10030705)
# DS-7300 & HER3-DXd (target B7-H3 & HER3, deruxtecan, CID 166041952)
# ABBV-011 (target SEZ6, calicheamicin, CID 10953353)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11294301

# SCLC
tissue = "Lung"
Anno_Cells %>% lookup_by_tissue(tissue)

sclc = "Small Cell Lung Carcinoma"
cells_sclc = Anno_Cells %>% subset(CANCER_TYPE==sclc & RNA) %>% pull(SANGER_MODEL_ID)   # 72
drugs = c("Etoposide", "Cisplatin", "Lurbinectedin", "Pyrrolobenzodiazepine", "Deruxtecan", "Calicheamicin")
cids = c(36462, 5702198, 57327016, 10030705, 166041952, 10953353)

Pred_SCLC = expand.grid(Cell=cells_sclc, Drug=cids)
idx = match(Pred_SCLC$Drug, cids)
Pred_SCLC$Drug_Name = drugs[idx]

# Run the code "test_sclc.sh" to predict the Pred_SCLC
dir = "../data/ic50_data"
file = sprintf("%s/IC50_SCLC.txt", dir)
write.table(Pred_SCLC, file=file, sep="\t", row.names=F, col.names=T, quote=F)


pattern = "pred_sclc_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_SCLC = read_pred(dir, pattern=pattern)   # 6480 [648*10]
Pred_SCLC$Drug = Pred_SCLC$Drug %>% as.character

Pred_SCLC = Pred_SCLC %>% group_by(Cell, Drug) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 609

idx1 = match(Pred_SCLC$Cell, Anno_Cells$SANGER_MODEL_ID)
idx2 = match(Pred_SCLC$Drug, cids)

Pred_SCLC = Pred_SCLC %>% 
  mutate(Cell_Name=Anno_Cells$MODEL_NAME[idx1], Drug_Name=drugs[idx2]) %>% 
  relocate(Prediction, .after=everything()) %>% arrange(Drug)


# SIRT3, RHBDF1, GSDME
# TUBA1A, TUBA1B, TUBA1C, TUBB1, TUBB2A

dir = mkdir("Case Study [SCLC]")
drugs = c("Etoposide", "Cisplatin", "Lurbinectedin", 
          "Pyrrolobenzodiazepine", "Deruxtecan", "Calicheamicin")

targets = c("TOP2B", "AKT1", "POLR2A", "MCM2", "TOP1", "MCM2")
pathways = c("BIOCARTA_MET_PATHWAY", "BIOCARTA_AKT_PATHWAY", "BIOCARTA_G1_PATHWAY", 
             "BIOCARTA_MCM_PATHWAY", "BIOCARTA_AKT_PATHWAY", "BIOCARTA_G1_PATHWAY")

Pred_List_SCLC1 = list()
Pred_List_SCLC2 = list()

for (i in 1:length(drugs)) {
  main = sprintf("%s/%s [%s]", dir, drugs[i], targets[i])
  sprintf("\n# %s [%s]", drugs[i], targets[i]) %>% cat
  Pred_List_SCLC1[[drugs[i]]] = Pred_SCLC %>% 
    subset(Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_TPM, gene=targets[i], main=main, save=T)
  
  main = sprintf("%s/%s [%s]", dir, drugs[i], pathways[i])
  sprintf("# %s [%s]", drugs[i], pathways[i]) %>% cat
  Pred_List_SCLC2[[drugs[i]]] = Pred_SCLC %>% 
    subset(Drug_Name==drugs[i]) %>% 
    plot_with_tpm(SANGER_RNA_GSVA, gene=pathways[i], main=main, save=T)
}



### Grad-CAM [1st Layer]
seed = 2021:2030
dir = "../results/IC50_GDSC/Normal/RGCN"
file = sprintf("%s/pred_sclc_seed%s_gcam.csv", dir, seed)

Imp_SCLC = data.table()
for (i in 1:length(file)) {
  Importance = fread(file[i])
  Imp_SCLC = Imp_SCLC %>% rbind(Importance)
}

Imp_SCLC = Imp_SCLC[, lapply(.SD, mean), by = .(Cell, Drug)]

col = c("Cell", "Drug")
Imp_SCLC$Drug = Imp_SCLC$Drug %>% as.character
Imp_SCLC = left_join(Imp_SCLC, Pred_SCLC, by=col)

Imp_SCLC = Imp_SCLC %>% 
  relocate(Cell_Name, .after=Cell) %>% 
  relocate(Drug_Name, Prediction, .after=Drug)

get_imp_by_drug = function(Importance, dname=NULL, n_cell=30, n_path=5, resistant=F) {
  Importance = Importance %>% 
    subset(Drug_Name==dname) %>% 
    mutate(Z_Pred=as.numeric(scale(Prediction))) %>% 
    relocate(Z_Pred, .after=Prediction) %>% as.data.frame
    
  if (!resistant) {
    Importance = Importance %>% 
      slice_min(Prediction, n=n_cell) %>% 
      subset(Z_Pred<0) %>% arrange(Prediction)
  } else {
    Importance = Importance %>% 
      slice_max(Prediction, n=n_cell) %>% 
      subset(Z_Pred>0) %>% arrange(desc(Prediction))
  }
  
  pathway_list = c()
  for (i in 1:n_cell) {
    pathways = Importance[i, 7:ncol(Importance)] %>% 
      unlist %>% sort(decreasing=T) %>% head(n_path)
    pathway_list = pathway_list %>% c(pathways)
  }
  
  col = c("Cell_Name", "Drug_Name", "Prediction", "Z_Pred")
  Info = Importance[rep(1:n_cell, each=n_path), col]
  Info$Pathway = names(pathway_list)
  rownames(Info) = NULL
  return(Info)
}

GCAM_SCLC_Eto = Imp_SCLC %>% get_imp_by_drug(dname="Etoposide")
# BIOCARTA_AGPCR_PATHWAY
# BIOCARTA_GH_PATHWAY
# > GHR Knockdown - MEK/ERK pathway down - poliferation down, apoptosis - sensitivity up (GBM)
# > Zhang, Hongmei, et al. "The inhibition of GHR enhanced cytotoxic effects of etoposide on neuroblastoma." Cellular Signalling 86 (2021): 110081.
# BIOCARTA_SLRP_PATHWAY
# BIOCARTA_ERBB3_PATHWAY
# BIOCARTA_BARRESTIN_PATHWAY 

GCAM_SCLC_Cis = Imp_SCLC %>% get_imp_by_drug(dname="Cisplatin")
# BIOCARTA_CIRCADIAN_PATHWAY
# BIOCARTA_TNFR1_PATHWAY
# > Release TNF-α - TNF-α & TNFR1/2 interaction - cell death, inflammation
# > https://synapse.koreamed.org/articles/1050552
# > TNFR1 Overexpression - Resistant [NSCLC]
# > Zhen, Jie, et al. "EDN1 facilitates cisplatin resistance of non-small cell lung cancer cells by regulating the TNF signaling pathway." World Journal of Surgical Oncology 23.1 (2025): 71.
# BIOCARTA_GATA3_PATHWAY
# > GATA3 - Hippo down - Resistant [NB]
# > Wang, Jing, Wang Dai, and Ming Zhang. "GATA3 positively regulates PAR1 to facilitate in vitro disease progression and decrease cisplatin sensitivity in neuroblastoma via inhibiting the hippo pathway." Anti-Cancer Drugs 34.1 (2023): 57-72.
# BIOCARTA_ATM_PATHWAY
# > ATM inhibition - JAK/STAT3, PD-L1 down - EMT reversal, metastasis in resistant cells [SCLC]
# Shen, Mingjing, et al. "Inhibition of ATM reverses EMT and decreases metastatic potential of cisplatin-resistant lung cancer cells through JAK/STAT3/PD-L1 pathway." Journal of Experimental & Clinical Cancer Research 38 (2019): 1-14.
# BIOCARTA_BCELLSURVIVAL_PATHWAY 

GCAM_SCLC_Lur = Imp_SCLC %>% get_imp_by_drug(dname="Lurbinectedin")
# BIOCARTA_CCR5_PATHWAY [drug-gene X, sclc-gene M]
# > CCL5/CCR5 axis in SCLC and various cancer - progression, including angiogenesis, cell migration, and metastasis
# > Zeng, Zhen, et al. "CCL5/CCR5 axis in human diseases and related treatments." Genes & diseases 9.1 (2022): 12-27.
# BIOCARTA_IL12_PATHWAY
# BIOCARTA_CARM_ER_PATHWAY [drug-gene X, sclc-gene O]
# > ER - modulate ECM remodeling, cell adhesion - tumor progression, metastasis in SCLC
# > Wang, Hong, et al. "Therapeutic targeting ERRγ suppresses metastasis via extracellular matrix remodeling in small cell lung cancer." Embo Molecular Medicine 16.9 (2024): 2043-2059.
# BIOCARTA_SARS_PATHWAY
# BIOCARTA_RHODOPSIN_PATHWAY

GCAM_SCLC_PBD = Imp_SCLC %>% get_imp_by_drug(dname="Pyrrolobenzodiazepine")
# BIOCARTA_SLRP_PATHWAY [drug-gene X, sclc-gene M]
# > 
# > Ao, Zhi, et al. "Tumor angiogenesis of SCLC inhibited by decreased expression of FMOD via downregulating angiogenic factors of endothelial cells." Biomedicine & Pharmacotherapy 87 (2017): 539-547.
# BIOCARTA_HER2_PATHWAY
# > 
# BIOCARTA_REELIN_PATHWAY
# BIOCARTA_IL7_PATHWAY
# BIOCARTA_CCR5_PATHWAY [drug-gene M, sclc-gene M]
# > CCL5 - 

GCAM_SCLC_DXd = Imp_SCLC %>% get_imp_by_drug(dname="Deruxtecan")
# BIOCARTA_CIRCADIAN_PATHWAY
# BIOCARTA_TNFR1_PATHWAY
# > "blocking TNFα may enhance the sensitivity of HER2-positive breast cancer to trastuzumab deruxtecan"
# > https://pmc.ncbi.nlm.nih.gov/articles/PMC10016294
# BIOCARTA_GABA_PATHWAY
# BIOCARTA_PDZS_PATHWAY
# BIOCARTA_RAC1_PATHWAY [drug-gene X, sclc-gene O]
# > "inhibiting RAC1 can decrease SCLC cell viability and tumorigenicity, potentially enhancing the effectiveness of chemotherapy."
# > RAC1 inhibition - Nur77 from nucleus to cytoplasm - Nur77 bind to BCL2 - apoptosis
# > https://www.sciencedirect.com/science/article/pii/S2211124721014583

GCAM_SCLC_Cal = Imp_SCLC %>% get_imp_by_drug(dname="Calicheamicin")
# BIOCARTA_CCR5_PATHWAY
# BIOCARTA_RHODOPSIN_PATHWAY
# BIOCARTA_WNT_LRP6_PATHWAY [drug-gene X, sclc-gene O]
# > Wnt is activate in chemoresistant SCLC...
# > https://www.nature.com/articles/s41467-018-06162-9
# BIOCARTA_CARM_ER_PATHWAY
# BIOCARTA_FEEDER_PATHWAY

GCAM_SCLC_Eto = GCAM_SCLC_Eto %>% subset(Cell_Name=="LB647-SCLC")   # BIOCARTA_P38MAPK_PATHWAY
GCAM_SCLC_Cis = GCAM_SCLC_Cis %>% subset(Cell_Name=="NCI-H748")     # BIOCARTA_MET_PATHWAY
GCAM_SCLC_Lur = GCAM_SCLC_Lur %>% subset(Cell_Name=="NCI-H847")     # X
GCAM_SCLC_PBD = GCAM_SCLC_PBD %>% subset(Cell_Name=="COR-L279")     # BIOCARTA_EPONFKB_PATHWAY
GCAM_SCLC_DXd = GCAM_SCLC_DXd %>% subset(Cell_Name=="NCI-H209")     # BIOCARTA_TNFR1_PATHWAY
GCAM_SCLC_Cal = GCAM_SCLC_Cal %>% subset(Cell_Name=="IST-SL1")      # BIOCARTA_EPONFKB_PATHWAY

GCAM_SCLC = Reduce(rbind, list(GCAM_SCLC_Eto, GCAM_SCLC_Cis, GCAM_SCLC_Lur, 
                               GCAM_SCLC_PBD, GCAM_SCLC_DXd, GCAM_SCLC_Cal))

file = "Grad-CAM [SCLC, Well-known & Promising].csv"
write.csv(GCAM_SCLC, file=file, row.names=F)

# drugs = Imp_SCLC$Drug %>% unique
# col = c("Cell", "Prediction")
# get_head_tail = function(df) {
#   rbind(head(df), tail(df))
# }
# 

for (drug in drugs) {
  Temp = Imp_SCLC %>% subset(Drug_Name==drug) %>% arrange(Prediction) %>% as.data.frame
  Anno_Col = data.frame(Prediction=scale(Temp$Prediction))
  rownames(Temp) = Temp$Cell
  rownames(Anno_Col) = Temp$Cell
  
  Temp[, 6:ncol(Temp)] %>% t %>%
    heatmap_def(show_row=F, show_col=F, Anno_Col=Anno_Col, clust_col=T, scale_col=F)
  # color_bottom="white", color_center="palegoldenrod", color_top="firebrick3"
}



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
  
  # From 139
  Pred_COREAD1[[1]]$Rest %>% sum      # 91
  Pred_COREAD1[[2]]$Rest %>% sum      # 91
  Pred_COREAD1[[3]]$Rest %>% sum      # 91
  Pred_COREAD1[[4]]$Rest %>% sum      # 91
  
  Pred_COREAD1[[1]]$Missing %>% sum   # 0
  Pred_COREAD1[[2]]$Missing %>% sum   # 0
  Pred_COREAD1[[3]]$Missing %>% sum   # 1
  Pred_COREAD1[[4]]$Missing %>% sum   # 0
  
  # From 61
  Pred_Breast1[[1]]$Rest %>% sum      # 10
  Pred_Breast1[[2]]$Rest %>% sum      # 10
  Pred_Breast1[[3]]$Rest %>% sum      # 10
  Pred_Breast1[[4]]$Rest %>% sum      # 10
  
  Pred_Breast1[[1]]$Missing %>% sum   # 1
  Pred_Breast1[[2]]$Missing %>% sum   # 0
  Pred_Breast1[[3]]$Missing %>% sum   # 0
  Pred_Breast1[[4]]$Missing %>% sum   # 0
  
  
  ### [Source Data] Fig. 7
  Pred_COREAD1_ = Reduce(rbind, Pred_COREAD1)
  Pred_COREAD2_ = Reduce(rbind, Pred_COREAD2)
  Pred_Breast1_ = Reduce(rbind, Pred_Breast1)
  Pred_Breast2_ = Reduce(rbind, Pred_Breast2)
  Pred_SCLC1_ = Reduce(rbind, Pred_List_SCLC1)
  Pred_SCLC2_ = Reduce(rbind, Pred_List_SCLC2)
  
  col = c("Cell", "Cell_Name", "Drug", "Drug_Name", 
          "LN_IC50", "Prediction", "Gene_Name", "Gene")
  
  Pred_COREAD1_ = Pred_COREAD1_[, col] %>% rename(Gene_Log2_TPM=Gene)
  Pred_Breast1_ = Pred_Breast1_[, col] %>% rename(Gene_Log2_TPM=Gene)
  Pred_COREAD2_ = Pred_COREAD2_[, col] %>% rename(Pathway_Name=Gene_Name, Pathway_GSVA_Score=Gene)
  Pred_Breast2_ = Pred_Breast2_[, col] %>% rename(Pathway_Name=Gene_Name, Pathway_GSVA_Score=Gene)
  
  col_ = c("Pathway_Name", "Pathway_GSVA_Score")
  identical(Pred_COREAD1_[, col[1:6]], Pred_COREAD2_[, col[1:6]])   # T
  identical(Pred_Breast1_[, col[1:6]], Pred_Breast2_[, col[1:6]])   # T
  Pred_COREAD_ = Pred_COREAD1_ %>% cbind(Pred_COREAD2_[, col_])
  Pred_Breast_ = Pred_Breast1_ %>% cbind(Pred_Breast2_[, col_])
  
  Temp = list(Pred_COREAD_, Pred_Breast_)
  Temp %>% save_for_nc(num=7, suppl=F)
  
  
  ### [Source Data] (Fig. 8)
  col_ = c("Cell", "Cell_Name", "Drug", "Drug_Name", 
           "Prediction", "Gene_Name", "Gene")
  
  by = c("Cell", "Drug")
  col2 = c("Cell", "Drug", "LN_IC50")
  
  IC50_GDSC$Drug = IC50_GDSC$Drug %>% as.character
  Pred_SCLC1_$Drug = Pred_SCLC1_$Drug %>% as.character
  Pred_SCLC2_$Drug = Pred_SCLC2_$Drug %>% as.character
  Pred_SCLC1_ = Pred_SCLC1_[, col_] %>% left_join(IC50_GDSC[, col2], by=by) %>% relocate(LN_IC50, .before=Prediction)
  Pred_SCLC2_ = Pred_SCLC2_[, col_] %>% left_join(IC50_GDSC[, col2], by=by) %>% relocate(LN_IC50, .before=Prediction)
  
  Pred_SCLC1_ = Pred_SCLC1_[, col] %>% rename(Gene_Log2_TPM=Gene)
  Pred_SCLC2_ = Pred_SCLC2_[, col] %>% rename(Pathway_Name=Gene_Name, Pathway_GSVA_Score=Gene)
  
  col_ = c("Pathway_Name", "Pathway_GSVA_Score")
  identical(Pred_SCLC1_[, col[1:6]], Pred_SCLC2_[, col[1:6]])   # T
  Pred_SCLC_ = Pred_SCLC1_ %>% cbind(Pred_SCLC2_[, col_])
  
  file1 = "SCLC_in_vitro.csv"
  write.csv(Pred_SCLC_, file=file1, row.names=F)
}

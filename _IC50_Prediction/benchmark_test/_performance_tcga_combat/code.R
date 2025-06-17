#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

source("../functions.R")
source("functions.R")
loadings()

# processed_data/cell_data/TCGA/TCGA_Drug_Info.csv
file = "TCGA_Drug_Info.csv"
TCGA_Drug_Info = read.csv(file)   # 16 x 8




##### 2-1. Open Prediction Files

pattern = "pred_tcga_seed([0-9]+)_combat.csv"

dir_list = "../RF/Results/IC50_GDSC/Normal"
Pred_RF = read_pred_all(dir_list, "RF", pattern)         # 4140 x 12

dir_list = "../HiDRA/Results/IC50_GDSC/Normal"
Pred_HiDRA = read_pred_all(dir_list, "HiDRA", pattern)         # 4140 x 12

dir_list = "../PaccMann_SANGER/Results/IC50_GDSC/Normal"
Pred_PaccMann_SG = read_pred_all(dir_list, "PaccMann_SG", pattern)   # 4140 x 12

dir_list = "../TGSA/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP = read_pred_all(dir_list, "TGDRP", pattern)   # 3630 x 12

dir_list = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP_SG = read_pred_all(dir_list, "TGDRP_SG", pattern)   # 3630 x 12

dir_list = "../DRPreter/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter = read_pred_all(dir_list, "DRPreter", pattern)      # 4140 x 12

dir_list = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG = read_pred_all(dir_list, "DRPreter_SG", pattern)   # 4140 x 12

dir_list = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath = read_pred_all(dir_list, "GCNPath", pattern)       # 4140 x 12

model_list = c("RF", "HiDRA", "PaccMann_SG", "TGDRP", "TGDRP_SG",
               "DRPreter", "DRPreter_SG", "GCNPath")

Pred_List = Reduce(rbind, list(Pred_RF, Pred_HiDRA, Pred_PaccMann_SG, 
                               Pred_TGDRP, Pred_TGDRP_SG, 
                               Pred_DRPreter, Pred_DRPreter_SG, Pred_GCNPath))

idx = match(Pred_List$Drug_CID, TCGA_Drug_Info$Drug_CID)

Pred_List = Pred_List %>% 
  mutate(Drug_Name=TCGA_Drug_Info$Drug_Name[idx], 
         GDSC=TCGA_Drug_Info$GDSC[idx])




##### 2-2. Responder vs Non-Responder [Boxplot]

dir = mkdir("Responder and Non-Responder [Boxplot]")
main = sprintf("%s/Drug & Effect_Size [%s]", dir, model_list)

for (i in 1:length(model_list)) {
  Pred_List %>% subset(Model==model_list[i]) %>% 
    boxplot_tcga(main=main[i], resp_class1=T, save=T)
}




##### 2-3. Tumor vs Non-Tumor [Boxplot]

Pred_List_TN = data.frame()
Summary_List_TN = data.frame()

for (i in 1:length(model_list)) {
  Pred = Pred_List %>% subset(Model==model_list[i])
  Pred_TN = Pred %>% select_pat_tn
  
  if (nrow(Pred_TN)==0) {
    sprintf("Matched normal samples were not examined in %s", model_list[i]) %>% print
  } else {
    dir = mkdir("Tumor and Normal [Boxplot]")
    main = sprintf("%s/Tumor vs Normal [%s]", dir, model_list[i])
    TN_Info = Pred_TN %>% boxplot_tn(resp_class1=T, main=main, save=T, width=9, height=15)
    
    diff_resp = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Responder"]
    diff_non = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Non-Responder"]
    wilcox.test(diff_resp, diff_non, alternative="less")$p.value
    
    Pred_List_TN = Pred_List_TN %>% rbind(Pred_TN)
    Summary_List_TN = Summary_List_TN %>% rbind(TN_Info)
  }
}

# Matched normal samples were not examined in TGDRP, TGDRP_SG




##### 2-4. Responder vs Non-Responder [Significance Test, Scatter-Plot]

Sig_List = data.frame()
for (i in 1:length(model_list)) {
  Pred = Pred_List %>% subset(Model==model_list[i])
  Sig = sig_test_all(Pred, resp_class1=T, tumor_only=T)
  idx = match(Sig$Drug_Name, TCGA_Drug_Info$Drug_Name)
  
  Sig = Sig %>% 
    mutate(Model=model_list[i],
           Freq=TCGA_Drug_Info$Total[idx], 
           GDSC=TCGA_Drug_Info$GDSC[idx])
  Sig_List = Sig_List %>% rbind(Sig)
  
  dir = mkdir("Performance Significance [Scatter Plot]")
  main = sprintf("%s/Pvalue & Effect_Size [%s]", dir, model_list[i])
  Sig %>% plot_effect_pval(main=main, save=T)
  
  col = c("Drug_Name", "Pval", "Effect_Size")
  sprintf("# %s", model_list[i]) %>% print
  Sig[, col] %>% subset(Effect_Size<0 & Pval<0.05) %>% arrange(Pval) %>% print
  cat("\n")
}

Sig_List = Sig_List %>% relocate(Model, .before=everything())

# [1] "# RF"
# Drug_Name         Pval Effect_Size
# 1  Gemcitabine 4.322947e-08 -0.48457931
# 2 Fluorouracil 3.399461e-06 -0.09983353
# 3  Carboplatin 1.431476e-03 -0.15672820
# 4   Leuprolide 3.782186e-02 -0.06042762
# 
# [1] "# HiDRA"
# [1] Drug_Name   Pval        Effect_Size
# <0 rows> (or 0-length row.names)
# 
# [1] "# PaccMann_SG"
# Drug_Name         Pval Effect_Size
# 1 Fluorouracil 0.0004387959  -0.5290353
# 
# [1] "# TGDRP"
# Drug_Name        Pval Effect_Size
# 1 Fluorouracil 0.006456483 -0.13296822
# 2  Doxorubicin 0.013741986 -0.16360444
# 3    Cisplatin 0.041837071 -0.02026168
# 4   Paclitaxel 0.042715488 -0.18750958
# 
# [1] "# TGDRP_SG"
# Drug_Name         Pval Effect_Size
# 1  Doxorubicin 4.939796e-05 -0.38795738
# 2 Fluorouracil 3.692906e-02 -0.07220503
# 
# [1] "# DRPreter"
# Drug_Name         Pval Effect_Size
# 1 Fluorouracil 7.417822e-08  -0.4327411
# 2  Gemcitabine 4.572216e-03  -0.3196752
# 3  Carboplatin 2.155464e-02  -0.2846355
# 
# [1] "# DRPreter_SG"
# Drug_Name         Pval Effect_Size
# 1 Fluorouracil 1.202691e-08  -0.3642870
# 2  Gemcitabine 6.644734e-04  -0.2974963
# 3  Dacarbazine 1.871366e-02  -0.1132759
# 
# [1] "# GCNPath"
# Drug_Name         Pval Effect_Size
# 1 Fluorouracil 3.942216e-15 -0.52194073
# 2 Temozolomide 1.408679e-06 -0.06788287
# 3    Sorafenib 1.279975e-03 -0.19799570
# 4  Gemcitabine 1.016131e-02 -0.25074725
# 5    Cisplatin 1.105756e-02 -0.14392705

file = "Prediction [TCGA].csv"
write.csv(Pred_List, file=file, row.names=F)

file = "Tumor vs Normal [TCGA].csv"
write.csv(Pred_List_TN, file=file, row.names=F)

file = "Tumor vs Normal [TCGA, Summary].csv"
write.csv(Summary_List_TN, file=file, row.names=F)




supplementary = F
if (supplementary) {
  save_for_nc = function(df_list, dir=".", num=1, num_fig=NULL, rowNames=F, suppl=T, source=T) {
    
    suppressMessages(library(openxlsx))
    is_list = inherits(df_list, "list")
    if (is_list & is.null(num_fig)) num_fig = letters[1:length(df_list)]
    
    # Source Data
    if (!suppl) {
      sheets = sprintf("Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Fig. %s", num)
      file = sprintf("%s/SourceData_Fig%s.xlsx", dir, num)
    } else {
      sheets = sprintf("Supplementary Fig. %s%s", num, num_fig)
      if (!is_list) sheets = sprintf("Supplementary Fig. %s", num)
      file = sprintf("%s/SourceData_SupplementaryFig%s.xlsx", dir, num)
    }
    
    # Supplementary Data
    if (!source) {
      if (length(num_fig)!=1) {
        stop("Supplementary Information supports only single sheet...")
      }
      sheets = sprintf("Supplementary Data %s", num)
      file = sprintf("%s/Supplementary Data %s.xlsx", dir, num)
    }
    
    write.xlsx(df_list, file=file, sheetName=sheets, rowNames=rowNames)
  }
  
  
  ### [Source Data] Supplementary Fig. 28
  Pred_List_ = Pred_List %>% 
    subset(select=-c(Dataset, Test_Type)) %>% 
    rename(Train_Seed=Seed, GDSC_Drug=GDSC, Cancer_Type=TCGA_Code) %>% 
    relocate(Cancer_Type, .after=Sample_Type) %>% 
    relocate(GDSC_Drug, .after=Drug_CID) %>% 
    relocate(Train_Seed, .after=everything()) %>% as.data.frame
  
  idx = match(Sig_List$Drug_Name, TCGA_Drug_Info$Drug_Name)
  
  Sig_List_ = Sig_List %>% 
    subset(select=-c(Test_Type, Hypothesis)) %>% 
    rename(GDSC_Drug=GDSC, Num_Total=Freq, 
           Minus_Log10_Pval=MLog10_Pval) %>% 
    mutate(Drug_CID=TCGA_Drug_Info$Drug_CID[idx], 
           Num_Pos=TCGA_Drug_Info$Resp_C1[idx], 
           Num_Neg=TCGA_Drug_Info$Non_C1[idx]) %>% 
    relocate(Drug_CID, GDSC_Drug, .after=Drug_Name) %>% 
    relocate(Num_Total, Num_Pos, Num_Neg, .after=everything())
  
  Summary_List_TN_ = Summary_List_TN %>% 
    subset(select=-c(Dataset, Test_Type)) %>% 
    rename(Pred_Normal=Normal, Pred_Tumor=Tumor, GDSC_Drug=GDSC,
           Train_Seed=Seed, Response_Class=Resp_Class, Cancer_Type=TCGA_Code) %>% 
    relocate(Cancer_Type, .after=Patient) %>% 
    relocate(GDSC_Drug, .after=Drug_CID) %>% 
    relocate(Train_Seed, .after=everything()) %>% as.data.frame
  
  Pred_TCGA = list(Pred_List_, Sig_List_, Summary_List_TN_)
  Pred_TCGA %>% save_for_nc(num=28, suppl=T)
}

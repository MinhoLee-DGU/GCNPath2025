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

pattern = "pred_tcga_seed([0-9]+).csv"

dir_list = "../HiDRA/Results/IC50_GDSC/Normal"
Pred_HiDRA = read_pred_all(dir_list, "HiDRA", pattern)         # 4140 x 12

dir_list = "../PaccMann_SANGER/Results/IC50_GDSC/Normal"
Pred_PaccMann_SG = read_pred_all(dir_list, "PaccMann_SG", pattern)   # 4140 x 12

dir_list = "../DRPreter/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter = read_pred_all(dir_list, "DRPreter", pattern)      # 4140 x 12

dir_list = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG = read_pred_all(dir_list, "DRPreter_SG", pattern)   # 4140 x 12

dir_list = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath = read_pred_all(dir_list, "GCNPath", pattern)       # 4140 x 12

model_list = c("HiDRA", "PaccMann_SG", "DRPreter", "DRPreter_SG", "GCNPath")

Pred_List = Reduce(rbind, list(Pred_HiDRA, Pred_PaccMann_SG, Pred_DRPreter, Pred_DRPreter_SG, Pred_GCNPath))
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
  
  dir = mkdir("Tumor and Normal [Boxplot]")
  main = sprintf("%s/Tumor vs Normal [%s]", dir, model_list[i])
  TN_Info = Pred_TN %>% boxplot_tn(resp_class1=T, main=main, save=T)
  
  diff_resp = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Responder"]
  diff_non = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Non-Responder"]
  wilcox.test(diff_resp, diff_non, alternative="less")$p.value
  
  Pred_List_TN = Pred_List_TN %>% rbind(Pred_TN)
  Summary_List_TN = Summary_List_TN %>% rbind(TN_Info)
}




##### 2-4. Responder vs Non-Responder [Significance Test, Scatter-Plot]

Sig_List = data.frame()
for (i in 1:length(model_list)) {
  Pred = Pred_List %>% subset(Model==model_list[i])
  Sig = sig_test_all(Pred, resp_class1=T, tumor_only=T)
  idx = match(Sig$Drug_Name, TCGA_Drug_Info$Drug_Name)
  
  Sig = Sig %>% 
    mutate(Freq=TCGA_Drug_Info$Total[idx], 
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

# GCNPath
# Drug_Name    Pval         Effect_Size
# Fluorouracil 6.305972e-15 -0.57408782
# Temozolomide 1.327700e-07 -0.07834665
# Sorafenib    4.902823e-04 -0.21196010
# Cisplatin    1.173591e-02 -0.15045287
# Gemcitabine  1.506262e-02 -0.22940618

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
  
  
  ### [Source Data] Fig. 6
  Pred_GCNPath_ = Pred_GCNPath %>% 
    subset(select=-c(Model, Dataset, Test_Type)) %>% 
    rename(Train_Seed=Seed, GDSC_Drug=GDSC, Cancer_Type=TCGA_Code) %>% 
    relocate(Cancer_Type, .after=Sample_Type) %>% 
    relocate(GDSC_Drug, .after=Drug_CID) %>% 
    relocate(Train_Seed, .after=everything()) %>% as.data.frame
  
  idx = match(Sig_GCNPath$Drug_Name, TCGA_Drug_Info$Drug_Name)
  
  Sig_GCNPath_ = Sig_GCNPath %>% 
    subset(select=-c(Test_Type, Hypothesis)) %>% 
    rename(GDSC_Drug=GDSC, Num_Total=Freq, 
           Minus_Log10_Pval=MLog10_Pval) %>% 
    mutate(Drug_CID=TCGA_Drug_Info$Drug_CID[idx], 
           Num_Pos=TCGA_Drug_Info$Resp_C1[idx], 
           Num_Neg=TCGA_Drug_Info$Non_C1[idx]) %>% 
    relocate(Drug_CID, GDSC_Drug, .after=Drug_Name) %>% 
    relocate(Num_Total, Num_Pos, Num_Neg, .after=everything())
  
  TN_Info_ = TN_Info %>% 
    subset(select=-c(Model, Dataset, Test_Type)) %>% 
    rename(Pred_Normal=Normal, Pred_Tumor=Tumor, GDSC_Drug=GDSC,
           Train_Seed=Seed, Response_Class=Resp_Class, Cancer_Type=TCGA_Code) %>% 
    relocate(Cancer_Type, .after=Patient) %>% 
    relocate(GDSC_Drug, .after=Drug_CID) %>% 
    relocate(Train_Seed, .after=everything()) %>% as.data.frame
  
  Pred_TCGA = list(Pred_GCNPath_, Sig_GCNPath_, TN_Info_)
  Pred_TCGA %>% save_for_nc(num=6, suppl=F)
  
  
  ### [Source Data] Supplementary Fig. 27
  col = c("Database", "TCGA_Code", "PC1", "PC2")
  PCA_Raw_ = PCA_Raw[, col] %>% rename(Cancer_Type=TCGA_Code)
  PCA_GSVA_ = PCA_GSVA[, col] %>% rename(Cancer_Type=TCGA_Code)
  
  num_fig = c("a-c", "d-f")
  PCA_Plot_ = list(PCA_Raw_, PCA_GSVA_)
  PCA_Plot_ %>% save_for_nc(num=27, num_fig=num_fig, suppl=T)
  
  
  ### Prediction
  Pred_GCNPath_ = Pred_GCNPath %>% 
    subset(select=-c(Model, Dataset, Test_Type)) %>% 
    rename(Train_Seed=Seed, GDSC_Drug=GDSC)
  
  Pred_GCNPath_ = Pred_GCNPath_ %>% 
    relocate(Train_Seed, .before=Prediction) %>% 
    relocate(TCGA_Code, .after=Sample_Type) %>% 
    relocate(GDSC_Drug, .after=Drug_CID) %>% 
    arrange(Train_Seed, Drug_CID, Patient, Sample) %>% as.data.frame
  
  file = "Prediction [TCGA].csv"
  fwrite(Pred_GCNPath_, file=file, row.names=F)
}

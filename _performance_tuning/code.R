#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))

source("../utils/functions.R")
source("functions.R")
loadings()



### GCNPath Ablation Test

dir_res = "../results"
dir_fig = mkdir("Ablation Test")
dir_test = c("Normal", "Strict_Blind")
dir_test = sprintf("IC50_GDSC2/%s", dir_test)

### Cell Layer
clayer = c("RGCN", "RGCN_KNN3", "RGCN_KNN7", 
           "Pert_Seed2021", "RGCN_STRING_Physical", "RGCN_STRING_Rest", 
           "KEGG", "WikiPathways", "C4CM", "Linear_GAT")

re_label = list("RGCN"="BIOCARTA_KNN5 [DE]", 
                "RGCN_KNN3"="BIOCARTA_KNN3", 
                "RGCN_KNN7"="BIOCARTA_KNN7", 
                "Pert_Seed2021"="BIOCARTA_KNN5_Pert",
                "RGCN_STRING_Physical"="BIOCARTA_KNN5_STR9P",
                "RGCN_STRING_Rest"="BIOCARTA_KNN5_STR9R",
                "KEGG"="KEGG_KNN5", 
                "WikiPathways"="WikiPathways_KNN5", 
                "C4CM"="C4CM_KNN5", 
                "Linear_GAT"="BIOCARTA_Linear")

dir_test_c = c("Normal", "Cell_Blind", "Strict_Blind")
dir_test_c = sprintf("IC50_GDSC2/%s", dir_test_c)
dir_list = to_dir_list(dir_res, dir_test_c, clayer)
Pred_Cell = summary_pred(dir_list, lvl_param=clayer, re_label=re_label)


group = c(rep("GCN", length(clayer)-1), "FCN")
df_col = data.frame(Model=clayer, Group=group)
for (i in 1:length(re_label)) {
  df_col$Model[df_col$Model==names(re_label)[i]] = (re_label)[[i]]
}

lvl_group = c("GCN", "FCN")
idx = match(Pred_Cell$Perf_Test_Fold$Model, df_col$Model)
Pred_Cell$Perf_Test_Fold$Group = df_col$Group[idx] %>% factor(levels=lvl_group)

color = c("royalblue3", "#66b3ed")
names(color) = lvl_group

labels = lvl_group
names(labels) = lvl_group
add = list(scale_fill_manual(values=color, labels=labels), 
           labs(fill="Architecture\n[Cell Layer]"))

main = sprintf("%s/Ablation [Cell Layer]", dir_fig)
Pred_Cell$Perf_Test_Fold %>% boxplot_perf(main=main, width=45, height=15, save=T, 
                                          add=add, fill="Group", cell_blind=T, hide_legend=F)
Perf_Cell = perf_avg_sd(Pred_Cell$Perf_Test_Fold, space_cols=T)


### Drug Layer
dlayer = c("RGCN", "RGCN_MF128", "RGCN_MF256", "RGCN_MF512", "RGCN_MF1024", 
           "RGCN_MF256_Radius1", "RGCN_MF256_Radius3", 
           "RGCN_SMILESVec_Ch23Word", "RGCN_SMILESVec_PChWord")

re_label = list("RGCN"="GAT [DE]", 
                "RGCN_MF128"="MF128", "RGCN_MF256"="MF256", 
                "RGCN_MF512"="MF512", "RGCN_MF1024"="MF1024", 
                "RGCN_MF256_Radius1"="MF256_Radius1", 
                "RGCN_MF256_Radius3"="MF256_Radius3", 
                "RGCN_SMILESVec_Ch23Word"="SMILESVec_ChEMBL", 
                "RGCN_SMILESVec_PChWord"="SMILESVec_PubChem")

dir_test_d = c("Normal", "Drug_Blind", "Strict_Blind")
dir_test_d = sprintf("IC50_GDSC2/%s", dir_test_d)
dir_list = to_dir_list(dir_res, dir_test_d, dlayer)
Pred_Drug = summary_pred(dir_list, lvl_param=dlayer, re_label=re_label)


group = c("GCN", rep("FCN", length(dlayer)-1))
df_col = data.frame(Model=dlayer, Group=group)
for (i in 1:length(re_label)) {
  df_col$Model[df_col$Model==names(re_label)[i]] = (re_label)[[i]]
}

lvl_group = c("GCN", "FCN")
idx = match(Pred_Drug$Perf_Test_Fold$Model, df_col$Model)
Pred_Drug$Perf_Test_Fold$Group = df_col$Group[idx] %>% factor(levels=lvl_group)

color = c("royalblue3", "#66b3ed")
names(color) = lvl_group

labels = lvl_group
names(labels) = lvl_group
add = list(scale_fill_manual(values=color, labels=labels), 
           labs(fill="Architecture\n[Drug Layer]"))

main = sprintf("%s/Ablation [Drug Layer]", dir_fig)
Pred_Drug$Perf_Test_Fold %>% boxplot_perf(main=main, width=45, height=15, save=T, 
                                          add=add, fill="Group", drug_blind=T, hide_legend=F)
Perf_Drug = perf_avg_sd(Pred_Drug$Perf_Test_Fold, space_cols=T)




### Save all files...

suppressMessages(library(openxlsx))
file = sprintf("%s/Ablation Summary.xlsx", dir_fig)
write.xlsx(Perf_GNN_Lin, file=file, rowNames=T, colNames=T)


supplementary = T
if (supplementary) {
  process_perf = function(Pred, cell_blind=F, drug_blind=F) {
    Perf_Test = Pred$Perf_Test_Fold %>% subset(select=-c(Test)) %>% 
      rename(Parameter=Model, Train_Fold=Fold, Num_Test=N_Test)
    
    Perf_Test_G2N = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Normal")
    Perf_Test_G2S = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Strict_Blind")
    if (cell_blind) Perf_Test_G2C = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Cell_Blind")
    if (drug_blind) Perf_Test_G2D = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Drug_Blind")
    
    Perf_List = list(Perf_Test_G2N)
    if (cell_blind) Perf_List[[length(Perf_List)+1]] = Perf_Test_G2C
    if (drug_blind) Perf_List[[length(Perf_List)+1]] = Perf_Test_G2D
    Perf_List[[length(Perf_List)+1]] = Perf_Test_G2S
    Perf_List = Perf_List %>% lapply(function(df) df %>% subset(select=-c(Dataset, Test_Type)))
    
    params = unique(as.character(Perf_Test_G2N$Parameter))
    # sprintf("# Num [Test] : %s, %s", nrow(Perf_Test_G2N), nrow(Perf_Test_G2S)) %>% print
    sprintf("# Param : %s", paste0(params, collapse=", ")) %>% print
    return(Perf_List)
  }
  
  save_for_nc = function(df_list, dir=".", num=1, num_fig=NULL, rowNames=F, suppl=T) {
    suppressMessages(library(openxlsx))
    if (is.null(num_fig)) num_fig = letters[1:length(df_list)]
    
    if (!suppl) {
      sheets = sprintf("Fig. %s%s", num, num_fig)
      file = sprintf("%s/SourceData_Fig%s.xlsx", dir, num)
    } else {
      sheets = sprintf("Supplementary Fig. %s%s", num, num_fig)
      file = sprintf("%s/SourceData_SupplementaryFig%s.xlsx", dir, num)
    }
    
    write.xlsx(df_list, file=file, sheetName=sheets, rowNames=rowNames)
  }
  
  ### [Source Data] Fig. 2
  Perf_GNN_Lin_ = Pred_GNN_Lin %>% process_perf
  Perf_GNN_Lin_[[1]]$Group = NULL
  Perf_GNN_Lin_[[2]]$Group = NULL
  Perf_GNN_Lin_ %>% save_for_nc(num=2, suppl=F)
  
  Perf_Cell_ = Pred_Cell %>% process_perf(cell_blind=T)
  Perf_Drug_ = Pred_Drug %>% process_perf(drug_blind=T)
  Perf_List_ = c(Perf_Cell_, Perf_Drug_)
  Perf_List_ %>% save_for_nc(num=2, suppl=F)
}

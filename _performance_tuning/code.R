#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))

source("../utils/functions.R")
loadings()

read_pred = function(dir, pattern, sep=",") {
  
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_.*([0-9]+$)\\.csv", "Pred_CCLE_.*([0-9]+)\\.csv"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      try({
        Pred_TP = fread(df_name, header=T, sep=sep)
        df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
        nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
        model = strsplit(dir, "/")[[1]] %>% tail(1)
        
        Pred_TP$Fold = nth
        Pred_TP$Model = model
        Pred = Pred %>% rbind(Pred_TP)
      })
    }
    return(Pred)
  }
}

read_pred_all = function(dir_list, pattern, sep=",") {
  
  Pred = data.frame()
  for (dir in dir_list) {
    Pred_TP = read_pred(dir, pattern, sep=sep)
    if (!is.null(Pred_TP)) {
      dir_detail = strsplit(dir, "/") %>% unlist
      
      if ("IC50_GDSC" %in% dir_detail) {
        dataset = "GDSC"
      } else if ("IC50_GDSC1" %in% dir_detail) {
        dataset = "GDSC1"
      } else if ("IC50_GDSC2" %in% dir_detail) {
        dataset = "GDSC2"
      } else if ("IC50_CCLE" %in% dir_detail) {
        dataset = "CCLE"
      } else print("Error!!!")
      
      if ("Normal" %in% dir_detail) {
        test_type = "Normal"
      } else if ("Cell_Blind" %in% dir_detail) {
        test_type = "Cell_Blind"
      } else if ("Drug_Blind" %in% dir_detail) {
        test_type = "Drug_Blind"
      } else if ("Strict_Blind" %in% dir_detail) {
        test_type = "Strict_Blind"
      } else print("Error!!!")
      
      Pred_TP$Dataset = dataset
      Pred_TP$Test_Type = test_type
      Pred_TP = Pred_TP %>% relocate(Model, Dataset, Test_Type, Fold, .before=everything())
      Pred = Pred %>% rbind(Pred_TP)
    }
  }
  return(Pred)
}

calc_perf = function(Pred) {
  
  Perf = Pred %>% group_by(Model, Dataset, Test_Type) %>%
    summarize(N_Test = n(),
              RMSE = RMSE(LN_IC50, Prediction),
              R2 = R2(LN_IC50, Prediction),
              Corr = cor(LN_IC50, Prediction)) %>% as.data.frame
  
  lvl = c("GDSC2 x Normal", "GDSC2 x Strict_Blind", "GDSC x Strict_Blind")
  Perf = Perf %>% mutate(Test = sprintf("%s x %s", Dataset, Test_Type))
  Perf$Test = Perf$Test %>% factor(levels=lvl)
  return(Perf)
}

calc_perf_fold = function(Pred) {
  
  Perf = Pred %>% group_by(Model, Dataset, Test_Type, Fold) %>%
    summarize(N_Test = n(),
              RMSE = RMSE(LN_IC50, Prediction),
              R2 = R2(LN_IC50, Prediction),
              Corr = cor(LN_IC50, Prediction)) %>% as.data.frame
  
  lvl = c("GDSC2 x Normal", "GDSC2 x Strict_Blind", "GDSC x Strict_Blind")
  Perf = Perf %>% mutate(Test=sprintf("%s x %s", Dataset, Test_Type))
  Perf$Test = Perf$Test %>% factor(levels=lvl)
  return(Perf)
}


### GCNPath Ablation Test

to_dir_list = function(dir, dir1, dir2) {
  dir_list = expand.grid(dir1, dir2) %>% arrange(Var1)
  dir_list = sprintf("%s/%s/%s", dir, dir_list$Var1, dir_list$Var2)
  return(dir_list)
}

# boxplot_perf = function(Perf_Fold, main=NULL, width=27, height=18, test=T, save=T) {
#   subtitle = ifelse(test, "Valid", "Test")
#   main = sprintf("%s [%s]", main, subtitle)
#   Perf_Fold %>% boxplot_def(Model, RMSE, fill=Test, main=main, width=width, height=height, save=save)
# }

boxplot_perf = function(Perf_Fold, main=NULL, re_label=NULL, 
                        width=48, height=13.5, axis_tl=30, axis_tx=18, save=T) {
  
  # legend = "Parameters"
  # p1 = Perf_Fold %>% subset(Test=="GDSC2 x Normal") %>% 
  #   boxplot_def(Model, RMSE, fill=Model, main=main, legend=legend, width=width/5.2, height=height, save=F)
  # p2 = Perf_Fold %>% subset(Test=="GDSC2 x Strict_Blind") %>% 
  #   boxplot_def(Model, RMSE, fill=Model, main=main, legend=legend, width=width/5.2, height=height, save=F)
  # p3 = Perf_Fold %>% subset(Test=="GDSC x Strict_Blind") %>% 
  #   boxplot_def(Model, RMSE, fill=Model, main=main, legend=legend, width=width/5.2, height=height, save=F)
  
  add.params = list(alpha=0.5)
  margin = margin(8, 8, 8, 8, unit="pt")
  font1 = font("ylab", color="grey30", size=axis_tl, margin=margin)
  font2 = font("axis.text", color="grey30", size=axis_tx, margin=margin)
  
  font = font1 + font2
  margin_pl = margin(0.2, 0.4, 0.2, 0.4, "cm")
  margin_pl = theme(plot.margin=margin_pl)
  
  p1 = Perf_Fold %>% subset(Test=="GDSC2 x Normal") %>% 
    ggboxplot("Model", y="RMSE", fill="Model", add="point", 
              xlab=F, add.params=add.params) + font + margin_pl + rotate_x_text(45)
  p2 = Perf_Fold %>% subset(Test=="GDSC2 x Strict_Blind") %>% 
    ggboxplot("Model", y="RMSE", fill="Model", add="point", 
              xlab=F, add.params=add.params) + font + margin_pl + rotate_x_text(45)
  p3 = Perf_Fold %>% subset(Test=="GDSC x Strict_Blind") %>% 
    ggboxplot("Model", y="RMSE", fill="Model", add="point", 
              xlab=F, add.params=add.params) + font + margin_pl + rotate_x_text(45)
  
  no_lg = function(pl) pl %>% ggpar(legend="none")
  pl = ggarrange(no_lg(p1), no_lg(p2), no_lg(p3), ncol=3)
  # common.legend=T, legend="right"
  
  if (save) {
    save_fig(pl, main=main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
}

summary_pred = function(dir_list, lvl_param=NULL, re_label=NULL) {
  
  # relabel : Re-label the parameter names
  # Ex. (vector)
  # Ex. (list)
  
  pattern_val = "pred_valid_([0-9]+).csv"
  pattern_test = "pred_test_([0-9]+).csv"
  Pred_Val = read_pred_all(dir_list, pattern_val)
  Pred_Test = read_pred_all(dir_list, pattern_test)
  
  if (!is.null(re_label)) {
    for (i in 1:length(re_label)) {
      re_label_ = ifelse(is.list(re_label), re_label[[i]], re_label[i])
      lvl_param[lvl_param==names(re_label)[i]] = re_label_
      Pred_Val$Model[Pred_Val$Model==names(re_label)[i]] = re_label_
      Pred_Test$Model[Pred_Test$Model==names(re_label)[i]] = re_label_
    }
  }
  
  if (!is.null(lvl_param)) {
    lvl_param = lvl_param %>% strsplit("/") %>% sapply(function(x) x[length(x)])
    Pred_Val$Model = Pred_Val$Model %>% factor(levels=lvl_param)
    Pred_Test$Model = Pred_Test$Model %>% factor(levels=lvl_param)
  }
  
  Perf_Val = Pred_Val %>% calc_perf
  Perf_Test = Pred_Test %>% calc_perf
  Perf_Val_Fold = Pred_Val %>% calc_perf_fold
  Perf_Test_Fold = Pred_Test %>% calc_perf_fold
  
  Pred_List = list(Pred_Val=Pred_Val, Perf_Val=Perf_Val, Perf_Val_Fold=Perf_Val_Fold,
                   Pred_Test=Pred_Test, Perf_Test=Perf_Test, Perf_Test_Fold=Perf_Test_Fold)
  
  return(Pred_List)
}

  perf_avg_sd = function(Perf_Fold, space_cols=F) {
  
  avg_plus_sd = function(x_mean, x_sd) sprintf("%.3f±%.3f", x_mean, x_sd)
  Perf_Fold_Avg = Perf_Fold %>% acast(Model~Test, value.var="RMSE", fun.aggregate=mean)
  Perf_Fold_SD = Perf_Fold %>% acast(Model~Test, value.var="RMSE", fun.aggregate=sd)
  
  Perf_Fold_Avg = Perf_Fold_Avg %>% as.data.frame
  Perf_Fold_SD = Perf_Fold_SD %>% as.data.frame
  Perf_Fold_Sum = mapply(avg_plus_sd, Perf_Fold_Avg, Perf_Fold_SD)
  if (space_cols) colnames(Perf_Fold_Sum) = gsub(" x ", "\n", colnames(Perf_Fold_Sum))
  
  Perf_Fold_Sum = Perf_Fold_Sum %>% as.data.frame
  rownames(Perf_Fold_Sum) = rownames(Perf_Fold_Avg)
  return(Perf_Fold_Sum)
}

valid_test = c("Valid", "Test")
dir_test = c("IC50_GDSC2/Normal", "IC50_GDSC2/Strict_Blind", "IC50_GDSC/Strict_Blind")

dir_res = "../results"
dir_fig = mkdir("Ablation Test")

Abl_Test = list()
Abl_Valid = list()



# Ablation Test [Network Type]
network_no = "Linear_GAT"
network_uni = c("STR7", "STR9", "Reg", "Corr", "Random")
network_uni = sprintf("GAT_%s", network_uni)
network_multi = c("", "_Inv", "_KNN3", "_KNN7", "_UD")
network_multi = sprintf("RGCN%s", network_multi)
network = c(network_no, network_uni, network_multi)

re_label = list("RGCN"="RGCN [DE]", "Linear_GAT"="Linear")
dir_list = to_dir_list(dir_res, dir_test, network)
Pred_Net = summary_pred(dir_list, lvl_param=network, re_label=re_label)

main = sprintf("%s/Ablation [Network, %s]", dir_fig, valid_test)
Pred_Net$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
Pred_Net$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)

Abl_Valid$Network = perf_avg_sd(Pred_Net$Perf_Val_Fold, space_cols=T)
Abl_Test$Network = perf_avg_sd(Pred_Net$Perf_Test_Fold, space_cols=T)


# Ablation Test [Architecture]
attn = c("Attn_v1", "Attn_v2")
attn_param = c("", "_C16", "_AL3", "_A64", "_AH8")
architecture = expand.grid(attn, attn_param) %>% arrange(Var1)
architecture = sprintf("%s%s", architecture$Var1, architecture$Var2)
architecture = "RGCN" %>% c(architecture)

re_label = list("RGCN"="No Attn [DE]")
dir_list = to_dir_list(dir_res, dir_test, architecture)
Pred_Arch = summary_pred(dir_list, lvl_param=architecture, re_label=re_label)

main = sprintf("%s/Ablation [Architecture, %s]", dir_fig, valid_test)
Pred_Arch$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, axis_tx=16.5, save=T)
Pred_Arch$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, axis_tx=16.5, save=T)

Abl_Valid$Architect = perf_avg_sd(Pred_Arch$Perf_Val_Fold, space_cols=T)
Abl_Test$Architect = perf_avg_sd(Pred_Arch$Perf_Test_Fold, space_cols=T)


# # Ablation Test [Architecture, SSI-DDI]
# attn = c("SSI_Attn")
# attn_param = c("", "_AL1", "_A32", "_AL1_A32", "_DL5", "_DL5_AL1_A32")
# architecture = expand.grid(attn, attn_param) %>% arrange(Var1)
# architecture = sprintf("%s%s", architecture$Var1, architecture$Var2)
# architecture = "RGCN" %>% c(architecture)
# 
# re_label = list("RGCN"="No Attn [DE]")
# dir_list = to_dir_list(dir_res, dir_test, architecture)
# Pred_SSI = summary_pred(dir_list, lvl_param=architecture, re_label=re_label)
# 
# main = sprintf("%s/Ablation [Architecture (SSI-DDI), %s]", dir_fig, valid_test)
# Pred_SSI$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
# Pred_SSI$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)
# 
# Abl_Valid$Architect_SSI = perf_avg_sd(Pred_SSI$Perf_Val_Fold, space_cols=T)
# Abl_Test$Architect_SSI = perf_avg_sd(Pred_SSI$Perf_Test_Fold, space_cols=T)


# # Ablation Test [Architecture, SA-DDI]
# attn = c("SA_Attn")
# attn_param = c("", "_A32", "_DL5", "_DL5_A32")
# architecture = expand.grid(attn, attn_param) %>% arrange(Var1)
# architecture = sprintf("%s%s", architecture$Var1, architecture$Var2)
# architecture = "RGCN" %>% c(architecture)
# 
# re_label = list("RGCN"="No Attn [DE]")
# dir_list = to_dir_list(dir_res, dir_test, architecture)
# Pred_SA = summary_pred(dir_list, lvl_param=architecture, re_label=re_label)
# 
# main = sprintf("%s/Ablation [Architecture (SA-DDI), %s]", dir_fig, valid_test)
# Pred_SA$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
# Pred_SA$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)
# 
# Abl_Valid$Architect_SA = perf_avg_sd(Pred_SA$Perf_Val_Fold, space_cols=T)
# Abl_Test$Architect_SA = perf_avg_sd(Pred_SA$Perf_Test_Fold, space_cols=T)


# Ablation Test [Pathway]
pathway = c("RGCN", "KEGG", "PID", "WikiPathways", "C4CM")

re_label = list("RGCN"="BIOCARTA [DE]")
dir_list = to_dir_list(dir_res, dir_test, pathway)
Pred_Pathway = summary_pred(dir_list, lvl_param=pathway, re_label=re_label)

main = sprintf("%s/Ablation [Pathway, %s]", dir_fig, valid_test)
Pred_Pathway$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=45, save=T)
Pred_Pathway$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=45, save=T)

Abl_Valid$Pathway = perf_avg_sd(Pred_Pathway$Perf_Val_Fold, space_cols=T)
Abl_Test$Pathway = perf_avg_sd(Pred_Pathway$Perf_Test_Fold, space_cols=T)



# Ablation Test [GNN vs Linear]
gnn_lin = c("RGCN", "RGAT", "RGCN_GIN", "RGCN_GINE", 
            "RGCN_MF256", "RGCN_MF512", "RGCN_SMILESVec", "Linear_GAT", "Linear_MF256")

re_label = list("RGCN"="RGCN_GAT [DE]", "RGAT"="RGAT_GAT")
dir_list = to_dir_list(dir_res, dir_test, gnn_lin)
Pred_GNN_Lin = summary_pred(dir_list, lvl_param=gnn_lin, re_label=re_label)

main = sprintf("%s/Ablation [GNN vs Linear, %s]", dir_fig, valid_test)
Pred_GNN_Lin$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
Pred_GNN_Lin$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)

Abl_Valid$GNN_Lin = perf_avg_sd(Pred_GNN_Lin$Perf_Val_Fold, space_cols=T)
Abl_Test$GNN_Lin = perf_avg_sd(Pred_GNN_Lin$Perf_Test_Fold, space_cols=T)



# # Ablation Test [GNN Stack]
# gnn_stack = c("RGCN", "Plain", "ResNet", "ResNet+", "JK_Concat")
# dim_drug = c("", "_D256")
# gnn_stack = expand.grid(gnn_stack, dim_drug) %>% arrange(Var1)
# gnn_stack = sprintf("%s%s", gnn_stack$Var1, gnn_stack$Var2)
# 
# re_label = list("RGCN"="DenseNet [DE]", "RGCN_D256"="DenseNet_D256")
# dir_list = to_dir_list(dir_res, dir_test, gnn_stack)
# Pred_GNN_Stack = summary_pred(dir_list, lvl_param=gnn_stack, re_label=re_label)
# 
# main = sprintf("%s/Ablation [GNN Stack, %s]", dir_fig, valid_test)
# Pred_GNN_Stack$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
# Pred_GNN_Stack$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)
# 
# Abl_Valid$GNN_Stack = perf_avg_sd(Pred_GNN_Stack$Perf_Val_Fold, space_cols=T)
# Abl_Test$GNN_Stack = perf_avg_sd(Pred_GNN_Stack$Perf_Test_Fold, space_cols=T)


# Ablation Test [Cell Dimension & Layer]
dim = c("", "_C4", "_C16", "_D64", "_D256", "_P256")
layer = c("_CL5", "_DL5", "_PL4")
dim_layer = sprintf("RGCN%s", c(dim, layer))

re_label = list("RGCN"="RGCN [DE]")
dir_list = to_dir_list(dir_res, dir_test, dim_layer)
Pred_DimLayer = summary_pred(dir_list, lvl_param=dim_layer, re_label=re_label)

main = sprintf("%s/Ablation [Dimension & Layer, %s]", dir_fig, valid_test)
Pred_DimLayer$Perf_Val_Fold %>% boxplot_perf(main=main[1], width=48, save=T)
Pred_DimLayer$Perf_Test_Fold %>% boxplot_perf(main=main[2], width=48, save=T)

Abl_Valid$Dim_Layer = perf_avg_sd(Pred_DimLayer$Perf_Val_Fold, space_cols=T)
Abl_Test$Dim_Layer = perf_avg_sd(Pred_DimLayer$Perf_Test_Fold, space_cols=T)


### Save all files...

suppressMessages(library(openxlsx))
file = sprintf("%s/Ablation Summary [%s].xlsx", dir_fig, valid_test)

write.xlsx(Abl_Valid, file=file[1], rowNames=T, colNames=T, sheetName=names(Abl_Valid))
write.xlsx(Abl_Test, file=file[2], rowNames=T, colNames=T, sheetName=names(Abl_Test))

file = "Ablation.RData"
save(Pred_Net, Pred_Arch, Pred_DimLayer, 
     Pred_GNN_Lin, Pred_Pathway, Abl_Valid, Abl_Test, file=file)

# rm(list=grep("Pred_", ls(), value=T))


supplementary = T
if (supplementary) {
  
  process_perf = function(Perf) {
    
    Perf_Val = Perf$Perf_Val_Fold %>% 
      subset(select=-c(Test, R2, Corr)) %>% 
      rename(Parameter=Model, Train_Fold=Fold, Num_Test=N_Test)
    
    Perf_Test = Perf$Perf_Test_Fold %>% 
      subset(select=-c(Test, R2, Corr)) %>% 
      rename(Parameter=Model, Train_Fold=Fold, Num_Test=N_Test)
    
    Perf_Val_G2N = Perf_Val %>% subset(Dataset=="GDSC2" & Test_Type=="Normal")
    Perf_Val_G2S = Perf_Val %>% subset(Dataset=="GDSC2" & Test_Type=="Strict_Blind") 
    Perf_Val_G0S = Perf_Val %>% subset(Dataset=="GDSC" & Test_Type=="Strict_Blind") 
    
    Perf_Test_G2N = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Normal")
    Perf_Test_G2S = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Strict_Blind") 
    Perf_Test_G0S = Perf_Test %>% subset(Dataset=="GDSC" & Test_Type=="Strict_Blind") 
    
    Perf_List = list(Perf_Val_G2N, Perf_Val_G2S, Perf_Val_G0S, 
                     Perf_Test_G2N, Perf_Test_G2S, Perf_Test_G0S) %>% 
      lapply(function(df) df %>% subset(select=-c(Dataset, Test_Type)))
    
    params = unique(as.character(Perf_Val_G2N$Parameter))
    sprintf("# Num [Val] : %s, %s, %s", nrow(Perf_Val_G2N), nrow(Perf_Val_G2S), nrow(Perf_Val_G0S)) %>% print
    sprintf("# Num [Test] : %s, %s, %s", nrow(Perf_Test_G2N), nrow(Perf_Test_G2S), nrow(Perf_Test_G0S)) %>% print
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
  
  
  ### [Source Data] Supplementary Fig. 29
  Perf_Net_ = Pred_Net %>% process_perf
  Perf_Net_ %>% save_for_nc(num=29, suppl=T)
  
  ### [Source Data] Supplementary Fig. 30
  Perf_GNN_Lin_ = Pred_GNN_Lin %>% process_perf
  Perf_GNN_Lin_ %>% save_for_nc(num=30, suppl=T)
  
  ### [Source Data] Supplementary Fig. 31
  Perf_DimLayer_ = Pred_DimLayer %>% process_perf
  Perf_DimLayer_ %>% save_for_nc(num=31, suppl=T)
  
  ### [Source Data] Supplementary Fig. 32
  Perf_Pathway_ = Pred_Pathway %>% process_perf
  Perf_Pathway_ %>% save_for_nc(num=32, suppl=T)
  
  ### [Source Data] Supplementary Fig. 33
  Perf_Attn_ = Pred_Arch %>% process_perf
  Perf_Attn_ %>% save_for_nc(num=33, suppl=T)
}

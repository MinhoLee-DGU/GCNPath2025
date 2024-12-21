#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))

source("../functions.R")
loadings()

options(dplyr.summarise.inform=F)

read_pred = function(dir, pattern, sep=",", seed=F) {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_.*([0-9]+$)\\.csv", "Pred_CCLE_.*([0-9]+)\\.csv"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      Pred_TP = fread(df_name, header=T, sep=sep)
      df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
      fold_nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
      Pred_TP$Fold = fold_nth
      
      if (seed) {
        seed_nth = strsplit(df_name, "/")[[1]] %>% tail(2) %>% head(1)
        seed_nth = gsub("Seed", "", seed_nth) %>% as.numeric
        Pred_TP$Seed = seed_nth
      }
      Pred = Pred %>% rbind(Pred_TP)
    }
    return(Pred)
  }
}

calc_perf = function(Pred, option=0, dataset=NULL, test_type=NULL, seed=F) {
  
  if (option==0) {
    col_by = "Fold"
  } else if (option==1) {
    col_by = c("Model", "Cell", "Dataset", "Test_Type")
  } else if (option==2) {
    col_by = c("Model", "Drug", "Dataset", "Test_Type")
  }
  
  if (seed) col_by = col_by %>% c("Seed")
  pred_inf = Pred$Prediction %>% is.infinite
  if (any(pred_inf)) {
    fold_inf = Pred[pred_inf, "Fold"] %>% 
      unique %>% unlist %>% sort %>% paste(collapse=" & ")
    
    Pred_Inf = Pred[pred_inf, ]
    Pred = Pred[!pred_inf, ]
    
    if (option==0) {
      text = sprintf("%s x %s", dataset, test_type)
      sprintf("Prediction Inf found in %s... [n=%s in Fold %s]", 
              text, sum(pred_inf), fold_inf) %>% print
    } else {
      sprintf("Prediction Inf found... [n=%s in Fold %s]", 
              sum(pred_inf), fold_inf) %>% print
    }
  }
  
  Perf = Pred %>% 
    group_by(across(all_of(col_by))) %>% 
    summarize(N_Test = n(),
              RMSE = RMSE(LN_IC50, Prediction),
              PCC = cor(LN_IC50, Prediction), 
              SCC = cor(LN_IC50, Prediction, method="spearman"))
  
  if (any(pred_inf)) {
    Perf_Inf = Pred_Inf %>% 
      group_by(across(all_of(col_by))) %>% 
      summarize(N_Test = n(), RMSE = NA, PCC = NA, SCC = NA)
    
    # Function distinct drop duplicatee rows except first one given a certain columns
    Perf = rbind(Perf, Perf_Inf) %>% as.data.frame %>% 
      distinct(across(all_of(col_by)), .keep_all=T) %>% as.data.frame
    if (option==0) Perf = Perf %>% arrange(Fold) %>% as.data.frame
  } else {
    Perf = Perf %>% as.data.frame
  }
  return(Perf)
}

read_pred_all = function(dir_list, model_name, pattern, sep=",", seed=F) {
  
  Pred = data.frame()
  Perf = data.frame()
  
  for (dir in dir_list) {
    Pred_TP = read_pred(dir, pattern, sep=sep, seed=seed)
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
      
      Perf_TP = Pred_TP %>% calc_perf(option=0, dataset, test_type, seed=seed)
      Perf_TP$Model = model_name
      Perf_TP$Dataset = dataset
      Perf_TP$Test_Type = test_type
      Perf_TP = Perf_TP %>% relocate(Model, Dataset, Test_Type, Fold, N_Test, .before=everything())
      Perf = Perf %>% rbind(Perf_TP)
      
      Pred_TP$Model = model_name
      Pred_TP$Dataset = dataset
      Pred_TP$Test_Type = test_type
      Pred_TP = Pred_TP %>% relocate(Model, Dataset, Test_Type, Fold, .before=everything())
      Pred = Pred %>% rbind(Pred_TP)
    }
  }
  
  Perf_Cell = Pred %>% calc_perf(option=1, seed=seed)
  Perf_Drug = Pred %>% calc_perf(option=2, seed=seed)
  
  if (seed) {
    Pred$Model = Pred %>% with(sprintf("%s_Seed%s", Model, Seed))
    Perf$Model = Perf %>% with(sprintf("%s_Seed%s", Model, Seed))
    Perf_Cell$Model = Perf_Cell %>% with(sprintf("%s_Seed%s", Model, Seed))
    Perf_Drug$Model = Perf_Drug %>% with(sprintf("%s_Seed%s", Model, Seed))
  }
  
  Pred = list(Pred=Pred, Perf=Perf, Perf_Cell=Perf_Cell, Perf_Drug=Perf_Drug)
  sprintf("# The number of Test : %s", nrow(Perf)) %>% print
  return(Pred)
}



##### 2. Get Performance Files [Overall]

is_inf_def = function(x) {
  inf_p = sum(x>0 & is.infinite(x))
  inf_m = sum(x<0 & is.infinite(x))
  if (inf_p!=0) sprintf("# Infinite (+) : %s", inf_p) %>% print
  if (inf_m!=0) sprintf("# Infinite (-) : %s", inf_m) %>% print
}

dir1 = sprintf("IC50_GDSC%s", c("", 1, 2))
dir2 = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")

Dir_List = expand.grid(dir1, dir2)
dir_list = Dir_List %>% apply(1, function(x) paste0(x, collapse="/"))

pattern = "pred_test_([0-9]+).csv"

# tCNNS [O]
dir = "../tCNNS/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_tCNNS = read_pred_all(dir_list_, "tCNNS", pattern)
# "Prediction Inf found in GDSC x Normal... [n=37180 in Fold 4]"
# "Prediction Inf found in GDSC2 x Normal... [n=19828 in Fold 3]"
# "Prediction Inf found in GDSC x Drug_Blind... [n=36970 in Fold 6]"
# "Prediction Inf found in GDSC2 x Drug_Blind... [n=36899 in Fold 6 & 8]"
# "Prediction Inf found in GDSC x Strict_Blind... [n=15086 in Fold 17]"
# "Prediction Inf found... [n=145963 in Fold 3 & 4 & 6 & 8 & 17]"

# Pred_tCNNS$Perf$RMSE %>% is.na %>% sum   # 6
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Normal" & Fold==4) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Drug_Blind" & Fold==6) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC" & Test_Type=="Strict_Blind" & Fold==17) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Normal" & Fold==3) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Drug_Blind" & Fold==6) %>% pull(Prediction) %>% is_inf_def
# Pred_tCNNS$Pred %>% subset(Dataset=="GDSC2" & Test_Type=="Drug_Blind" & Fold==8) %>% pull(Prediction) %>% is_inf_def

# HiDRA [O]
dir = "../HiDRA/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_HiDRA = read_pred_all(dir_list_, "HiDRA", pattern)

# PaccMann [O]
dir = "../PaccMann/Results"

dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_PaccMann = read_pred_all(dir_list_, "PaccMann", pattern)

# PaccMann_SANGER [O]
dir = "../PaccMann_SANGER/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_PaccMann_SG = read_pred_all(dir_list_, "PaccMann_SG", pattern)

# GraphDRP [O]
dir = "../GraphDRP/Results"
dir_list_ = sprintf("%s/%s", dir, dir_list)
Pred_GraphDRP = read_pred_all(dir_list_, "GraphDRP", pattern)

# TGDRP & TGSA [O/O]
dir = "../TGSA/Results"
dir_list_tgdrp = sprintf("%s/%s/TGDRP", dir, dir_list)
dir_list_tgsa = sprintf("%s/%s/TGSA", dir, dir_list)
Pred_TGDRP = read_pred_all(dir_list_tgdrp, "TGDRP", pattern)
Pred_TGSA = read_pred_all(dir_list_tgsa, "TGSA", pattern)

# TGDRP_SANGER & TGSA_SANGER [O/O]
dir = "../TGSA_SANGER/Results"
dir_list_tgdrp = sprintf("%s/%s/TGDRP", dir, dir_list)
dir_list_tgsa = sprintf("%s/%s/TGSA", dir, dir_list)
Pred_TGDRP_SG = read_pred_all(dir_list_tgdrp, "TGDRP_SG", pattern)
Pred_TGSA_SG = read_pred_all(dir_list_tgsa, "TGSA_SG", pattern)

# DRPreter & SA [O/O]
dir = "../DRPreter/Results"
dir_list_drp = sprintf("%s/%s/DRPreter", dir, dir_list)
dir_list_drpsa = sprintf("%s/%s/DRPreter_SA", dir, dir_list)
Pred_DRPreter = read_pred_all(dir_list_drp, "DRPreter", pattern)
Pred_DRPreter_SA = read_pred_all(dir_list_drpsa, "DRPreter_SA", pattern)

# DRPreter & SA [O/△]
dir = "../DRPreter_SANGER/Results"
dir_list_drp = sprintf("%s/%s/DRPreter", dir, dir_list)
dir_list_drpsa = sprintf("%s/%s/DRPreter_SA", dir, dir_list)
Pred_DRPreter_SG = read_pred_all(dir_list_drp, "DRPreter_SG", pattern)
Pred_DRPreter_SA_SG = read_pred_all(dir_list_drpsa, "DRPreter_SA_SG", pattern)

# GCNPath [O]
dir = "../GCNPath/results"
dir_list_ = sprintf("%s/%s/RGCN", dir, dir_list)
Pred_GCNPath = read_pred_all(dir_list_, "GCNPath", pattern)

# cf. Examine cell lines
cells_sanger = Pred_GCNPath$Perf_Cell$Cell %>% unique      # 972
cells_drpreter = Pred_DRPreter$Perf_Cell$Cell %>% unique   # 696
cells_tgsa = Pred_TGSA$Perf_Cell$Cell %>% unique           # 700
cells_paccmann = Pred_PaccMann$Perf_Cell$Cell %>% unique   # 395
cells_graphdrp = Pred_GraphDRP$Perf_Cell$Cell %>% unique   # 969
cells_hidra = Pred_HiDRA$Perf_Cell$Cell %>% unique         # 947

all(cells_drpreter %in% cells_sanger)   # T
all(cells_tgsa %in% cells_sanger)       # T
all(cells_paccmann %in% cells_sanger)   # F [394/395]
all(cells_graphdrp %in% cells_sanger)   # F [964/969]
all(cells_hidra %in% cells_sanger)      # F [944/947]



##### 3-1. Compare Performances [Overall]

perf_avg_sd = function(Perf_Fold, stat="RMSE", space_cols=F) {
  
  object = function(x, envir=NULL) {
    if (is.null(envir)) envir=parent.frame()
    eval(parse(text=x), envir=envir)
  }
  
  avg_plus_sd = function(x_mean, x_sd) sprintf("%.3f±%.3f", x_mean, x_sd)
  Perf_Fold = Perf_Fold %>% subset(!is.na(object(stat)))
  Perf_Fold_Avg = Perf_Fold %>% reshape2::acast(Model~Test, value.var=stat, fun.aggregate=mean)
  Perf_Fold_SD = Perf_Fold %>% reshape2::acast(Model~Test, value.var=stat, fun.aggregate=sd)
  
  Perf_Fold_Avg = Perf_Fold_Avg %>% as.data.frame
  Perf_Fold_SD = Perf_Fold_SD %>% as.data.frame
  Perf_Fold_Sum = mapply(avg_plus_sd, Perf_Fold_Avg, Perf_Fold_SD)
  if (space_cols) colnames(Perf_Fold_Sum) = gsub(" x ", "\n", colnames(Perf_Fold_Sum))
  
  Perf_Fold_Sum = Perf_Fold_Sum %>% as.data.frame
  rownames(Perf_Fold_Sum) = rownames(Perf_Fold_Avg)
  return(Perf_Fold_Sum)
}

rbind_perf = function(Perf1, Perf2, col=NULL) {
  if (is.null(col)) col = intersect(colnames(Perf1), colnames(Perf2))
  if (is.data.table(Perf1) & is.data.table(Perf2)) {
    Perf = rbind(Perf1[, col, with=F], Perf2[, col, with=F])
  } else Perf = rbind(Perf1[, col], Perf2[, col])
  return(Perf)
}

model_names = c("tCNNS", "HiDRA", "PaccMann", "PaccMann_SG",
                "GraphDRP", "TGDRP", "TGSA", "TGDRP_SG", "TGSA_SG",
                "DRPreter", "DRPreter_SA", "DRPreter_SG", "DRPreter_SA_SG", "GCNPath")

Perf_List = list()
pred_names = sprintf("Pred_%s", model_names)

for (i in 1:length(model_names)) Perf_List[[i]] = get(pred_names[i])$Perf
Perf_List = Reduce(rbind_perf, Perf_List)   # 2310 x 11

dataset = c("GDSC", "GDSC1", "GDSC2")
test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")

Perf_List = Perf_List %>% 
  mutate(Model = Model %>% factor(levels=model_names), 
         Dataset = Dataset %>% factor(levels=dataset), 
         Test_Type = Test_Type %>% factor(levels=test_type))

Test = expand.grid(dataset, test_type) %>% arrange(Var2)
test = sprintf("%s x %s", Test$Var1, Test$Var2)

# RMSE Performances [Mean & SD]
Perf_List$Test = paste(Perf_List$Dataset, Perf_List$Test_Type, sep=" x ")
Perf_List$Test = Perf_List$Test %>% factor(levels=test)

Perf_RMSE_Avg = Perf_List %>% perf_avg_sd(stat="RMSE")
Perf_PCC_Avg = Perf_List %>% perf_avg_sd(stat="PCC")
Perf_SCC_Avg = Perf_List %>% perf_avg_sd(stat="SCC")

# The number of IC50, Cell, Drug
IC50_Number = Perf_List %>% subset(Test_Type=="Normal") %>% 
  group_by(Model, Dataset) %>% summarise(n=sum(N_Test)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

Perf_Cell = list()
Perf_Drug = list()

for (i in 1:length(pred_names)) Perf_Cell[[i]] = get(pred_names[i])$Perf_Cell
for (i in 1:length(pred_names)) Perf_Drug[[i]] = get(pred_names[i])$Perf_Drug

Perf_Cell = Reduce(rbind_perf, Perf_Cell)   # 142028 x 8
Perf_Drug = Reduce(rbind_perf, Perf_Drug)   # 55384 x 8

# Number of cells, drugs
ulen = function(x) x %>% unique %>% length
Perf_Cell$Model = Perf_Cell$Model %>% factor(levels=model_names)
Perf_Drug$Model = Perf_Drug$Model %>% factor(levels=model_names)

Cell_Number = Perf_Cell %>% group_by(Model, Dataset) %>% summarise(n=ulen(Cell)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame
Drug_Number = Perf_Drug %>% group_by(Model, Dataset) %>% summarise(n=ulen(Drug)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

colnames(Cell_Number) = sprintf("Cell_%s", colnames(Cell_Number))
colnames(Drug_Number) = sprintf("Drug_%s", colnames(Drug_Number))
colnames(IC50_Number) = sprintf("IC50_%s", colnames(IC50_Number))

identical(rownames(Cell_Number), rownames(Drug_Number))   # T
identical(rownames(Cell_Number), rownames(IC50_Number))   # T
Stat_Number = Reduce(cbind, list(Cell_Number, Drug_Number, IC50_Number))
Stat_Number = Stat_Number[, c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

file = "Performance_Summary.xlsx"
sheet = c("RMSE_Test", "PCC_Test", "SCC_Test", "Number")
df_list = list(Perf_RMSE_Avg, Perf_PCC_Avg, Perf_SCC_Avg, Stat_Number)
write.xlsx(df_list, file=file, rowNames=T, sheetName=sheet)


# Prediction
save_pred = T
if (save_pred) {
  Pred_List = list()
  for (i in 1:length(model_names)) Pred_List[[i]] = get(pred_names[i])$Pred
  Pred_List = Reduce(rbind_perf, Pred_List)   # 40848552 x 8
  Pred_List = Pred_List %>% rename(Train_Fold=Fold)
  
  Pred_List_ = rbind(Pred_tCNNS$Pred, Pred_GraphDRP$Pred)
  Pred_List_ = Pred_List_ %>% rename(Train_Fold=Fold) %>% 
    subset(select=-c(Cell_BROAD, Cell_COSMIC))
  
  file = "Prediction [GDSC].csv"
  fwrite(Pred_List, file=file, row.names=F)
  
  file = "Prediction [GDSC (tCNNS, GraphDRP)].csv"
  fwrite(Pred_List_, file=file, row.names=F)
  
  rm(Pred_List, Pred_List_)
  gc()
}


# Scatter plots of prediction-actual ln(IC50) values
# Be careful that it takes too much times to draw all figures...

plot_pred = function(Pred, model=NULL, test_type="Normal", dir=NULL, 
                     width=15, height=15, return_info=T) {
  
  xlab = bquote(Actual~ln(IC[50]))
  ylab = bquote(Predicted~ln(IC[50]))
  main = sprintf("%s/Prediction [%s, %s]", dir, model, test_type)
  
  Pred_ = Pred %>% subset(!is.na(Prediction) & !is.infinite(Prediction))
  rmse = Pred_ %>% with(RMSE(LN_IC50, Prediction)) %>% round(3)
  corr = Pred_ %>% with(cor(LN_IC50, Prediction)) %>% round(3)
  info = c(nrow(Pred), rmse, corr)
  sprintf("# N=%s, RMSE=%.3f, PCC=%.3f", info[1], info[2], info[3]) %>% print
  if (nrow(Pred_)!=0) sprintf("# Inf found... [%s, n=%s]", model, nrow(Pred)-nrow(Pred_)) %>% print
  
  Pred_ %>% plot_def(LN_IC50, Prediction, main=main, xlab=xlab, ylab=ylab,
                     size=1.5, alpha=0.25, axis_tl=30, axis_tx=24, dpi=1200,
                     width=width, height=height, xy_line=T, raster=T, save=T, save_svg=T)
  
  if (return_info) return(info)
}

draw_pred = T
if (draw_pred) {
  Info_N = data.frame()
  Info_C = data.frame()
  Info_D = data.frame()
  Info_S = data.frame()
  
  dir_n = mkdir("Prediction/Normal")
  dir_c = mkdir("Prediction/Cell_Blind")
  dir_d = mkdir("Prediction/Drug_Blind")
  dir_s = mkdir("Prediction/Strict_Blind")
  
  for (i in 1:length(pred_names)) {
    Pred_Temp = get(pred_names[i])$Pred
    model = gsub("Pred_", "", pred_names[i])
    
    info_n = Pred_Temp %>% 
      subset(Test_Type=="Normal" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_n, model=model, test_type="Normal")
    info_c = Pred_Temp %>% 
      subset(Test_Type=="Cell_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_c, model=model, test_type="Cell-Blind")
    info_d = Pred_Temp %>% 
      subset(Test_Type=="Drug_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_d, model=model, test_type="Drug-Blind")
    info_s = Pred_Temp %>% 
      subset(Test_Type=="Strict_Blind" & Dataset=="GDSC") %>% 
      plot_pred(dir=dir_s, model=model, test_type="Strict-Blind")
    
    Info_N = Info_N %>% rbind(c(model, info_n))
    Info_C = Info_C %>% rbind(c(model, info_c))
    Info_D = Info_D %>% rbind(c(model, info_d))
    Info_S = Info_S %>% rbind(c(model, info_s))
  }
  
  col = c("Model", "N_Test", "RMSE", "PCC")
  Info_N = Info_N %>% setNames(col)
  Info_C = Info_C %>% setNames(col)
  Info_D = Info_D %>% setNames(col)
  Info_S = Info_S %>% setNames(col)
}


# Statistic Analysis per Dataset & Test-Type
# Non-parametric U-Test [Mann-Whitney]

wilcox_test_sub = function(Perf_List, model_h1, model_h0, dataset=NULL, test_type=NULL) {
  
  if (!is.null(dataset)) Perf_List = Perf_List %>% subset(Dataset %in% dataset)
  if (!is.null(test_type)) Perf_List = Perf_List %>% subset(Test_Type %in% test_type)
  Perf_Model_H1 = Perf_List %>% subset(Model==model_h1)
  Perf_Model_H0 = Perf_List %>% subset(Model==model_h0)
  
  utest_rmse = wilcox.test(Perf_Model_H1$RMSE, Perf_Model_H0$RMSE, alternative="less")
  utest_pcc = wilcox.test(Perf_Model_H1$PCC, Perf_Model_H0$PCC, alternative="greater")
  utest_scc = wilcox.test(Perf_Model_H1$SCC, Perf_Model_H0$SCC, alternative="greater")
  dataset = ifelse(!is.null(dataset), dataset, "Total")
  test_type = ifelse(!is.null(test_type), test_type, "Total")
  
  utest_res = c(model_h0, model_h1, dataset, test_type, 
                utest_rmse$alternative, utest_pcc$alternative, utest_scc$alternative, 
                utest_rmse$statistic, utest_pcc$statistic, utest_scc$statistic,
                utest_rmse$p.value, utest_pcc$p.value, utest_scc$p.value)
  
  return(utest_res)
}

wilcox_test = function(Perf_List, fdr_adjust=T) {
  
  Test_Res = data.frame()
  dataset = c("GDSC", "GDSC1", "GDSC2")
  test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
  
  model_list = Perf_List$Model %>% unique
  Model_Pair = expand.grid(Model_H1=model_list, Model_H0=model_list)
  Model_Pair = Model_Pair %>% subset(Model_H1!=Model_H0)
  
  for (i in 1:nrow(Model_Pair)) {
    model_h1 = Model_Pair$Model_H1[i] %>% as.character
    model_h0 = Model_Pair$Model_H0[i] %>% as.character
    
    for (dataset_ in dataset) {
      for (test_type_ in test_type) {
        test_res = wilcox_test_sub(Perf_List, model_h1, model_h0, dataset=dataset_, test_type=test_type_)
        Test_Res = Test_Res %>% rbind(test_res)
      }
    }
    
    for (test_type_ in test_type) {
      test_res = wilcox_test_sub(Perf_List, model_h1, model_h0, test_type=test_type_)
      Test_Res = Test_Res %>% rbind(test_res)
    }
    
    for (dataset_ in dataset) {
      test_res = wilcox_test_sub(Perf_List, model_h1, model_h0, dataset=dataset_)
      Test_Res = Test_Res %>% rbind(test_res)
    }
    
    test_res = wilcox_test_sub(Perf_List, model_h1, model_h0)
    Test_Res = Test_Res %>% rbind(test_res)
  }
  
  colnames(Test_Res) = c("Model_H0", "Model_H1", "Dataset", "Test_Type", 
                         "Hypothesis_RMSE", "Hypothesis_PCC", "Hypothesis_SCC", 
                         "W_RMSE", "W_PCC", "W_SCC", "Pval_RMSE", "Pval_PCC", "Pval_SCC")
  
  Test_Res[, 8:13] = Test_Res[, 8:13] %>% sapply(as.numeric)
  fdr_adj = function(x) p.adjust(x, "fdr")
  mlog10 = function(x) -log10(x)
  signif = function(x) x<0.05
  
  metrics = c("RMSE", "SCC", "PCC")
  col_fdr = sprintf("FDR_%s", metrics)
  col_pval = sprintf("Pval_%s", metrics)
  
  if (fdr_adjust) {
    Test_Res = Test_Res %>% 
      group_by(Dataset, Test_Type) %>%
      mutate(across(all_of(col_pval), fdr_adj, .names="{sub('Pval', 'FDR', col)}")) %>%
      mutate(across(all_of(col_fdr), mlog10, .names = "MLog10_{col}")) %>%
      mutate(across(all_of(col_fdr), signif, .names = "Signif_{col}"))
  } else {
    Test_Res = Test_Res %>% 
      mutate(across(all_of(col_pval), mlog10, .names = "MLog10_{col}")) %>%
      mutate(across(all_of(col_pval), signif, .names = "Signif_{col}"))
  }
  
  Test_Res = Test_Res %>% as.data.frame
  return(Test_Res)
}

wilcox_model_grid = function(Perf_UTest, fdr_adjust=T, models=NULL, main=NULL, lvl_model=NULL,
                             axis_tx=18, legend_tl=16, legend_tx=16, size=4.2, angle=36,
                             axis_face="plain", legend_face="plain", text_face="plain", 
                             round=2, width=30.6, height=24, asterisk=T, save=F) {
  
  if (!is.null(models)) {
    Perf_UTest = Perf_UTest %>% 
      subset(Model_H0 %in% models & Model_H1 %in% models)
  }
  
  if (!is.null(lvl_model)) {
    Perf_UTest = Perf_UTest %>% 
      mutate(Model_H0=factor(Model_H0, levels=lvl_model), 
             Model_H1=factor(Model_H1, levels=lvl_model))
  }
  
  if (fdr_adjust) {
    legend1 = bquote(atop("-Log"["10"](FDR), p(RMSE["y"]<RMSE["x"])))
    legend2 = bquote(atop("-Log"["10"](FDR), p(PCC["y"]>PCC["x"])))
    legend3 = bquote(atop("-Log"["10"](FDR), p(SCC["y"]>SCC["x"])))
    
    Perf_UTest = Perf_UTest %>% 
      mutate(Pval_RMSE_Log=-log(FDR_RMSE, 10), 
             Pval_PCC_Log=-log(FDR_PCC, 10), 
             Pval_SCC_Log=-log(FDR_SCC, 10))
  } else {
    legend1 = bquote(atop("-Log"["10"](Pval), p(RMSE["y"]<RMSE["x"])))
    legend2 = bquote(atop("-Log"["10"](Pval), p(PCC["y"]>PCC["x"])))
    legend3 = bquote(atop("-Log"["10"](Pval), p(SCC["y"]>SCC["x"])))
    
    Perf_UTest = Perf_UTest %>% 
      mutate(Pval_RMSE_Log=-log(Pval_RMSE, 10), 
             Pval_PCC_Log=-log(Pval_PCC, 10), 
             Pval_SCC_Log=-log(Pval_SCC, 10))
  }
  
  if (asterisk) {
    if (fdr_adjust) {
      idx1 = Perf_UTest$FDR_RMSE<0.05
      idx2 = Perf_UTest$FDR_PCC<0.05
      idx3 = Perf_UTest$FDR_SCC<0.05
    } else {
      idx1 = Perf_UTest$Pval_RMSE<0.05
      idx2 = Perf_UTest$Pval_PCC<0.05
      idx3 = Perf_UTest$Pval_SCC<0.05
    }
    
    Perf_UTest = Perf_UTest %>% 
      mutate(Pval_RMSE_Log_=round(Pval_RMSE_Log, round), 
             Pval_PCC_Log_=round(Pval_PCC_Log, round), 
             Pval_SCC_Log_=round(Pval_SCC_Log, round))
    
    Perf_UTest$Pval_RMSE_Log_[idx1] = paste0(Perf_UTest$Pval_RMSE_Log_[idx1], "*")
    Perf_UTest$Pval_PCC_Log_[idx2] = paste0(Perf_UTest$Pval_PCC_Log_[idx2], "*")
    Perf_UTest$Pval_SCC_Log_[idx3] = paste0(Perf_UTest$Pval_SCC_Log_[idx3], "*")
    
    label1 = aes(label=Pval_RMSE_Log_)
    label2 = aes(label=Pval_PCC_Log_)
    label3 = aes(label=Pval_SCC_Log_)
  }
  
  Perf_UTest = Perf_UTest %>% as.data.frame
  color = scale_fill_gradient(low="white", high="firebrick3")
  add = list(theme(legend.key.width=unit(0.8, "cm"), 
                   legend.key.height=unit(1, "cm")))
  
  Perf_UTest %>% grid_def(Model_H0, Model_H1, fill=Pval_RMSE_Log, main=main[1], color=color,
                          legend=legend1, label=label1, round=round, size=size, angle=angle,
                          axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx, margin_lg=0.5,
                          axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                          add=add, width=width+1.2, height=height, mean_summary=F, save=save)
  
  Perf_UTest %>% grid_def(Model_H0, Model_H1, fill=Pval_PCC_Log, main=main[2], color=color,
                          legend=legend2, label=label2, round=round, size=size, angle=angle,
                          axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx, margin_lg=0.5,
                          axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                          add=add, width=width+0.5, height=height, mean_summary=F, save=save)
  
  Perf_UTest %>% grid_def(Model_H0, Model_H1, fill=Pval_SCC_Log, main=main[3], color=color,
                          legend=legend3, label=label3, round=round, size=size, angle=angle,
                          axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx, margin_lg=0.5,
                          axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                          add=add, width=width+0.5, height=height, mean_summary=F, save=save)
}

Perf_UTest = Perf_List %>% wilcox_test(fdr_adjust=T)

# Confirmation of FDR adjust p.val
utest_is_fdr = c()
for (dataset_ in dataset) {
  for (test_type_ in test_type) {
    Temp = Perf_UTest %>% subset(Dataset==dataset_ & Test_Type==test_type_)
    utest_is_fdr_ = identical(Temp$FDR_RMSE, p.adjust(Temp$Pval_RMSE, "fdr"))
    utest_is_fdr = utest_is_fdr %>% c(utest_is_fdr_)
  }
}

all(utest_is_fdr)   # T
test_type_ = gsub("_", "-", test_type)
Perf_UTest_N = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_UTest_C = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_UTest_D = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_UTest_S = Perf_UTest %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])

dir = mkdir("Performance [Wilcox Test, Grid]")
main1 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[1])
main2 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[2])
main3 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[3])
main4 = sprintf("%s/%s Wilcox of Models [GDSC, %s]", dir, metrics, test_type_[4])

Perf_UTest_N %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main1, save=T)
Perf_UTest_C %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main2, save=T)
Perf_UTest_D %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main3, save=T)
Perf_UTest_S %>% wilcox_model_grid(models=model_names, lvl_model=model_names, main=main4, save=T)

plot_perf = function(Perf_List, score="RMSE", test_type="Normal", 
                     dir=NULL, axis_tl=36, axis_tx=24, legend_tl=25, legend_tx=25, alpha=0.9, 
                     margin=0.6, margin_lg=0.5, pos_dodge=0.75, width=36, height=22.5, save=F) {
  
  # [score] RMSE, PCC, SCC
  pos = position_dodge(pos_dodge)
  test_type_ = gsub("_", "-", test_type)
  if (test_type_=="Normal") test_type_ = "Unblinded"
  ylab = sprintf("%s [%s Test]", score, test_type_)
  file = sprintf("%s/Performances [%s, %s]", dir, score, test_type_)
  
  Perf_List = Perf_List %>% subset(Test_Type==test_type)
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
  # add = list(scale_fill_brewer(palette="Reds"))
  # add = list(scale_fill_manual(values=c("#610000", "#ff1e1e", "#feb9b9")))
  
  theme_lg = theme(legend.key.size=unit(1.6, 'cm'))
  add = list(scale_fill_manual(values=color), theme_lg)
  
  Perf_List %>% boxplot_def(Model, object(score), fill=Dataset,
                            main=file, legend="Dataset", ylab=ylab, pos=pos,
                            axis_tl=axis_tl, axis_tx=axis_tx, vjust=1, hjust=1, add=add,
                            legend_tl=legend_tl, legend_tx=legend_tx, margin=margin, pos_legend="bottom",
                            margin_lg=margin_lg, alpha=alpha, width=width, height=height, save=save)
}

plot_perf_bar = function(Perf_List, score="RMSE", dataset=NULL, dir=NULL,
                         test_type="Normal", col_highlight="#206ba3", 
                         width=18, height=13.5, axis_tl=27, axis_tx=13.5,  
                         reorder_score=T, model_except=NULL, include_null=F, save=T) {
  
  # width=30, height=18
  suppressMessages(library(ggtext))
  
  highlight = function(x, pattern, color="black", family="") {
    suppressMessages(library(glue))
    suppressMessages(library(ggtext))
    ifelse(grepl(pattern, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
  }
  
  font1 = font("ylab", size=axis_tl, margin=margin(r=0.25, unit="cm"))
  font2 = font("x.text", color="grey30", size=axis_tx, margin=margin(t=0.25, unit="cm"))
  font3 = font("y.text", color="grey30", size=axis_tx, margin=margin(r=0.25, unit="cm"))
  
  pos = position_dodge(0.8)
  font = font1 + font2 + font3
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
  main = sprintf("%s/Performance [%s, %s]", dir, score, test_type)
  
  Perf_List = Perf_List %>% subset(Test_Type==test_type)
  # levels(Perf_List$Model)[levels(Perf_List$Model)=="GCNPath"] = "<b>GCNPath</b>"
  
  if (!is.null(dataset)) {
    Perf_List = Perf_List %>% subset(Dataset==dataset)
    main = sprintf("%s/Performance [%s, %s, %s]", dir, score, test_type, dataset)
    if (include_null) main = sprintf("%s/Performance [%s, %s, %s (+Null Model)]", dir, score, test_type, dataset)
  }
  
  if (!is.null(model_except)) {
    Perf_List = Perf_List %>% subset(!(Model %in% model_except))
    Perf_List$Model = Perf_List$Model %>% base::droplevels()
    if (is.null(dataset)) {
      main = sprintf("%s/Performance [%s, %s (Except)]", dir, score, test_type)
    } else main = sprintf("%s/Performance [%s, %s, %s (Except)]", dir, score, test_type, dataset)
  }
  
  if (reorder_score) {
    lvl_models = Perf_List %>% group_by(Model) %>%
      summarise(Mean=mean(object(score), na.rm=T)) %>%
      arrange(desc(Mean)) %>% pull(Model)
    if (!(score %in% c("RMSE"))) {
      lvl_models = lvl_models %>% rev
    }
    Perf_List$Model = Perf_List$Model %>% factor(levels=lvl_models)
  }
  
  yval = Perf_List[[score]]
  highlight_x = function(x) highlight(x, "GCNPath", color=col_highlight, family="bold")
  # sd_max = Perf_List %>% group_by(Model) %>%
  #   summarise(Mean=mean(object(score), na.rm=T),
  #             SD=sd(object(score)), na.rm=T) %>%
  #   filter(Mean==max(Mean)) %>% pull(SD)
  
  if (score %in% c("RMSE")) {
    ymin = min(yval, na.rm=T)-0.002
    ymax = max(yval, na.rm=T)+0.002
  } else {
    ymin = max(min(yval, na.rm=T)-0.002, -1)
    ymax = min(max(yval, na.rm=T)+0.002, 1)
  }
  
  if (is.null(dataset)) {
    pl = Perf_List %>%
      ggbarplot(x="Model", y=score, color="black", fill="Dataset",
                xlab=F, add="mean_se", position=pos) +
      scale_fill_manual(values=color)
  } else {
    pl = Perf_List %>%
      ggbarplot(x="Model", y=score, color="black", fill="Model",
                xlab=F, add="mean_se", position=pos)
    
    models = Perf_List$Model %>% unique %>% as.character
    fill = c("royalblue3", rep("#66b3ed", length(models)-1))
    names(fill) = c("GCNPath", models[models!="GCNPath"])
    pl = pl + scale_fill_manual(values=fill)
    pl = pl %>% ggpar(legend="none")
    # cf. #539ed6
  }
  
  pl = pl + font + rotate_x_text(30)
  pl = pl %>% ggpar(ylim=c(ymin, ymax))
  pl = pl + geom_point(alpha=0.5) +
    scale_x_discrete(labels=highlight_x) + 
    theme(axis.text.x=element_markdown())
  # theme(axis.text.x=element_markdown(size=15, color="grey30", margin=margin(t=10, unit="pt")))
  
  if (save) {
    pl %>% save_fig_ggpubr(main=main, width=width, height=height, svg=T)
  } else print(pl)
}

dir = mkdir("Performance [Overall, Boxplot]")

Perf_List %>% plot_perf("RMSE", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("RMSE", "Strict_Blind", dir=dir, save=T)

Perf_List %>% plot_perf("PCC", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("PCC", "Strict_Blind", dir=dir, save=T)

Perf_List %>% plot_perf("SCC", "Normal", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Cell_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Drug_Blind", dir=dir, save=T)
Perf_List %>% plot_perf("SCC", "Strict_Blind", dir=dir, save=T)

dir = mkdir("Performance [Overall, Boxplot (Except)]")
models_except = c("tCNNS", "HiDRA", "GraphDRP", "PaccMann", "PaccMann_SG")

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("RMSE", "Normal", dir=dir, width=30, save=T)

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("PCC", "Normal", dir=dir, width=36, save=T)

Perf_List %>% 
  subset(!(Model %in% models_except)) %>% 
  mutate(Model=droplevels(Model)) %>%
  plot_perf("SCC", "Normal", dir=dir, width=36, save=T)


# Barplot [GDSC1+2 Test]
dir = mkdir("Performance [Overall, Barplot]")
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[1], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[2], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[3], dir=dir, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[4], dir=dir, save=T)

model_except_n = c("tCNNS", "HiDRA", "PaccMann", "PaccMann_SG", "GraphDRP")
model_except_c = c("tCNNS", "GraphDRP")
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[1],
                            model_except=model_except_n, dir=dir, width=20, axis_tx=16.5, save=T)
Perf_List %>% plot_perf_bar(score="RMSE", dataset="GDSC", test_type=test_type[2],
                            model_except=model_except_c, dir=dir, width=20, axis_tx=16.5, save=T)



##### 3-2. Compare Performances [Cell & Drug]

compare_model_grid = function(Perf, models=NULL, by="Cell", main=NULL, lvl_model=NULL,
                              axis_tx=16.5, legend_tl=18, legend_tx=18, size=4.5,
                              axis_face="plain", legend_face="plain", text_face="plain", 
                              round=2, width=30, height=24, save=F) {
  
  # axis_tx=20, legend_tl=18, legend_tx=16.5, size=5.25
  # round=2, width=36, height=30
  
  if (!is.null(models)) {
    Perf = Perf %>% subset(Model %in% models)
  }
  
  if (by=="Cell") {
    Perf_RMSE = Perf %>% reshape2::acast(Cell~Model, value.var="RMSE")
    Perf_PCC = Perf %>% reshape2::acast(Cell~Model, value.var="PCC")
    Perf_SCC = Perf %>% reshape2::acast(Cell~Model, value.var="SCC")
  } else {
    Perf_RMSE = Perf %>% reshape2::acast(Drug~Model, value.var="RMSE")
    Perf_PCC = Perf %>% reshape2::acast(Drug~Model, value.var="PCC")
    Perf_SCC = Perf %>% reshape2::acast(Drug~Model, value.var="SCC")
  }
  
  Perf_RMSE_Corr = Perf_RMSE %>% cor(use="pairwise.complete.obs")
  Perf_PCC_Corr = Perf_PCC %>% cor(use="pairwise.complete.obs")
  Perf_SCC_Corr = Perf_SCC %>% cor(use="pairwise.complete.obs")
  
  Perf_RMSE_Corr = Perf_RMSE_Corr %>% as.matrix %>% reshape2::melt() %>% as.data.frame
  Perf_PCC_Corr = Perf_PCC_Corr %>% as.matrix %>% reshape2::melt() %>% as.data.frame
  Perf_SCC_Corr = Perf_SCC_Corr %>% as.matrix %>% reshape2::melt() %>% as.data.frame
  
  if (!is.null(lvl_model)) {
    by_ = c("Var1"="Var1", "Var2"="Var2")
    full_join_ = function(df1, df2) full_join(df1, df2, by=by_)
    Perf_Corr = Reduce(full_join_, list(Perf_RMSE_Corr, Perf_PCC_Corr, Perf_SCC_Corr))
    Perf_Corr = Perf_Corr %>% mutate(Var1=factor(Var1, levels=lvl_model),
                                     Var2=factor(Var2, levels=lvl_model))
  }
  
  colnames(Perf_Corr) = c("Model1", "Model2", "RMSE_Corr", "PCC_Corr", "SCC_Corr")
  color = scale_fill_gradient(low="white", high="firebrick3")
  add = list(theme(legend.key.width=unit(0.8, "cm"), 
                   legend.key.height=unit(1, "cm")))
  
  legend_subs = if (by=="Cell") "C" else "D"
  legend1 = bquote(PCC(RMSE[.(legend_subs)]))
  legend2 = bquote(PCC(PCC[.(legend_subs)]))
  legend3 = bquote(PCC(SCC[.(legend_subs)]))
  
  Perf_Corr %>% grid_def(Model1, Model2, fill=RMSE_Corr, main=main[1], add=add,
                         legend=legend1, round=round, size=size, color=color, 
                         axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx,
                         axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                         margin_lg=0.5, width=width, height=height, mean_summary=F, save=save)
  
  Perf_Corr %>% grid_def(Model1, Model2, fill=PCC_Corr, main=main[2], add=add,
                         legend=legend2, round=round, size=size, color=color, 
                         axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx,
                         axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                         margin_lg=0.5, width=width, height=height, mean_summary=F, save=save)
  
  Perf_Corr %>% grid_def(Model1, Model2, fill=SCC_Corr, main=main[3], add=add,
                         legend=legend3, round=round, size=size, color=color, 
                         axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx,
                         axis_face=axis_face, legend_face=legend_face, text_face=text_face,
                         margin_lg=0.5, width=width, height=height, mean_summary=F, save=save)
  
  return(Perf_Corr)
}

# Cell Annotation from SANGER Cell Passports
file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)

# Drug Annotation from GDSC
file = "Anno_Drugs.csv"
Anno_Drugs = read.csv(file)

idx = match(Perf_Cell$Cell, Anno_Cells$SANGER_MODEL_ID)
Perf_Cell = Perf_Cell %>% mutate(Cell_TCGA=Anno_Cells$TCGA_CODE[idx]) %>% relocate(Cell_TCGA, .after=Cell)
Perf_Cell$Cell_TCGA %>% is.na %>% sum   # 0
Perf_Cell$Cell_TCGA[is.na(Perf_Cell$Cell_TCGA)] = "UNCLASSIFIED"

idx = match(Perf_Drug$Drug, Anno_Drugs$Drug_CID)
Perf_Drug = Perf_Drug %>% 
  mutate(Drug_Pathway=Anno_Drugs$Target_Pathway[idx]) %>% 
  relocate(Drug_Pathway, .after=Drug) %>% as.data.frame

Perf_Drug$Drug_Pathway %>% is.na %>% sum   # 224
Perf_Drug$Drug_Pathway[is.na(Perf_Drug$Drug_Pathway)] = "Unclassified"

Perf_Cell_N = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_Cell_C = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_Cell_D = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_Cell_S = Perf_Cell %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])

Perf_Drug_N = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[1])
Perf_Drug_C = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[2])
Perf_Drug_D = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[3])
Perf_Drug_S = Perf_Drug %>% subset(Dataset=="GDSC" & Test_Type==test_type[4])


dir = mkdir("Performance [Cell and Drug, Grid]")
metrics = c("RMSE", "PCC", "SCC")

test_type_ = gsub("_", "-", test_type)
file1 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[1])
file2 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[2])
file3 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[3])
file4 = sprintf("%s/%s Corr of Models [Cell, %s]", dir, metrics, test_type_[4])

Perf_Cell_Corr_N = Perf_Cell_N %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file1, save=T)
Perf_Cell_Corr_C = Perf_Cell_C %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file2, save=T)
Perf_Cell_Corr_D = Perf_Cell_D %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file3, save=T)
Perf_Cell_Corr_S = Perf_Cell_S %>% 
  compare_model_grid(models=model_names, by="Cell", lvl_model=model_names, main=file4, save=T)


file1 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[1])
file2 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[2])
file3 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[3])
file4 = sprintf("%s/%s Corr of Models [Drug, %s]", dir, metrics, test_type_[4])

Perf_Drug_Corr_N = Perf_Drug_N %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file1, save=T)
Perf_Drug_Corr_C = Perf_Drug_C %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file2, save=T)
Perf_Drug_Corr_D = Perf_Drug_D %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file3, save=T)
Perf_Drug_Corr_S = Perf_Drug_S %>% 
  compare_model_grid(models=model_names, by="Drug", lvl_model=model_names, main=file4, save=T)

# Number of cells, drugs
ulen = function(x) x %>% unique %>% length
Perf_Cell$Model = Perf_Cell$Model %>% factor(levels=model_names)
Perf_Drug$Model = Perf_Drug$Model %>% factor(levels=model_names)

Cell_Number = Perf_Cell %>% group_by(Model, Dataset) %>% summarise(n=ulen(Cell)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame
Drug_Number = Perf_Drug %>% group_by(Model, Dataset) %>% summarise(n=ulen(Drug)) %>% 
  as.data.frame %>% acast(Model~Dataset, value.var="n") %>% as.data.frame

colnames(Cell_Number) = sprintf("Cell_%s", colnames(Cell_Number))
colnames(Drug_Number) = sprintf("Drug_%s", colnames(Drug_Number))
colnames(IC50_Number) = sprintf("IC50_%s", colnames(IC50_Number))

identical(rownames(Cell_Number), rownames(Drug_Number))   # T
identical(rownames(Cell_Number), rownames(IC50_Number))   # T
Stat_Number = Reduce(cbind, list(Cell_Number, Drug_Number, IC50_Number))
Stat_Number = Stat_Number[, c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

Perf_Cell_Corr = list(Perf_Cell_Corr_N, Perf_Cell_Corr_C, Perf_Cell_Corr_D, Perf_Cell_Corr_S)
Perf_Drug_Corr = list(Perf_Drug_Corr_N, Perf_Drug_Corr_C, Perf_Drug_Corr_D, Perf_Drug_Corr_S)

names(Perf_Cell_Corr) = test_type
names(Perf_Drug_Corr) = test_type

file = sprintf("%s/Performance_Cell.csv", dir)
write.csv(Perf_Cell, file=file, row.names=F)

file = sprintf("%s/Performance_Drug.csv", dir)
write.csv(Perf_Drug, file=file, row.names=F)



##### Performance Examination [Detailed]
### Which factors affect the RMSE (based on GCNPath)?

test = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
Pred_GCNPath_N = Pred_GCNPath$Pred %>% as.data.frame %>% 
  subset(Test_Type==test[1] & Dataset=="GDSC")

# Anno_Cells, Anno_Drugs
idx1 = match(Pred_GCNPath_N$Cell, Anno_Cells$SANGER_MODEL_ID)
idx2 = match(Pred_GCNPath_N$Drug, Anno_Drugs$Drug_CID)

Pred_GCNPath_N = Pred_GCNPath_N %>% 
  mutate(Cell_TCGA=Anno_Cells$TCGA_CODE[idx1], 
         Drug_Pathway=Anno_Drugs$Target_Pathway[idx2])

Pred_GCNPath_N$Cell_TCGA %>% is.na %>% sum      # 0
Pred_GCNPath_N$Drug_Pathway %>% is.na %>% sum   # 1426
Pred_GCNPath_N$Drug_Pathway[is.na(Pred_GCNPath_N$Drug_Pathway)] = "Unclassified"



### 1. IC50 Mean & Variation

dir = mkdir("Analysis of Error")
col = c("Cell", "Drug", "LN_IC50", "Prediction", "Cell_TCGA", "Drug_Pathway")

IC50_Stat_Cell = Pred_GCNPath_N[, col] %>% group_by(Cell, Cell_TCGA) %>% 
  summarise(Num=n(), RMSE=RMSE(LN_IC50, Prediction), 
            IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50)) %>% as.data.frame   # 972

IC50_Stat_Drug = Pred_GCNPath_N[, col] %>% group_by(Drug, Drug_Pathway) %>% 
  summarise(Num=n(), RMSE=RMSE(LN_IC50, Prediction), 
            IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50)) %>% as.data.frame   # 432

IC50_Stat_Cell_Group = Pred_GCNPath_N[, col] %>% group_by(Cell_TCGA) %>% 
  summarise(Num=n(), RMSE=RMSE(LN_IC50, Prediction), 
            IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50)) %>% as.data.frame   # 38

IC50_Stat_Drug_Group = Pred_GCNPath_N[, col] %>% group_by(Drug_Pathway) %>% 
  summarise(Num=n(), RMSE=RMSE(LN_IC50, Prediction), 
            IC50_Mean=mean(LN_IC50), IC50_SD=sd(LN_IC50)) %>% as.data.frame   # 24

# by = c("Cell", "Drug")
# ylab = sprintf("RMSE [%s]", by)

ylab_c = bquote(RMSE[C])
ylab_d = bquote(RMSE[D])
xlab = bquote(Mean~ln(IC[50]))
main = sprintf("%s/IC50_Mean & RMSE Relation [%s]", dir, by)
add = list(stat_smooth(method="lm"))

IC50_Stat_Cell %>% with(cor(IC50_Mean, RMSE)) %>% round(3)   # PCC=0.220
IC50_Stat_Cell %>% plot_def(IC50_Mean, RMSE, main=main[1], xlab=xlab, ylab=ylab_c, add=add, 
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

IC50_Stat_Drug %>% with(cor(IC50_Mean, RMSE)) %>% round(3)   # PCC=-0.305
IC50_Stat_Drug %>% plot_def(IC50_Mean, RMSE, main=main[2], xlab=xlab, ylab=ylab_d, add=add, 
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

xlab = bquote(SD~ln(IC[50]))
main = sprintf("%s/IC50_SD & RMSE Relation [%s]", dir, by)

IC50_Stat_Cell %>% with(cor(IC50_SD, RMSE)) %>% round(3)   # PCC=-0.003
IC50_Stat_Cell %>% plot_def(IC50_SD, RMSE, main=main[1], xlab=xlab, ylab=ylab_c, add=add, 
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

IC50_Stat_Drug %>% with(cor(IC50_SD, RMSE)) %>% round(3)   # PCC=0.776
IC50_Stat_Drug %>% plot_def(IC50_SD, RMSE, main=main[2], xlab=xlab, ylab=ylab_d, add=add, 
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)


### 2. Num of Cell & Drug 

subtitle = c("Cell", "Drug", "Cell TCGA", "Drug Pathway")
main = sprintf("%s/IC50 Number & RMSE [%s]", dir, subtitle)

xlab = bquote(Number~of~ln(IC[50]))
# ylab = c("RMSE [Cell]", "RMSE [Drug]", "RMSE [Cell, TCGA Code]", "RMSE [Drug, Target Pathway]")

ylab_c = bquote(RMSE[C])
ylab_d = bquote(RMSE[D])
ylab_c_ = "RMSE per Cancer Type"
ylab_d_ = "RMSE per Target Pathway"

IC50_Stat_Cell %>% with(cor(Num, RMSE)) %>% round(3)   # PCC=-0.198
IC50_Stat_Cell %>% plot_def(Num, RMSE, main=main[1], xlab=xlab, ylab=ylab_c,
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

IC50_Stat_Drug %>% with(cor(Num, RMSE)) %>% round(3)   # PCC=-0.115
IC50_Stat_Drug %>% plot_def(Num, RMSE, main=main[2], xlab=xlab, ylab=ylab_d,
                            size=2, alpha=0.5, axis_tl=30, axis_tx=24, save=T)

add = list(theme(axis.title.x=element_text(size=30)))

IC50_Stat_Cell_Group %>% with(cor(Num, RMSE)) %>% round(3)   # PCC=-0.11
IC50_Stat_Cell_Group %>% plot_def(Num, RMSE, main=main[3], xlab=xlab, ylab=ylab_c_,
                                  size=2, alpha=0.5, axis_tl=25, axis_tx=18, add=add, save=T)

IC50_Stat_Drug_Group %>% with(cor(Num, RMSE)) %>% round(3)   # PCC=0.067
IC50_Stat_Drug_Group %>% plot_def(Num, RMSE, main=main[4], xlab=xlab, ylab=ylab_d_,
                                  size=2, alpha=0.5, axis_tl=25, axis_tx=18, add=add, save=T)

# Cell 972, Drug 432
# Large sample size does not always guarantee the low RMSE, especially drugs...



### 3. IC50 Variation [GDSC-CCLE]

# IC50 CCLE
# GDSC_Last/processed_data/ic50_data/IC50_CCLE.csv
IC50_CCLE = read.csv("IC50_CCLE.csv")
IC50_CCLE = IC50_CCLE %>% rename(LN_IC50_CCLE=LN_IC50)
# IC50_CCLE, IC50_CCLE_FT

col = c("Cell_SANGER_ID", "Drug_CID", "LN_IC50_CCLE")
by = c("Cell"="Cell_SANGER_ID", "Drug"="Drug_CID")
Pred_GCNPath_N_CCLE = Pred_GCNPath_N %>% inner_join(IC50_CCLE[, col], by=by)   # 5567

min_ccle = min(Pred_GCNPath_N_CCLE$LN_IC50_CCLE)
max_ccle = max(Pred_GCNPath_N_CCLE$LN_IC50_CCLE)
sum(Pred_GCNPath_N_CCLE$LN_IC50_CCLE>=max_ccle)   # 2967 [53.29%]
sum(Pred_GCNPath_N_CCLE$LN_IC50_CCLE<=min_ccle)   # 38 [0.68%]

Pred_GCNPath_N_CCLE$Cell %>% unique %>% length    # 374
Pred_GCNPath_N_CCLE$Drug %>% unique %>% length    # 18

Pred_GCNPath_N_CCLE$LN_IC50 %>% is.na %>% sum         # 0
Pred_GCNPath_N_CCLE$LN_IC50_CCLE %>% is.na %>% sum    # 0
Pred_GCNPath_N_CCLE$LN_IC50_CCLE %>% hist             # Almost 2.07944154 [Max 8uM]


xlab = bquote(GDSC~ln(IC[50]))
ylab = bquote(CCLE~ln(IC[50]))

main = sprintf("%s/LN_IC50 [GDSC-CCLE] & Pred_Error (color)", dir)

add = list(scale_color_gradient(low="beige", high="firebrick3"))
midpoint = Pred_GCNPath_N_CCLE %>% with(abs(Prediction-LN_IC50)) %>% quantile(0.01)
# add = list(scale_color_gradient2(low="beige", mid="lightyellow", high="firebrick3"))

Pred_GCNPath_N_CCLE %>% 
  mutate(Pred_Error=abs(Prediction-LN_IC50)) %>% 
  plot_def(LN_IC50, LN_IC50_CCLE, color=Pred_Error, 
           main=main, xlab=xlab, ylab=ylab, size=2.5, alpha=0.8, 
           legend="|Error|", color_line="red", add=add,
           axis_tl=30, axis_tx=24, legend_tl=20, legend_tx=20, margin_lg=0.4,
           width=20, height=16.5, dpi=1200, xy_line=T, raster=T, save=T, save_svg=T)

xlab = bquote("|"~GDSC~ln(IC[50]) - " " * CCLE~ln(IC[50])~"|")
ylab = bquote("|"~GDSC~ln(IC[50]) - " " * Prediction~"|")
main = sprintf("%s/LN_IC50 [GDSC-CCLE] & Pred_Error (no color)", dir)

Pred_GCNPath_N_CCLE %>% 
  mutate(Pred_Error=abs(Prediction-LN_IC50), 
         LN_IC50_Var=abs(LN_IC50-LN_IC50_CCLE)) %>% 
  plot_def(LN_IC50_Var, Pred_Error, main=main, xlab=xlab, ylab=ylab, 
           size=2.5, alpha=0.5, axis_tl=22.5, axis_tx=22.5, 
           width=16, height=16, dpi=1500, raster=T, save=T, save_svg=T)

Pred_GCNPath_N_CCLE %>% 
  mutate(Pred_Error=abs(Prediction-LN_IC50), 
         LN_IC50_Var=abs(LN_IC50-LN_IC50_CCLE)) %>% 
  with(cor(Pred_Error, LN_IC50_Var)) %>% round(3)   # 0.144 [n=5567]



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
  
  process_perf = function(Perf_List, score="RMSE") {
    Perf_List_ = Perf_List %>% 
      rename(Train_Fold=Fold, Num_Test=N_Test) %>% 
      mutate(Dataset=recode(Dataset, "GDSC"="GDSC1+2")) %>% 
      subset(select=c(Model, Dataset, Test_Type, Train_Fold, Num_Test, object(score)))   # 2310
    
    test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
    Perf_List_N = Perf_List_ %>% subset(Test_Type==test_type[1]) %>% subset(select=-Test_Type)
    Perf_List_C = Perf_List_ %>% subset(Test_Type==test_type[2]) %>% subset(select=-Test_Type)
    Perf_List_D = Perf_List_ %>% subset(Test_Type==test_type[3]) %>% subset(select=-Test_Type)
    Perf_List_S = Perf_List_ %>% subset(Test_Type==test_type[4]) %>% subset(select=-Test_Type)
    
    Perf_List_ = list(Perf_List_N, Perf_List_C, Perf_List_D, Perf_List_S)
    return(Perf_List_)
  }
  
  process_utest = function(Perf_UTest, score="RMSE") {
    
    score_ = sprintf("_%s", score)
    col_rm = c("RMSE", "PCC", "SCC") %>% setdiff(score)
    
    Perf_UTest_ = Perf_UTest %>% 
      subset(Dataset=="GDSC") %>% subset(select=-Dataset) %>% 
      select(!contains(col_rm)) %>% rename_with(~gsub(score_, "", .x)) %>% 
      rename(Minus_Log10_FDR=MLog10_FDR) %>% 
      relocate(Minus_Log10_FDR, .after=FDR) %>% as.data.frame
    
    test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
    Perf_UTest_N = Perf_UTest_ %>% subset(Test_Type==test_type[1]) %>% subset(select=-Test_Type)
    Perf_UTest_C = Perf_UTest_ %>% subset(Test_Type==test_type[2]) %>% subset(select=-Test_Type)
    Perf_UTest_D = Perf_UTest_ %>% subset(Test_Type==test_type[3]) %>% subset(select=-Test_Type)
    Perf_UTest_S = Perf_UTest_ %>% subset(Test_Type==test_type[4]) %>% subset(select=-Test_Type)
    
    Perf_UTest_List = list(Perf_UTest_N, Perf_UTest_C, Perf_UTest_D, Perf_UTest_S)
    return(Perf_UTest_List)
  }
  
  ### [Source Data] Supplementary Fig. 8
  Perf_RMSE_ = Perf_List %>% process_perf(score="RMSE")
  Perf_RMSE_ %>% save_for_nc(num=8, suppl=T)
  
  ### [Source Data] Supplementary Fig. 9
  Perf_PCC_ = Perf_List %>% process_perf(score="PCC")
  Perf_PCC_ %>% save_for_nc(num=9, suppl=T)
  
  ### [Source Data] Supplementary Fig. 10
  Perf_SCC_ = Perf_List %>% process_perf(score="SCC")
  Perf_SCC_ %>% save_for_nc(num=10, suppl=T)
  
  ### [Source Data] Supplementary Fig. 11
  Perf_UTest_RMSE_ = Perf_UTest %>% process_utest(score="RMSE")
  Perf_UTest_RMSE_ %>% save_for_nc(num=11, suppl=T)
  
  ### [Source Data] Supplementary Fig. 12
  Perf_UTest_PCC_ = Perf_UTest %>% process_utest(score="PCC")
  Perf_UTest_PCC_ %>% save_for_nc(num=12, suppl=T)
  
  ### [Source Data] Supplementary Fig. 13
  Perf_UTest_SCC_ = Perf_UTest %>% process_utest(score="SCC")
  Perf_UTest_SCC_ %>% save_for_nc(num=13, suppl=T)
  
  ### [Source Data] Supplementary Fig. 17
  Perf_Cell_Corr_ = Perf_Cell_Corr %>% 
    lapply(function(df) df[, 1:3] %>% 
             rename(Model_X=Model1, Model_Y=Model2, PCC_of_RMSE=RMSE_Corr))
  Perf_Cell_Corr_ %>% save_for_nc(num=17, suppl=T)
  
  ### [Source Data] Supplementary Fig. 18
  Perf_Drug_Corr_ = Perf_Drug_Corr %>% 
    lapply(function(df) df[, 1:3] %>% 
             rename(Model_X=Model1, Model_Y=Model2, PCC_of_RMSE=RMSE_Corr))
  Perf_Drug_Corr_ %>% save_for_nc(num=18, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 19
  IC50_Stat_Cell_ = IC50_Stat_Cell %>% 
    rename(Cancer_Type=Cell_TCGA, Num_Test=Num, 
           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD)
  IC50_Stat_Drug_ = IC50_Stat_Drug %>% 
    rename(Target_Pathway=Drug_Pathway, Num_Test=Num, 
           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD)
  
  IC50_Stat_Cell_Group_ = IC50_Stat_Cell_Group %>% 
    rename(Cancer_Type=Cell_TCGA, Num_Test=Num, 
           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD)
  IC50_Stat_Drug_Group_ = IC50_Stat_Drug_Group %>% 
    rename(Target_Pathway=Drug_Pathway, Num_Test=Num, 
           LN_IC50_Mean=IC50_Mean, LN_IC50_SD=IC50_SD)
  
  IC50_Stat_ = list(IC50_Stat_Drug_, IC50_Stat_Drug_Group_, 
                    IC50_Stat_Cell_, IC50_Stat_Cell_Group_)
  
  num_fig = c("a-c", "d", "e-g", "h")
  IC50_Stat_ %>% save_for_nc(num=19, num_fig=num_fig, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 20
  Pred_GCNPath_N_CCLE_ = Pred_GCNPath_N_CCLE %>% 
    subset(select=-c(Model, Dataset, Test_Type, Cell_COSMIC, Cell_TCGA, Drug_Pathway)) %>% 
    rename(Train_Fold=Fold, LN_IC50_GDSC=LN_IC50) %>% 
    mutate(Diff_GDSC_Pred=abs(Prediction-LN_IC50_GDSC), 
           Diff_GDSC_CCLE=abs(LN_IC50_CCLE-LN_IC50_GDSC)) %>% 
    relocate(Train_Fold, .after=Drug) %>% 
    relocate(Prediction, .after=LN_IC50_CCLE) %>% as.data.frame
  
  Pred_GCNPath_N_CCLE_ %>% save_for_nc(num=20, suppl=T)
  
  
  ### [Source Data] Fig. 3
  subset_db = function(Perf) Perf %>% subset(Dataset=="GDSC1+2") %>% subset(select=-Dataset)
  Perf_RMSE_GDSC_ = Perf_RMSE_ %>% lapply(subset_db)
  Perf_RMSE_GDSC_ %>% sapply(nrow)   # 140, 140, 140, 350
  Perf_RMSE_GDSC_ %>% save_for_nc(num=3, suppl=F)
}

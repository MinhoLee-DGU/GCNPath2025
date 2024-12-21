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
    summarize(N_Test = n(), RMSE = RMSE(LN_IC50, Prediction)) %>% as.data.frame
  
  lvl = c("GDSC2 x Normal", "GDSC2 x Strict_Blind")
  Perf = Perf %>% mutate(Test = sprintf("%s x %s", Dataset, Test_Type))
  Perf$Test = Perf$Test %>% factor(levels=lvl)
  return(Perf)
}

calc_perf_fold = function(Pred) {
  Perf = Pred %>% group_by(Model, Dataset, Test_Type, Fold) %>%
    summarize(N_Test = n(), RMSE = RMSE(LN_IC50, Prediction)) %>% as.data.frame
  
  lvl = c("GDSC2 x Normal", "GDSC2 x Strict_Blind")
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

boxplot_perf = function(Perf_Fold, x="Model", y="RMSE", fill="Model", 
                        main=NULL, mark_def=T, color_def="royalblue3",
                        add=NULL, hide_legend=T, width=32, height=12, 
                        axis_tl=30, axis_tx=18, legend_tl=20, legend_tx=18, save=T) {
  
  # legend = "Parameters"
  # p1 = Perf_Fold %>% subset(Test=="GDSC2 x Normal") %>% 
  #   boxplot_def(Model, RMSE, fill=Model, main=main, legend=legend, width=width/5.2, height=height, save=F)
  # p2 = Perf_Fold %>% subset(Test=="GDSC2 x Strict_Blind") %>% 
  #   boxplot_def(Model, RMSE, fill=Model, main=main, legend=legend, width=width/5.2, height=height, save=F)
  
  add.params = list(alpha=0.2, color="black")
  margin = margin(8, 8, 8, 8, unit="pt")
  font1 = font("ylab", color="grey30", size=axis_tl, margin=margin)
  font2 = font("axis.text", color="grey30", size=axis_tx, margin=margin)
  
  font = font1 + font2
  margin_pl = margin(0.2, 0.4, 0.2, 0.4, "cm")
  margin_pl = theme(plot.margin=margin_pl)
  # fill = ifelse(mark_def, "Model", "white")
  
  # color="Model"
  p1 = Perf_Fold %>% subset(Test=="GDSC2 x Normal") %>% 
    ggboxplot(x, y, fill=fill, add="point", outlier.shape=NA,
              xlab=F, add.params=add.params) + font + margin_pl + rotate_x_text(45)
  p2 = Perf_Fold %>% subset(Test=="GDSC2 x Strict_Blind") %>% 
    ggboxplot(x, y, fill=fill, add="point", outlier.shape=NA,
              xlab=F, add.params=add.params) + font + margin_pl + rotate_x_text(45)
  
  if (mark_def) {
    param_list = Perf_Fold$Model %>% unique %>% as.character
    param_def = grep("\\[DE\\]", param_list, value=T)
    fill = c(color_def, rep("#66b3ed", length(param_list)-1))
    names(fill) = c(param_def, param_list[param_list!=param_def])
    p1 = p1 + scale_fill_manual(values=fill)
    p2 = p2 + scale_fill_manual(values=fill)
  } 
  
  if (!is.null(add)) {
    for (i in 1:length(add)) {
      p1 = p1 + add[[i]]
      p2 = p2 + add[[i]]
    }
  }
  
  if (hide_legend) {
    no_lg = function(pl) pl %>% ggpar(legend="none")
    pl = ggarrange(no_lg(p1), no_lg(p2), ncol=2)
    # common.legend=T, legend="right"
  } else {
    margin_lg = 0.4
    margin_lgl = margin(b=margin_lg, unit="cm")
    margin_lgx = margin(b=margin_lg, l=margin_lg, unit="cm")
    theme_ = theme(legend.title = element_text(size=legend_tl, margin=margin_lgl),
                   legend.text = element_text(size=legend_tx, margin=margin_lgx))
    
    add_lg = function(pl) pl + theme_
    pl = ggarrange(add_lg(p1), add_lg(p2), ncol=2, common.legend=T, legend="right")
  }
  
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

dir_res = "../results"
dir_fig = mkdir("Ablation Test")
dir_test = c("IC50_GDSC2/Normal", "IC50_GDSC2/Strict_Blind")


# Ablation Test [GNN vs Linear]
gnn_lin = c("RGCN", "RGCN_MF256", "RGCN_MF512", 
            "RGCN_SMILESVec", "Linear_GAT", "Linear_MF256")

re_label = list("RGCN"="RGCN_GAT")
dir_list = to_dir_list(dir_res, dir_test, gnn_lin)
Pred_GNN_Lin = summary_pred(dir_list, lvl_param=gnn_lin, re_label=re_label)

group = c("GCN_C_D", "GCN_C", "GCN_C", "GCN_C", "GCN_D", "NULL")
df_col = data.frame(Model=gnn_lin, Group=group)
df_col$Model[df_col$Model=="RGCN"] = "RGCN_GAT"

lvl_group = c("GCN_C_D", "GCN_C", "GCN_D", "NULL")
idx = match(Pred_GNN_Lin$Perf_Test_Fold$Model, df_col$Model)
Pred_GNN_Lin$Perf_Test_Fold$Group = df_col$Group[idx] %>% factor(levels=lvl_group)

# color = RColorBrewer::brewer.pal(8, "Blues")[c(8, 6, 4, 2)]
color1 = c("royalblue3", "#66b3ed")
color2 = RColorBrewer::brewer.pal(8, "Reds")[c(7, 4)]
color = c(color1, color2)
names(color) = lvl_group

labels = c(bquote(GCN[C] ~ "×" ~ GCN[D]), 
           bquote(GCN[C] ~ "×" ~ FCN[D]), 
           bquote(FCN[C] ~ "×" ~ GCN[D]), 
           bquote(FCN[C] ~ "×" ~ FCN[D]))

names(labels) = lvl_group
add = list(scale_fill_manual(values=color, labels=labels))

main = sprintf("%s/Ablation [GNN vs Linear]", dir_fig)
Pred_GNN_Lin$Perf_Test_Fold %>% boxplot_perf(main=main, width=32, save=T, fill="Group", add=add, hide_legend=F)
Perf_GNN_Lin = perf_avg_sd(Pred_GNN_Lin$Perf_Test_Fold, space_cols=T)



### Save all files...

suppressMessages(library(openxlsx))
file = sprintf("%s/Ablation Summary.xlsx", dir_fig)
write.xlsx(Perf_GNN_Lin, file=file, rowNames=T, colNames=T)


supplementary = T
if (supplementary) {
  
  process_perf = function(Perf) {
    Perf_Test = Perf$Perf_Test_Fold %>% 
      subset(select=-c(Test, R2, Corr)) %>% 
      rename(Parameter=Model, Train_Fold=Fold, Num_Test=N_Test)
    
    Perf_Test_G2N = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Normal")
    Perf_Test_G2S = Perf_Test %>% subset(Dataset=="GDSC2" & Test_Type=="Strict_Blind") 
    
    Perf_List = list(Perf_Test_G2N, Perf_Test_G2S) %>% 
      lapply(function(df) df %>% subset(select=-c(Dataset, Test_Type)))
    
    params = unique(as.character(Perf_Test_G2N$Parameter))
    sprintf("# Num [Test] : %s, %s", nrow(Perf_Test_G2N), nrow(Perf_Test_G2S)) %>% print
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
  Perf_GNN_Lin_ = Pred_GNN_Lin_ %>% process_perf
  Perf_GNN_Lin_[[1]]$Group = NULL
  Perf_GNN_Lin_[[2]]$Group = NULL
  Perf_GNN_Lin_ %>% save_for_nc(num=2, suppl=F)
}

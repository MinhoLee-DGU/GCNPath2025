#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

source("../functions.R")
loadings()


read_pred = function(dir, pattern, model=NULL, idx_weight=1, sep=",") {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_[^CCLE]", "Pred_CCLE"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      try({
        Pred_TP = fread(df_name, header=T, sep=sep)
        weight = strsplit(dir, "/")[[1]] %>% tail(idx_weight)
        Pred_TP$Model = model
        Pred_TP$Weight = weight[1]
        Pred = Pred %>% rbind(Pred_TP)
      })
    }
    return(Pred)
  }
}

read_pred_all = function(dir_list, pattern, model=NULL, weight=NULL, idx_weight=1, sep=",") {
  
  Pred = data.frame()
  # weight : [character] Rename model having the following weight
  # weight = c("SA", "SA_pre")
  # names(weight) = c("TGSA", "TGSA_Pre")
  # Ex. SA > TGSA, SA_pre > TGSA_Pre
  
  for (dir in dir_list) {
    Pred_TP = read_pred(dir, pattern, model=model, idx_weight=idx_weight, sep=sep)
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
      
      Pred_TP$Dataset = dataset
      Pred_TP = Pred_TP %>% relocate(Model, Weight, Dataset, .before=everything())
      Pred = Pred %>% rbind(Pred_TP)
    }
  }
  
  if (!is.null(weight) & !is.null(names(weight))) {
    for (i in 1:length(weight)) {
      Pred$Model[Pred$Weight==weight[i]] = names(weight)[i]
    }
  }
  
  Perf = Pred %>% calc_perf
  Perf_Cell = Pred %>% calc_perf(option=1)
  Perf_Drug = Pred %>% calc_perf(option=2)
  Pred = list(Pred=Pred, Perf=Perf, Perf_Cell=Perf_Cell, Perf_Drug=Perf_Drug)
  return(Pred)
}

calc_perf = function(Pred, option=0, ignore_weight=F) {
  
  if (option==0) {
    col_by = c("Model", "Dataset", "Weight")
  } else if (option==1) {
    col_by = c("Model", "Dataset", "Weight", "Cell")
  } else  {
    col_by = c("Model", "Dataset", "Weight", "Drug")
  }
  
  if (ignore_weight) col_by = col_by[col_by!="Weight"]
  
  Perf = Pred %>% 
    group_by(across(all_of(col_by))) %>% 
    summarize(N_Test = n(),
              RMSE = RMSE(LN_IC50, Prediction),
              MAE = MAE(LN_IC50, Prediction),
              R2 = R2(LN_IC50, Prediction),
              PCC = cor(LN_IC50, Prediction), 
              SCC = cor(LN_IC50, Prediction, method="spearman"))
  
  Perf = Perf %>% as.data.frame
  return(Perf)
}

to_dir_list = function(dir, ...) {
  dir_list = expand.grid(...) %>% arrange(Var1)
  dir_list = dir_list %>% apply(1, function(...) paste0(..., collapse="/"))
  dir_list = sprintf("%s/%s", dir, dir_list)
  return(dir_list)
}

# barplot_def = function(df, x, y, fill=NULL, main=NULL, xlab=NULL, ylab=NULL, 
#                        text=NULL, add=NULL, size=1, alpha=1, color="black", 
#                        plot_tl=22.5, axis_tl=18, axis_tx=15, legend_tl=16.5, legend_tx=15, 
#                        margin=0.25, margin_pl=0.25, plot_face="plain", axis_face="plain", legend_face="plain",
#                        width=15, height=15, show_title=F, reorder=T, force_bold=F, save=T) {
# 
#   title = NULL
#   e_text = element_text
#   # if (!is.null(text)) text = sprintf("\n%s", text)
#   
#   if (show_title) { 
#     title = sprintf("%s\n", main)
#     title = title %>% strsplit("/") %>% unlist
#     title = title[length(title)]
#   }
#   
#   mapping = aes(x={{x}}, y={{y}})
#   fill_string = deparse(substitute(fill))
#   reorder_string = deparse(substitute(reorder))
#   
#   margin_lg = margin(b=0.2, unit="cm")
#   margin_x = margin(t=margin, b=margin, unit="cm")
#   margin_y = margin(l=margin, r=margin, unit="cm")
#   
#   xlab = ifelse_def(!is.null(xlab), xlab, deparse(substitute(x)))
#   ylab = ifelse_def(!is.null(ylab), ylab, deparse(substitute(y)))
#   
#   if (force_bold) {
#     plot_face = "bold"; axis_face = "bold"; 
#     legend_face = "bold" ; text_face = "bold"; 
#   }
#   
#   if (fill_string=="NULL" & !reorder) {
#     mapping = aes(x={{x}}, y={{y}}, fill={{x}})
#   } else if (fill_string!="NULL" & !reorder) {
#     mapping = aes(x={{x}}, y={{y}}, fill={{fill}})
#   } else if (fill_string=="NULL" & reorder) {
#     mapping = aes(x=reorder({{x}}, {{y}}), y={{y}}, fill=reorder({{x}}, {{y}}))
#   } else {
#     mapping = aes(x=reorder({{x}}, {{y}}), y={{y}}, fill=reorder({{fill}}, {{y}}))
#   }
#   
#   pl = ggplot(df, mapping) + 
#     theme_classic() + labs(x=xlab, y=ylab, title=title) + 
#     geom_bar(stat="identity", width=0.9) +
#     theme(plot.margin = unit(rep(margin_pl, 4), units="cm"),
#           plot.title = e_text(size=plot_tl, face=plot_face, hjust=0.5),
#           axis.title = e_text(size=axis_tl, face=axis_face, hjust=0.5), 
#           legend.title = e_text(size=legend_tl, face=legend_face, margin=margin_lg),
#           legend.text = e_text(size=legend_tx, face=legend_face, margin=margin_lg),
#           axis.text.x = e_text(size=axis_tx, face=axis_face, margin=margin_x), 
#           axis.text.y = e_text(size=axis_tx, face=axis_face, margin=margin_y))
#   
#   if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
#   
#   if (save) {
#     main = gsub("\n", " ", main)
#     file = sprintf("%s.png", main)
#     ggsave(file, pl, width=width, height=height, units="cm")
#   } else print(pl)
# }



##### 2. Read Prediction

pattern = "pred_test.csv"
dataset = sprintf("IC50_GDSC%s", c("", 1, 2))

# PaccMann_V2_MSE/PCC
model = "PaccMann"
dir = "../PaccMann/Results"
weight = c("best_mse_paccmann_v2", "best_pearson_paccmann_v2")

names(weight) = c("PaccMann_MSE", "PaccMann_PCC")
dir_list = to_dir_list(dir=dir, dataset, weight)
Pred_PaccMann = read_pred_all(dir_list, pattern, model=model, weight=weight)
# Pred_PaccMann$Pred = Pred_PaccMann$Pred[, -c("V1")]

# GraphDRP
model = "GraphDRP"
dir = "../GraphDRP/Results"
weight = "GINConv"

names(weight) = "GraphDRP"
dir_list = to_dir_list(dir=dir, dataset, weight)
Pred_GraphDRP = read_pred_all(dir_list, pattern, model=model, weight=weight)

# TGSA [TGDRP]
model = "TGDRP"
dir = "../TGSA/Results"
weight = c("TGDRP", "TGDRP_pre")

names(weight) = c("TGDRP", "TGDRP_Pre")
dir_list = to_dir_list(dir=dir, dataset, weight)
Pred_TGDRP = read_pred_all(dir_list, pattern, model=model, weight=weight)

# TGSA [TGSA]
model = "TGSA"
dir = "../TGSA/Results"
weight = c("SA", "SA_pre")

names(weight) = c("TGSA", "TGSA_Pre")
dir_list = to_dir_list(dir=dir, dataset, weight)
Pred_TGSA = read_pred_all(dir_list, pattern, model=model, weight=weight)

# DRPreter
model = "DRPreter"
dir = "../DRPreter/Results"
weight = c(2, 16, 33, 61, 79, 100, 220, 653, 1004, 4001)
weight = sprintf("Seed%s/DRPreter", weight)

dir_list = to_dir_list(dir=dir, dataset, weight)
Pred_DRPreter = read_pred_all(dir_list, pattern, model=model, idx_weight=2)



### Barplot + SD
dataset = c("GDSC", "GDSC1", "GDSC2")
Perf_Pretrain = Reduce(rbind, list(Pred_PaccMann$Perf, Pred_GraphDRP$Perf, 
                                   Pred_TGDRP$Perf, Pred_TGSA$Perf, Pred_DRPreter$Perf))

plot_perf_pre = function(Perf_Pretrain, y="RMSE", axis_tl=24, axis_tx=18,
                         legend_tl=20, legend_tx=20, width=20, height=15, save=T) {
  
  # width=30, height=18
  margin1 = margin(r=0.25, l=0.25, unit="cm")
  margin2 = margin(0.25, 0.25, 0.25, 0.25, unit="cm")
  margin3 = margin(r=1, unit="cm")
  margin4 = margin(r=0.25, l=0.25, unit="cm")
  
  font1 = font(object="ylab", size=axis_tl, margin=margin1)
  font2 = font(object="axis.text", size=axis_tx, margin=margin2, color="grey30")
  font3 = font(object="legend.title", size=legend_tl, margin=margin3)
  font4 = font(object="legend.text", size=legend_tx, margin=margin4)
  font = font1 + font2 + font3 + font4
  
  font_text = element_text(size=7.2)
  # ymax = max(Perf_Pretrain[[y]]) + 0.2
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
  main = sprintf("Performance [Barplot, %s]", y)
  
  pl = Perf_Pretrain %>% 
    ggbarplot(x="Model", y=y, color="black", fill="Dataset",
              xlab=F, add="mean_se", position=position_dodge(0.8), legend="right") + 
    theme(text=font_text) + font + rotate_x_text(36) + scale_fill_manual(values=color)
  
  pl = pl %>% ggpar(legend="bottom")
  # pl = Perf_Pretrain %>% 
  #   ggbarplot(x="Model", y=y, color="black", fill="Dataset",
  #             xlab=F, add="mean_se", label=T, lab.vjust=-1.2, lab.nb.digits=2, 
  #             fontface="bold", ylim=c(0, ymax), position=position_dodge(0.8), legend="right") + 
  #   theme(text=font_text) + font + rotate_x_text(40) + scale_fill_manual(values=color)
  
  if (save) {
    pl %>% save_fig_ggpubr(main=main, width=width, height=height, svg=T)
  } else print(pl)
}

perf_avg_sd = function(Perf_Avg, stat="RMSE", space_cols=F) {
  
  avg_plus_sd = function(x_mean, x_sd) sprintf("%.3f±%.3f", x_mean, x_sd)
  Perf_Fold_Avg = Perf_Fold %>% acast(Model~Test, value.var=stat, fun.aggregate=mean)
  Perf_Fold_SD = Perf_Fold %>% acast(Model~Test, value.var=stat, fun.aggregate=sd)
  
  Perf_Fold_Avg = Perf_Fold_Avg %>% as.data.frame
  Perf_Fold_SD = Perf_Fold_SD %>% as.data.frame
  Perf_Fold_Sum = mapply(avg_plus_sd, Perf_Fold_Avg, Perf_Fold_SD)
  if (space_cols) colnames(Perf_Fold_Sum) = gsub(" x ", "\n", colnames(Perf_Fold_Sum))
  
  Perf_Fold_Sum = Perf_Fold_Sum %>% as.data.frame
  rownames(Perf_Fold_Sum) = rownames(Perf_Fold_Avg)
  return(Perf_Fold_Sum)
}

Perf_Pretrain %>% plot_perf_pre(y="RMSE")
Perf_Pretrain %>% plot_perf_pre(y="PCC")
Perf_Pretrain %>% plot_perf_pre(y="SCC")

col = c("RMSE", "MAE", "R2", "PCC", "SCC")
avg_plus_sd_ = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
avg_plus_sd = function(x) ifelse(length(x)!=1, avg_plus_sd_(x), sprintf("%.3f", x))

Perf_Pretrain_Avg = Perf_Pretrain %>% group_by(Model, Dataset) %>% 
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame

file = "Performances Summary [Pretrained Model].xlsx"
write.xlsx(Perf_Pretrain_Avg, file=file, rowNames=F)




### Scatter Plot

plot_pred = function(Pred, model=NULL, dir=NULL, width=15, height=15) {
  
  xlab = bquote(Actual~ln(IC[50]))
  ylab = bquote(Predicted~ln(IC[50]))
  main = sprintf("%s/Prediction [%s, Pretrained]", dir, model)
  
  rmse = Pred %>% with(RMSE(LN_IC50, Prediction)) %>% round(3)
  corr = Pred %>% with(cor(LN_IC50, Prediction)) %>% round(3)
  sprintf("# N=%s, RMSE=%.3f, PCC=%.3f", nrow(Pred), rmse, corr) %>% print
  
  Pred %>% plot_def(LN_IC50, Prediction, main=main, xlab=xlab, ylab=ylab,
                    size=1.5, alpha=0.25, force_bold=F, axis_tl=30, axis_tx=24, dpi=1200,
                    width=width, height=height, xy_line=T, raster=T, save=T, save_svg=T)
}

dir = mkdir("Prediction [GSDC]")

# PaccMann_MSE
# N=150204, RMSE=2.074, PCC=0.642
Pred_PaccMann$Pred %>%
  subset(Dataset=="GDSC" & Weight=="best_mse_paccmann_v2") %>% 
  plot_pred(model="PaccMann_MSE", dir=dir)

# PaccMann_PCC
# N=150204, RMSE=2.089, PCC=0.641
Pred_PaccMann$Pred %>%
  subset(Dataset=="GDSC" & Weight=="best_pearson_paccmann_v2") %>% 
  plot_pred(model="PaccMann_PCC", dir=dir)

# GraphDRP
# N=371803, RMSE=2.074, PCC=0.679
Pred_GraphDRP$Pred %>%
  subset(Dataset=="GDSC") %>% 
  plot_pred(model="GraphDRP", dir=dir)

# TGDRP
# N=269519, RMSE=2.447, PCC=0.510
Pred_TGDRP$Pred %>%
  subset(Dataset=="GDSC" & Model=="TGDRP") %>%
  plot_pred(model="TGDRP", dir=dir)

# TGDRP_Pre
# N=269519, RMSE=2.404, PCC=0.501
Pred_TGDRP$Pred %>%
  subset(Dataset=="GDSC" & Model=="TGDRP_Pre") %>%
  plot_pred(model="TGDRP_Pre", dir=dir)

# TGSA
# N=269519, RMSE=2.386, PCC=0.497
Pred_TGSA$Pred %>%
  subset(Dataset=="GDSC" & Model=="TGSA") %>%
  plot_pred(model="TGSA", dir=dir)

# TGSA_Pre
# N=269519, RMSE=2.408, PCC=0.526
Pred_TGSA$Pred %>%
  subset(Dataset=="GDSC" & Model=="TGSA_Pre") %>%
  plot_pred(model="TGSA_Pre", dir=dir)

# DRPreter
# N=267965, RMSE=2.176, PCC=0.635
seed_best = Pred_DRPreter$Perf %>% 
  subset(Dataset=="GDSC") %>% 
  subset(RMSE==min(RMSE)) %>% pull(Weight)   # Seed 4001

Pred_DRPreter$Pred %>%
  subset(Dataset=="GDSC" & Weight==seed_best) %>%
  plot_pred(model="DRPreter", dir=dir)


# Prediction
save_pred = T
if (save_pred) {
  rbind_perf = function(Perf1, Perf2, col=NULL) {
    if (is.null(col)) col = intersect(colnames(Perf1), colnames(Perf2))
    if (is.data.table(Perf1) & is.data.table(Perf2)) {
      Perf = rbind(Perf1[, col, with=F], Perf2[, col, with=F])
    } else Perf = rbind(Perf1[, col], Perf2[, col])
    return(Perf)
  }
  
  Pred_List = list()
  model_names = c("PaccMann", "GraphDRP", "TGDRP", "TGSA", "DRPreter")
  
  pred_names = sprintf("Pred_%s", model_names)
  for (i in 1:length(model_names)) Pred_List[[i]] = get(pred_names[i])$Pred
  Pred_List = Reduce(rbind_perf, Pred_List)   # 9955304 x 7
  
  file = "Prediction [GDSC, Pre-trained Model].csv"
  fwrite(Pred_List, file=file, row.names=F)
  
  rm(Pred_List) ; gc()
}



### Compare Test Prediction [Pre-train vs Train, GDSC1+2]

compare_pretrain_multi = function(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                  model="DRPreter", pretrain="DRPreter", 
                                  dataset="GDSC", dir=NULL, save=T) {
  
  corr_cell = c()
  corr_drug = c()
  test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
  
  for (i in 1:length(test_type)) {
    corr_cell_ = compare_pretrain(Perf_Cell, Perf_Cell_Pre, model=model, pretrain=pretrain, 
                                  test_type=test_type[i], by="Cell", dir=dir, save=save)
    corr_drug_ = compare_pretrain(Perf_Drug, Perf_Drug_Pre, model=model, pretrain=pretrain, 
                                  test_type=test_type[i], by="Drug", dir=dir, save=save)
    
    corr_cell = corr_cell %>% c(corr_cell_)
    corr_drug = corr_drug %>% c(corr_drug_)
  }
  
  Corr_Test = data.frame(Corr_Cell_RMSE=corr_cell, Corr_Drug_RMSE=corr_drug)
  rownames(Corr_Test) = test_type
  return(Corr_Test)
}

compare_pretrain = function(Perf, Perf_Pre, model="DRPreter", pretrain="DRPreter", 
                            dataset="GDSC", test_type="Normal", by="Cell", dir=NULL, save=T) {
  
  Perf = Perf %>% subset(Model==model & Dataset==dataset & Test_Type==test_type)
  Perf_Pre = Perf_Pre %>% subset(Test_Type==pretrain & Dataset==dataset)
  col = intersect(colnames(Perf), colnames(Perf_Pre))
  Perf = rbind(Perf[, col], Perf_Pre[, col])
  
  by_group = ifelse(by=="Cell", "Cell_TCGA", "Drug_Pathway")
  Perf[[by_group]][is.na(Perf[[by_group]])] = "Unknown"
  N_By_Group = Perf %>% group_by(across(by_group)) %>% summarise(Num=length(unique(Test_Type)))
  groups = N_By_Group %>% subset(Num==2) %>% pull(by_group)
  Perf = Perf[Perf[[by_group]] %in% groups, ]
  
  lvl = Perf %>% subset(Test_Type==pretrain) %>% 
    group_by(across(by_group)) %>% 
    summarise(RMSE_Median=median(RMSE)) %>% 
    arrange(RMSE_Median) %>% pull(by_group)
  
  Perf[[by_group]] = Perf[[by_group]] %>% factor(levels=lvl)
  # Perf %>% compare_boxplot_pre(model, pretrain=pretrain, test_type=test_type, by=by, dir=dir, save=save)
  corr = Perf %>% compare_scatter_pre(model, pretrain=pretrain, test_type=test_type, by=by, dir=dir, save=save)
  return(corr)
}

compare_scatter_pre = function(Perf, model="DRPreter", pretrain="DRPreter",
                               test_type="Normal", by="Cell", dir=NULL,
                               axis_tl=20, axis_tx=18, width=16, height=16, dpi=1500, save=T) {

  xlab = sprintf("%s RMSE [Pretrained Model]", by)
  ylab = sprintf("%s RMSE [%s Test]", by, gsub("_", "-", test_type))
  add = list(stat_smooth(method="lm"))

  main = sprintf("%s/Scatter Comparison of %s RMSE [%s, Pretrain & %s]", dir, by, pretrain, test_type)
  Perf = Perf %>% acast(as.formula(paste(by, "Test_Type", sep="~")), value.var="RMSE")
  Perf = Perf[complete.cases(Perf), ] %>% as.data.frame

  corr = cor(as.numeric(Perf[, 1]), as.numeric(Perf[, 2]))
  sprintf("[Pretrain & %s Test] Corr : %s", test_type, round(corr, 3)) %>% print

  # Perf %>% plot_def(object(pretrain), object(test_type),
  #                   main=main, xlab=xlab, ylab=ylab,
  #                   alpha=0.5, axis_tl=axis_tl, axis_tx=axis_tx, add=add,
  #                   force_bold=F, width=width, height=height, dpi=dpi, save=save)

  return(corr)
}

# compare_boxplot_pre = function(Perf, model="DRPreter", pretrain="DRPreter", 
#                                test_type="Normal", by="Cell", dir=NULL, 
#                                axis_tl=20, axis_tx=20, legend_tl=18, legend_tx=18,
#                                width=54, height=20, dpi=400, save=T) {
#   
#   legend = "Test Type"
#   Perf$Test_Type[Perf$Test_Type==pretrain] = "Pretrain"
#   Perf$Test_Type = Perf$Test_Type %>% factor(levels=c("Pretrain", test_type))
#   main = sprintf("%s/Boxplot Comparison of %s RMSE [%s, Pretrain & %s]", dir, by, pretrain, test_type)
#   
#   if (by=="Cell") {
#     Perf %>% boxplot_def(Cell_TCGA, RMSE, fill=Test_Type, main=main, 
#                          xlab=NULL, ylab="Cell RMSE", axis_tl=axis_tl, axis_tx=axis_tx, 
#                          legend_tl=legend_tl, legend_tx=legend_tx, legend=legend, 
#                          width=width, height=height, reorder=F, force_bold=F, dpi=dpi, save=save)
#   } else {
#     Perf %>% boxplot_def(Drug_Pathway, RMSE, fill=Test_Type, main=main, 
#                          xlab=NULL, ylab="Drug RMSE", axis_tl=axis_tl, axis_tx=axis_tx, 
#                          legend_tl=legend_tl, legend_tx=legend_tx, legend=legend, 
#                          width=width, height=height, reorder=F, force_bold=F, dpi=dpi, save=save)
#   }
# }


dir = "../_performance/Performance [Cell, Drug]"

file= sprintf("%s/Performance_Cell.csv", dir)
Perf_Cell = read.csv(file)

file= sprintf("%s/Performance_Drug.csv", dir)
Perf_Drug = read.csv(file)

file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)

file = "Anno_Drugs.csv"
Anno_Drugs = read.csv(file)

Perf_Cell_DRPreter = Pred_DRPreter$Pred %>% 
  subset(Weight==seed_best) %>% 
  calc_perf(option=1) %>% mutate(Weight="DRPreter")

Perf_Drug_DRPreter = Pred_DRPreter$Pred %>% 
  subset(Weight==seed_best) %>% 
  calc_perf(option=2) %>% mutate(Weight="DRPreter")

test_type = c("Normal", "Cell_Blind", "Drug_Blind", "Strict_Blind")
Perf_Cell_Pre = Reduce(rbind, list(Pred_PaccMann$Perf_Cell, Pred_GraphDRP$Perf_Cell, 
                                   Pred_TGDRP$Perf_Cell, Pred_TGSA$Perf_Cell, Perf_Cell_DRPreter))
Perf_Drug_Pre = Reduce(rbind, list(Pred_PaccMann$Perf_Drug, Pred_GraphDRP$Perf_Drug, 
                                   Pred_TGDRP$Perf_Drug, Pred_TGSA$Perf_Drug, Perf_Drug_DRPreter))

colnames(Perf_Cell_Pre)[colnames(Perf_Cell_Pre)=="Weight"] = "Test_Type"
colnames(Perf_Drug_Pre)[colnames(Perf_Drug_Pre)=="Weight"] = "Test_Type"
Perf_Cell_Pre$Test_Type[Perf_Cell_Pre$Test_Type=="GINConv"] = "GraphDRP"
Perf_Drug_Pre$Test_Type[Perf_Drug_Pre$Test_Type=="GINConv"] = "GraphDRP"
Perf_Cell_Pre$Test_Type[Perf_Cell_Pre$Test_Type=="best_mse_paccmann_v2"] = "PaccMann_MSE"
Perf_Drug_Pre$Test_Type[Perf_Drug_Pre$Test_Type=="best_mse_paccmann_v2"] = "PaccMann_MSE"
Perf_Cell_Pre$Test_Type[Perf_Cell_Pre$Test_Type=="best_pearson_paccmann_v2"] = "PaccMann_PCC"
Perf_Drug_Pre$Test_Type[Perf_Drug_Pre$Test_Type=="best_pearson_paccmann_v2"] = "PaccMann_PCC"

# Annotate Cell TCGA Code
idx = match(Perf_Cell_Pre$Cell, Anno_Cells$SANGER_MODEL_ID)
Perf_Cell_Pre = Perf_Cell_Pre %>% mutate(Cell_TCGA=Anno_Cells$TCGA_CODE[idx])
Perf_Cell_Pre$Cell_TCGA[is.na(Perf_Cell_Pre$Cell_TCGA)] = "Unknown"

# Annotate Drug Target Pathway
idx = match(Perf_Drug_Pre$Drug, Anno_Drugs$Drug_CID)
Perf_Drug_Pre = Perf_Drug_Pre %>% mutate(Drug_Pathway=Anno_Drugs$Target_Pathway[idx])
Perf_Drug_Pre$Drug_Pathway[is.na(Perf_Drug_Pre$Drug_Pathway)] = "Unknown"

Perf_Cell_Pre$Test_Type = Perf_Cell_Pre$Model
Perf_Drug_Pre$Test_Type = Perf_Drug_Pre$Model


model = "PaccMann"
pretrain = "PaccMann_MSE"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_PaccMann_MSE = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre,
                                           model=model, pretrain=pretrain, dir=dir, save=T)

model = "PaccMann"
pretrain = "PaccMann_PCC"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_PaccMann_PCC = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre,
                                           model=model, pretrain=pretrain, dir=dir, save=T)

model = "GraphDRP"
pretrain = "GraphDRP"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_GraphDRP = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                       model=model, pretrain=pretrain, dir=dir, save=T)

model = "TGDRP"
pretrain = "TGDRP"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_TGDRP = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                    model=model, pretrain=pretrain, dir=dir, save=T)

model = "TGDRP"
pretrain = "TGDRP_Pre"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_TGDRP_Pre = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                        model=model, pretrain=pretrain, dir=dir, save=T)

model = "TGSA"
pretrain = "TGSA"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_TGSA = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                   model=model, pretrain=pretrain, dir=dir, save=T)

model = "TGSA"
pretrain = "TGSA_Pre"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_TGSA_Pre = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                       model=model, pretrain=pretrain, dir=dir, save=T)

model = "DRPreter"
pretrain = "DRPreter"
# dir = mkdir(sprintf("RMSE [%s, Pretrain & Train]", pretrain))
Test_DRPreter = compare_pretrain_multi(Perf_Cell, Perf_Cell_Pre, Perf_Drug, Perf_Drug_Pre, 
                                       model=model, pretrain=pretrain, dir=dir, save=T)


test_type = c("Normal", "Cell-Blind", "Drug-Blind", "Strict-Blind")
model_names = c("PaccMann_MSE", "PaccMann_PCC", "GraphDRP", 
                "TGDRP", "TGDRP_Pre", "TGSA", "TGSA_Pre", "DRPreter")

Test_Corr_Pre = list(Test_PaccMann_MSE, Test_PaccMann_PCC, Test_GraphDRP,
                     Test_TGDRP, Test_TGDRP_Pre, Test_TGSA, Test_TGSA_Pre, Test_DRPreter)

names(Test_Corr_Pre) = model_names

Test_Corr_Pre_ = Reduce(rbind, Test_Corr_Pre)
rownames(Test_Corr_Pre_) = NULL

Test_Corr_Pre_ = Test_Corr_Pre_ %>% 
  mutate(Model = rep(model_names, each=4), 
         Test_Type = rep(test_type, times=length(model_names))) %>% 
  mutate(Model = Model %>% factor(levels=model_names), 
         Test_Type = Test_Type %>% factor(levels=test_type)) %>% 
  relocate(Model, Test_Type, .before=everything()) %>% as.data.frame


Test_Corr_Cell = Test_Corr_Pre_ %>% 
  reshape2::acast(Model~Test_Type, value.var="Corr_Cell_RMSE") %>% 
  as.data.frame %>% relocate(all_of(test_type))

Test_Corr_Drug = Test_Corr_Pre_ %>% 
  reshape2::acast(Model~Test_Type, value.var="Corr_Drug_RMSE") %>% 
  as.data.frame %>% relocate(all_of(test_type))


main1 = "PCC of Cell RMSE [Pretrain-Train]"
main2 = "PCC of Drug RMSE [Pretrain-Train]"
color = scale_fill_gradient(low="white", high="firebrick3")
add = list(theme(legend.key.width=unit(0.8, "cm"), 
                 legend.key.height=unit(1, "cm")))

legend_c = bquote(PCC(RMSE[C]))
legend_d = bquote(PCC(RMSE[D]))

Test_Corr_Pre_ %>% grid_def(Model, Test_Type, fill=PCC_Cell_RMSE, main=main1, 
                            legend=legend_c, color=color, add=add,
                            round=2, size=7.2, axis_tx=20, legend_tl=20, legend_tx=20,
                            margin_lg=0.4, width=30, height=14.4, mean_summary=F, save=T, save_svg=T)

Test_Corr_Pre_ %>% grid_def(Model, Test_Type, fill=PCC_Drug_RMSE, main=main2, 
                            legend=legend_d, color=color, add=add,
                            round=2, size=7.2, axis_tx=20, legend_tl=20, legend_tx=20,
                            margin_lg=0.4, width=30, height=14.4, mean_summary=F, save=T, save_svg=T)


supplementary = T
if (supplementary) {
  # sheets = c("Cell", "Drug")
  # file = "RMSE Corr of Pretrain & Trained Model.xlsx"
  # write.xlsx(list(Test_Corr_Cell[model_names, ], 
  #                 Test_Corr_Drug[model_names, ]), 
  #            rowNames=T, sheetName=sheets, file=file)
  
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
  
  
  ### [Source Data] Supplementary Fig. 19
  Test_Corr_Pre_ = Test_Corr_Pre_ %>% 
    rename(PCC_Cell_RMSE=Corr_Cell_RMSE, PCC_Drug_RMSE=Corr_Drug_RMSE)
  Test_Corr_Pre_ %>% save_for_nc(num=19, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 20
  Perf_Pretrain_ = Perf_Pretrain %>% rename(Num_Test=N_Test)
  Perf_Pretrain_ %>% save_for_nc(num=20, suppl=T)
}


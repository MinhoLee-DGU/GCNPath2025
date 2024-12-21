#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

source("../utils/functions.R")
loadings()

# processed_data/cell_data/TCGA/TCGA_Drug_Info.csv
file = "TCGA_Drug_Info.csv"
TCGA_Drug_Info = read.csv(file)   # 16 x 8



##### 2-1. Open Prediction Files

read_pred = function(dir, pattern, sep=",") {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_.*([0-9]+$)\\.csv", "Pred_CCLE_.*([0-9]+)\\.csv"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      Pred_TP = fread(df_name, header=T, sep=sep)
      df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
      nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
      Pred_TP$Seed = nth
      Pred = Pred %>% rbind(Pred_TP)
    }
    return(Pred)
  }
}

read_pred_all = function(dir_list, model_name, pattern, sep=",") {
  
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
      
      Pred_TP$Model = model_name
      Pred_TP$Dataset = dataset
      Pred_TP$Test_Type = test_type
      Pred_TP = Pred_TP %>% mutate(Model=model_name, Dataset=dataset, Test_Type=test_type) %>% 
        relocate(Model, Dataset, Test_Type, Seed, .before=everything()) %>% as.data.frame
      Pred = Pred %>% rbind(Pred_TP)
    }
  }
  return(Pred)
}

test_type = "Normal"
pattern = "pred_tcga_seed([0-9]+).csv"
dir_list = sprintf("../results/IC50_GDSC/%s/RGCN", test_type)

Pred_GCNPath = read_pred_all(dir_list, "GCNPath", pattern)   # 4140 x 12
idx = match(Pred_GCNPath$Drug_CID, TCGA_Drug_Info$Drug_CID)

Pred_GCNPath = Pred_GCNPath %>% 
  mutate(Drug_Name=TCGA_Drug_Info$Drug_Name[idx], 
         GDSC=TCGA_Drug_Info$GDSC[idx])



##### 2-2. Responder vs Non-Responder [Boxplot]

boxplot_tcga = function(Pred, test_type="Normal", main=NULL, 
                        resp_class1=T, tumor_only=T, drug_gdsc=T, axis_tl=25, axis_tx=18, 
                        legend_tl=22.5, legend_tx=22.5, width=22.5, height=18, save=F) {
  
  normal_type = c("Solid Tissue Normal")
  resp_class = c("Complete Response", "Partial Response", 
                 "Stable Disease", "Clinical Progressive Disease")
  
  if (drug_gdsc) {
    Pred = Pred %>% subset(Test_Type %in% test_type) %>% 
      mutate(Drug_Name=sprintf("%s [%s]", Drug_Name, ifelse(GDSC, "O", "X")))
  }
  
  if (tumor_only) {
    Pred = Pred %>% subset(!(Sample_Type %in% normal_type))
  }
  
  if (resp_class1) {
    labs = c("Responder", "Non-Responder")
    Pred$Resp_Class = ifelse(Pred$Response %in% resp_class[1:2], labs[1], labs[2])
  } else {
    # labs = c("Complete\nResponse", "The Rest")
    labs = c("Complete Response", "The Rest")
    Pred$Resp_Class = ifelse(Pred$Response %in% resp_class[1], labs[1], labs[2])
  }
  
  ylab = bquote(Predicted~ln(IC[50]))
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
  Pred$Resp_Class = Pred$Resp_Class %>% factor(levels=labs)
  Pred$Drug_Name = Pred$Drug_Name %>% factor
  
  # add = list(scale_fill_manual(values=color))
  # add[[2]] = theme(legend.key.width=unit(0.8, "cm"),
  #                  legend.key.height=unit(1.0, "cm"))
  # 
  # pos_point = position_dodge(width=0.75)
  # pos = position_dodge2(width=0.75, preserve="single")
  # 
  # Pred %>% boxplot_def(Drug_Name, Prediction, fill=Resp_Class, main=main,
  #                      ylab=ylab, add=add, legend="Class", alpha=0.9, alpha_point=0.25,
  #                      vjust=1, hjust=1, margin=0.5, pos=pos, pos_point=pos_point,
  #                      axis_tl=axis_tl, axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx,
  #                      width=width, height=height, force_bold=F, save=save)
  
  Resp_Num = Pred %>% reshape2::acast(Drug_CID~Resp_Class, value.var="Resp_Class", fun.aggregate=length)
  drugs_except = rownames(Resp_Num)[rowSums(Resp_Num==0)>0]
  Pred = Pred %>% subset(!(Drug_CID %in% drugs_except))

  add = list(scale_fill_manual(values=color))
  add[[2]] = theme(legend.key.width=unit(0.8, "cm"),
                   legend.key.height=unit(1.0, "cm"))

  pos_point = position_dodge(width=0.7)
  pos = position_dodge2(width=0.8, preserve="single")
  pos_pwc = position_dodge2(width=0.8, preserve="single")

  margin1 = margin(l=0.25, r=0.25, unit="cm")
  margin2 = margin(0.25, 0.25, 0.15, 0.25, unit="cm")
  margin3 = margin(r=1, unit="cm")
  margin4 = margin(l=0.25, r=0.5, unit="cm")

  font1 = font(object="ylab", size=axis_tl, margin=margin1)
  font2 = font(object="axis.text", size=axis_tx, margin=margin2, color="grey30")
  font3 = font(object="legend.title", size=legend_tl, margin=margin3)
  font4 = font(object="legend.text", size=legend_tx, margin=margin4)
  font = font1 + font2 + font3 + font4

  pl = Pred %>%
    ggboxplot(x="Drug_Name", y="Prediction", fill="Resp_Class",
              outlier.shape=NA, size=0.5, alpha=0.9, xlab=F)

  pl$layers[[1]]$position = pos
  pl_point = geom_point(data=Pred, aes(group=Resp_Class), size=1.5, alpha=0.2, position=pos_point)

  pl = pl + pl_point +
    labs(y=ylab, fill="Response") + font +
    rotate_x_text(angle=30, hjust=1, vjust=1) +
    scale_fill_manual(labels=labs, values=color)

  # Comparison within drug [Predicted vs Actual]
  # Not know why set "greater", but it produced the correct plots...
  method_args = list(alternative="greater")
  pl = pl + stat_pwc(aes(group=Resp_Class), data=Pred, label="p.signif", 
                     label.size=5, ref.group=labs[2], method.args=method_args,
                     tip.length=0, bracket.nudge.y=0.02, step.increase=0, hide.ns=F)

  pl = pl %>% ggpar(legend="bottom") +
    theme(legend.key.size=unit(1, "cm"))

  if (save) {
    save_fig(pl, main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
}

dir = mkdir("Responder and Non-Responder [Boxplot]")
main = sprintf("%s/GCNPath [Drug & Effect_Size]", dir)
Pred_GCNPath %>% boxplot_tcga(test_type=test_type, main=main, resp_class1=T, save=T)



##### 2-3. Tumor vs Non-Tumor [Boxplot]

select_pat_tn = function(Pred) {
  
  col = c("Patient", "Sample", "Sample_Type", "Drug_CID", "TCGA_Code")
  TCGA_TN_Info = Pred[, col] %>% distinct %>% as.data.frame   # 414
  
  normal_type = c("Solid Tissue Normal")
  pat_normal = TCGA_TN_Info %>% subset(Sample_Type %in% normal_type) %>% pull(Patient)   # 13
  TCGA_TN_Info = TCGA_TN_Info %>% subset(Patient %in% pat_normal)   # 13
  TCGA_TN_Info$Tumor_Normal = ifelse(TCGA_TN_Info$Sample_Type %in% normal_type, "Normal", "Tumor")
  
  # Patients having both normal and tumor samples
  pat_group1 = TCGA_TN_Info %>% group_by(Patient) %>% 
    summarise(n=n_distinct(Tumor_Normal)) %>% subset(n==2) %>% pull(Patient)
  
  # Patients whose normal and tumor samples were from the same tissue
  pat_group2 = TCGA_TN_Info %>% group_by(Patient) %>% 
    summarise(n=n_distinct(TCGA_Code)) %>% subset(n==1) %>% pull(Patient)
  
  patient = intersect(pat_group1, pat_group2)
  Pred = Pred %>% subset(Patient %in% patient)
  sprintf("%s patients have both normal-tumor samples from the same tissue...", length(patient)) %>% print
  
  return(Pred)
}

boxplot_tn = function(Pred_TN, test_type="Normal", resp_class1=T, main=NULL, 
                      axis_tl=22.5, axis_tx=18, width=7.2, height=15, save=F) {
  
  normal_type = c("Solid Tissue Normal")
  resp_class = c("Complete Response", "Partial Response", 
                 "Stable Disease", "Clinical Progressive Disease")
  
  Pred_TN = Pred_TN %>% subset(Test_Type %in% test_type) %>% 
    mutate(Tumor_Normal=ifelse(Sample_Type %in% normal_type, "Normal", "Tumor"))
  
  if (resp_class1) {
    labs = c("Responder", "Non-Responder")
    Pred_TN$Resp_Class = ifelse(Pred_TN$Response %in% resp_class[1:2], labs[1], labs[2])
  } else {
    labs = c("Complete\nResponse", "The Rest")
    Pred_TN$Resp_Class = ifelse(Pred_TN$Response %in% resp_class[1], labs[1], labs[2])
  }
  
  # if ("Seed" %in% colnames(Pred_TN)) {
  #   Pred_TN = Pred_TN %>%
  #     group_by(Patient, Tumor_Normal, Seed, Resp_Class) %>%
  #     summarise(Pred_Mean=mean(Prediction)) %>% as.data.frame
  # } else {
  #   Pred_TN = Pred_TN %>%
  #     group_by(Patient, Tumor_Normal, Resp_Class) %>%
  #     summarise(Pred_Mean=mean(Prediction)) %>% as.data.frame
  # }
  
  Pred_TN = Pred_TN %>% subset(select=-c(Sample_Type, Sample)) %>% 
    reshape2::dcast(...~Tumor_Normal, value.var=c("Prediction")) %>% 
    as.data.frame %>% mutate(Diff_Pred=Tumor-Normal)
  
  ylab = "△Pred [Tumor-Normal]"
  color = c("royalblue3", "#66b3ed")
  # add = list(scale_fill_manual(values=color))
  Pred_TN$Resp_Class = Pred_TN$Resp_Class %>% factor(levels=labs)
  
  # Pred_TN %>% boxplot_def(Resp_Class, Diff_Pred, main=main, add=add, 
  #                         alpha=0.9, vjust=1.0, hjust=1.0, margin=0.4,
  #                         ylab=ylab, legend=F, axis_tl=axis_tl, axis_tx=axis_tx, 
  #                         force_bold=F, width=width, height=height, save=save)
  # 
  # Pred_Fine = Pred_TN %>% subset(Resp_Class==labs[1])
  # Pred_Ill = Pred_TN %>% subset(Resp_Class==labs[2])
  # stat_diff = wilcox.test(Pred_Fine$Diff_Pred, Pred_Ill$Diff_Pred, alternative="less")
  # 
  # tn_info = c(mean(Pred_Fine$Diff_Pred), mean(Pred_Ill$Diff_Pred), 
  #             median(Pred_Fine$Diff_Pred), median(Pred_Ill$Diff_Pred), stat_diff$p.value)
  # 
  # Pred_TN = list(Pred_TN=Pred_TN, TN_Info=tn_info)
  # return(Pred_TN)
  
  pos_point = position_dodge(width=0.8)
  pos = position_dodge2(width=0.8, preserve="single")
  
  margin1 = margin(5, 5, 5, 5, unit="pt")
  margin2 = margin(10, 10, 5, 10, unit="pt")
  
  font1 = font(object="ylab", size=axis_tl, margin=margin1)
  font2 = font(object="axis.text", size=axis_tx, margin=margin2, color="grey30")
  font = font1 + font2
  
  pl = Pred_TN %>%
    ggboxplot(x="Resp_Class", y="Diff_Pred", fill="Resp_Class",
              outlier.shape=NA, size=0.5, alpha=0.9, xlab=F)
  
  pl$layers[[1]]$position = pos
  pl_point = geom_point(size=2, alpha=0.25, position=pos_point)
  
  pl = pl + pl_point + labs(y=ylab) + font +
    rotate_x_text(angle=30, hjust=1, vjust=1) +
    scale_fill_manual(labels=labs, values=color)
  
  method_args = list(alternative="greater")
  pl = pl + stat_pwc(aes(group=Resp_Class), label="p.signif", label.size=5, 
                     ref.group=labs[2], method.args=method_args,
                     tip.length=0, bracket.nudge.y=0.1, hide.ns=F)
  
  pl = pl %>% ggpar(legend="None")
  
  if (save) {
    save_fig(pl, main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
  
  return(Pred_TN)
}


Pred_GCNPath_TN = Pred_GCNPath %>% select_pat_tn

dir = mkdir("Tumor and Normal [Boxplot]")
main = sprintf("%s/GCNPath [Tumor & Normal]", dir)
TN_Info = Pred_GCNPath_TN %>% boxplot_tn(test_type=test_type, resp_class1=T, main=main, save=T)

diff_resp = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Responder"]
diff_non = TN_Info$Diff_Pred[TN_Info$Resp_Class=="Non-Responder"]
wilcox.test(diff_resp, diff_non, alternative="less")$p.value

file = sprintf("%s/GCNPath_TN.csv", dir)
write.csv(Pred_GCNPath_TN, file=file, row.names=F)

file = sprintf("%s/GCNPath_TN.RData", dir)
save(TN_Info, file=file)



##### 2-4. Responder vs Non-Responder [Significance Test, Scatter-Plot]

sig_test = function(Pred, resp_class1=T, tumor_only=T) {
  
  normal_type = c("Solid Tissue Normal")
  resp_class = c("Complete Response", "Partial Response", 
                 "Stable Disease", "Clinical Progressive Disease")
  
  tryCatch({
    if (resp_class1) {
      labs = c("Responder", "Non-Responder")
      Pred$Resp_Class = ifelse(Pred$Response %in% resp_class[1:2], labs[1], labs[2])
    } else {
      labs = c("Complete\nResponse", "The Rest")
      Pred$Resp_Class = ifelse(Pred$Response %in% resp_class[1], labs[1], labs[2])
    }
    
    if (tumor_only) {
      Pred = Pred %>% subset(!(Sample_Type %in% normal_type))
    }
    
    Pred_Resp_O = Pred %>% subset(Resp_Class==labs[1])
    Pred_Resp_X = Pred %>% subset(Resp_Class==labs[2])
    
    mean_pos = ifelse(nrow(Pred_Resp_O)!=0, mean(Pred_Resp_O$Prediction), NA)
    mean_neg = ifelse(nrow(Pred_Resp_X)!=0, mean(Pred_Resp_X$Prediction), NA)
    effect_size = mean_pos - mean_neg
    
    stat_resp = wilcox.test(Pred_Resp_O$Prediction, Pred_Resp_X$Prediction, alternative="less")
    stat_info = c(stat_resp$alternative, mean_pos, mean_neg, 
                  effect_size, stat_resp$statistic, stat_resp$p.value)
    
    return(stat_info)
  }, error = function(e) {
    return_error = c("less", mean_pos, mean_neg, rep(NA, 3))
    return(return_error)
  })
}

sig_test_all = function(Pred, resp_class1=T, tumor_only=T) {
  
  Sig_Test = data.frame()
  drug_list = Pred$Drug_Name %>% unique
  test_type = Pred$Test_Type %>% unique
  
  for (test_type_ in test_type) {
    for (drug in drug_list) {
      try({
        Pred_Temp = Pred %>% subset(Drug_Name==drug & Test_Type==test_type_)
        stat_info = Pred_Temp %>% sig_test(resp_class1=resp_class1, tumor_only=tumor_only)
        Sig_Test =  Sig_Test %>% rbind(c(drug, test_type_, stat_info, as.numeric(stat_info[6])<0.05))
      })
    }
  }
  
  colnames(Sig_Test) = c("Drug_Name", "Test_Type", "Hypothesis", 
                         "Pred_Mean_Pos", "Pred_Mean_Neg", 
                         "Effect_Size", "W", "Pval", "Signif_Pval")
  
  Sig_Test = Sig_Test %>% 
    mutate_at(colnames(Sig_Test)[4:8], as.numeric) %>% 
    arrange(Test_Type, Effect_Size) %>% as.data.frame
  
  Sig_Test = Sig_Test %>% 
    mutate(MLog10_Pval=-log10(Pval)) %>% 
    relocate(MLog10_Pval, .after=Pval) %>% as.data.frame
  
  return(Sig_Test)
}

plot_effect_pval = function(Pred, main=NULL, axis_tl=32, axis_tx=27, 
                            legend_tl=20, legend_tx=18, inner_tx=6, 
                            width=18.6, height=16, drug_gdsc=T, save=F) {
  
  if (drug_gdsc) {
    Pred = Pred %>% mutate(GDSC = ifelse(GDSC, "O", "X")) %>% 
      mutate(Drug_Name = sprintf("%s [%s]", Drug_Name, GDSC))
  }
  
  line_h = geom_hline(aes(yintercept=0), col="red", lty=2)
  line_v = geom_vline(aes(xintercept=-log(0.05, 10)), col="red", lty=2)
  color = scale_color_gradient2(low="beige", mid="yellow", high="red")
  repel = text_repel_def(Pred, label=Drug_Name, fontface="plain", size=inner_tx, force=4)
  
  ylab = "Effect Size"
  xlab = bquote("-Log"[10](Pval))
  add_list = list(line_h, line_v, repel)
  # add_list = list(color, line_h, line_v, repel)
  
  Pred %>% plot_def(MLog10_Pval, Effect_Size, size=Freq, alpha=0.5,
                    main=main, xlab=xlab, ylab=ylab, margin=0.4, 
                    axis_tl=axis_tl, axis_tx=axis_tx, legend_tl=legend_tl, legend_tx=legend_tx,
                    add=add_list, force_bold=F, width=width, height=height, save=save)
  
  # pl = Pred %>% ggscatter("MLog10_Pval", "Effect_Size", color="Freq", size="Freq", 
  #                         xlab="-Log10(Pval)", ylab="Effect Size", legend="right")
}

Sig_GCNPath = sig_test_all(Pred_GCNPath, resp_class1=T, tumor_only=T)
idx = match(Sig_GCNPath$Drug_Name, TCGA_Drug_Info$Drug_Name)

Sig_GCNPath = Sig_GCNPath %>% 
  mutate(Freq=TCGA_Drug_Info$Total[idx], 
         GDSC=TCGA_Drug_Info$GDSC[idx])

dir = mkdir("Performance Significance [Scatter Plot]")
main = sprintf("%s/GCNPath [Pvalue & Effect_Size]", dir)
Sig_GCNPath %>% plot_effect_pval(main=main, save=T)

col = c("Drug_Name", "Pval", "Effect_Size")
Sig_GCNPath[, col] %>% subset(Effect_Size<0 & Pval<0.05)
Sig_GCNPath[, col] %>% subset(Effect_Size<0 & Pval<0.05) %>% arrange(Pval)

# Drug_Name    Pval         Effect_Size
# Fluorouracil 6.305972e-15 -0.57408782
# Temozolomide 1.327700e-07 -0.07834665
# Sorafenib    4.902823e-04 -0.21196010
# Cisplatin    1.173591e-02 -0.15045287
# Gemcitabine  1.506262e-02 -0.22940618

file = "Prediction [TCGA].csv"
write.csv(Pred_GCNPath, file=file, row.names=F)



##### 2-5. Conformation of batch correction
# Does GSVA virtually do batch correction?

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

plot_def_ <- function(df, x, y, main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, legend=NULL, 
                     text=NULL, add=NULL, color="black", shape=NULL, color_line="red", size=1, alpha=1, stroke=1,
                     plot_tl=22.5, axis_tl=18, axis_tx=15, legend_tl=16.5, legend_tx=15, 
                     margin=0.5, margin_pl=0.25, margin_lg=0.4, width=15, height=15, dpi=400, text_ratio=1,
                     plot_face="plain", axis_face="plain", legend_face="plain", raster_dev="ragg_png", pos="identity", 
                     pos_legend=NULL, show_title=F, xy_line=F, force_bold=F, raster=F, save=F, save_svg=T, ...) {
  
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  margin_x = margin(t=margin, b=margin, unit="cm")
  margin_y = margin(l=margin, r=margin, unit="cm")
  
  if (force_bold) {
    plot_face = "bold"; axis_face = "bold"; legend_face = "bold" 
  }
  
  e_text = element_text
  pl <- ggscatter(df, x = x, y = y, 
                  color = color, shape = shape, size = size, alpha = alpha, 
                  xlab = xlab, ylab = ylab, ...) + theme_classic() +
    theme(
      plot.title = e_text(size = plot_tl, face = plot_face, hjust = 0.5),
      axis.title = e_text(size = axis_tl, face = axis_face, hjust = 0.5), 
      axis.text.x = e_text(size = axis_tx, face = axis_face, margin = margin_x), 
      axis.text.y = e_text(size = axis_tx, face = axis_face, margin = margin_y),
      legend.title = e_text(size = legend_tl, face = legend_face),
      legend.text = e_text(size = legend_tx, face = legend_face)
    )
  
  if (!is.null(xlim)) pl <- pl + xlim(xlim)
  if (!is.null(ylim)) pl <- pl + ylim(ylim)
  if (xy_line) pl <- pl + geom_abline(slope = 1, intercept = 0, color = color_line, lty = 2)
  
  if (!is.null(add)) for (i in 1:length(add)) pl <- pl + add[[i]]
  if (!is.null(legend) && legend==F) pl <- pl + theme(legend.position="none")
  if (!is.null(pos_legend)) pl <- pl + theme(legend.position = pos_legend)
  
  # Save the plot if requested
  if (save) {
    pl %>% save_fig_ggpubr(main=main, width=width, height=height, svg=save_svg)
  } else print(pl)
}

plot_pca_batch = function(Omics_Source, Omics_Target, tcga_source, tcga_target, 
                          dir=NULL, db_source="SANGER", db_target="TCGA", 
                          tcga_except="UNCLASSIFIED", size=1, size_s=2, size_t=1, alpha=0.5,
                          stroke=1.2, tcga_topk=7, axis_tl=24, axis_tx=20,  
                          width=16, height=12, add_db=NULL, add_tcga=NULL, save=T) {
  
  # Annotation Label [TCGA Code]
  tcga_code = c(tcga_source, tcga_target)
  cond1 = length(tcga_source)==nrow(Omics_Source)
  cond2 = length(tcga_target)==nrow(Omics_Target)
  if (!cond1 | !cond2) stop("Check the TCGA Codes...")
  
  # Annotation Label [DB]
  db1 = rep(db_source, nrow(Omics_Source))
  db2 = rep(db_target, nrow(Omics_Target))
  db = c(db1, db2)
  
  # Compress Omics Data by PCA
  db_lvl = c(db_source, db_target)
  cond = identical(colnames(Omics_Source), colnames(Omics_Target))   # T
  if (!cond) Omics_Source = Omics_Source[, colnames(Omics_Target)]
  Omics_PCA = rbind(Omics_Source, Omics_Target) %>% 
    pca_exp(db=db, tcga_code=tcga_code, db_lvl=db_lvl)
  
  # Select most frequent TCGA Codes from target omics data
  # Except TCGA Codes not existing or UNCLASSIFIED in source omics data
  tcga_source_ = tcga_source[!(tcga_source %in% tcga_except)]
  tcga_target_ = tcga_target[tcga_target %in% tcga_source_]
  
  tcga_main = tcga_target_
  cond_topk = !is.null(tcga_topk) && tcga_topk!=0
  
  if (cond_topk) {
    tcga_main = tcga_main %>% table %>% sort(decreasing=T) %>% head(tcga_topk) %>% names
  } else {
    tcga_main = tcga_main %>% table %>% sort(decreasing=T) %>% names
  }
  
  color_tcga = add_tcga[[1]]$palette(1)
  if (length(color_tcga)!=length(tcga_main)) {
    add_tcga[1:2] = NULL
  }
  
  # Set X and Y Ranges
  xlim = Omics_PCA$PC1 %>% range
  ylim = Omics_PCA$PC2 %>% range
  xlim = c(xlim[1]-0.1, xlim[2]+0.1)
  ylim = c(ylim[1]-0.1, ylim[2]+0.1)
  
  # Plot PCA [Source & Target]
  main = sprintf("%s/PC2 Plot", dir)
  Omics_PCA %>% arrange(Database) %>%
    plot_def_("PC1", "PC2", color="Database", shape="Database", main=main,
              legend="Database", xlim=xlim, ylim=ylim, size=size, stroke=stroke,
              alpha=alpha, add=add_db, axis_tl=axis_tl, axis_tx=axis_tx,  
              width=width, height=height, save=save, save_svg=save)
  
  # Plot PCA [Source]
  main = sprintf("%s/PC2 Plot [%s, Top %s]", dir, db_source, tcga_topk)
  if (!cond_topk) main = sprintf("%s/PC2 Plot [%s]", dir, db_source)
  
  Omics_PCA %>%
    subset(Database==db_source & TCGA_Code %in% tcga_main) %>%
    mutate(`TCGA Code`=TCGA_Code %>% factor(levels=tcga_main)) %>%
    plot_def_("PC1", "PC2", color="TCGA Code", shape="TCGA Code",
              main=main, xlim=xlim, ylim=ylim, size=size_s, stroke=stroke,
              alpha=alpha, add=add_tcga, axis_tl=axis_tl, axis_tx=axis_tx,  
              width=width, height=height, save=save, save_svg=save)
  
  # Plot PCA [Target]
  main = sprintf("%s/PC2 Plot [%s, Top %s]", dir, db_target, tcga_topk)
  if (!cond_topk) main = sprintf("%s/PC2 Plot [%s]", dir, db_target)
  
  Omics_PCA %>%
    subset(Database==db_target & TCGA_Code %in% tcga_main) %>%
    mutate(`TCGA Code`=TCGA_Code %>% factor(levels=tcga_main)) %>%
    plot_def_("PC1", "PC2", color="TCGA Code", shape="TCGA Code", 
              main=main, xlim=xlim, ylim=ylim, size=size_t, stroke=stroke,
              alpha=alpha, add=add_tcga, axis_tl=axis_tl, axis_tx=axis_tx,  
              width=width, height=height, save=save, save_svg=save)
  
  return(Omics_PCA)
}

dir = "../processed/cell_data_biocarta"
file1 = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
file2 = sprintf("%s/TCGA_RNA_GSVA.csv", dir)

SANGER_RNA_GSVA = fread_def(file1)
TCGA_RNA_GSVA = fread_def(file2)

dir = "../data/cell_data"
file1 = sprintf("%s/SANGER_RNA_TPM.csv", dir)
file2 = sprintf("%s/TCGA_RNA_TPM.csv", dir)

SANGER_RNA_TPM = fread_def(file1, col_numeric=T)
TCGA_RNA_TPM = fread_def(file2, col_numeric=T)

suppressMessages(library(cogena))
identical(rownames(SANGER_RNA_TPM), rownames(SANGER_RNA_GSVA))   # T
identical(rownames(TCGA_RNA_TPM), rownames(TCGA_RNA_GSVA))       # T

dir = "../data/path_data"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.entrez.gmt", dir)
Path_List = gmt2list(file)
geneset = Path_List %>% unlist %>% unique   # 1509
geneset_ = Reduce(intersect, list(geneset, colnames(SANGER_RNA_TPM), colnames(TCGA_RNA_TPM)))   # 1503


# Cell Annotation File from
# GDSC_Last/processed/cell_data/GDSC/Anno_Cells.csv
file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)

# TCGA Annotation File from
# GDSC_Last/processed_data/cell_data/TCGA/TCGA_Info.csv
file = "TCGA_Info.csv"
TCGA_Info = read.csv(file)


# Annotation Label [TCGA Code from SANGER]
idx = match(rownames(SANGER_RNA_GSVA), Anno_Cells$SANGER_MODEL_ID)
tcga_code1 = Anno_Cells$TCGA_CODE[idx]
tcga_code1 = ifelse(!is.na(tcga_code1), tcga_code1, "UNCLASSIFIED")

# Annotation Label [TCGA Code from TCGA]
all(rownames(TCGA_RNA_GSVA) %in% TCGA_Info$cases)   # T
idx = match(rownames(TCGA_RNA_GSVA), TCGA_Info$cases)
tcga_code2 = gsub("^TCGA-", "", TCGA_Info$project[idx])

# # GSVA, All Samples
# dir = mkdir("Batch Correction [GSVA, Total Patients]")
# plot_pca_batch(SANGER_RNA_GSVA, TCGA_RNA_GSVA, tcga_code1, tcga_code2, 
#                dir=dir, db_source="SANGER", db_target="TCGA", tcga_topk=5, save=T)

samples = Pred_GCNPath$Sample %>% unique   # 414
idx_tcga = rownames(TCGA_RNA_GSVA) %in% samples
tcga_code2[idx_tcga] %>% unique %>% length   # 26
tcga_code2[idx_tcga] %>% table %>% sort(decreasing=T) %>% head(10)
# LGG CESC HNSC STAD PAAD LIHC COAD PRAD SKCM ESCA 
# 87  59   47   47   41   17   15   14   12   11

idx_sanger = tcga_code1 %in% tcga_code2[idx_tcga]   # 546/1431
tcga_code1[idx_sanger] %>% unique %>% length   # 21
tcga_code1[idx_sanger] %>% table %>% sort(decreasing=T) %>% head(10)
# LUAD SKCM BRCA HNSC PAAD ESCA GBM STAD BLCA LUSC 
# 73   57   52   42   42   35   33  32   26   23 

shape_db = c(21, 24)
color_db = c("brown1", "royalblue1")
names(color_db) = c("SANGER", "TCGA")

shape_tcga = 0:6
color_tcga = c("brown1", "coral", "gold", "seagreen3", "royalblue1", "mediumorchid", "black")

add_db = list(scale_color_manual(values=color_db), 
              scale_shape_manual(values=shape_db), 
              theme(legend.key.size=unit(1, "cm")), 
              guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
 
add_tcga = list(scale_color_manual(values=color_tcga), 
                scale_shape_manual(values=shape_tcga), 
                theme(legend.key.size=unit(1, "cm")), 
                guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))


# GSVA, Subset of Samples
dir = mkdir("Batch Correction [GSVA]")
PCA_GSVA = plot_pca_batch(SANGER_RNA_GSVA[idx_sanger, ], 
                          TCGA_RNA_GSVA[idx_tcga, ], 
                          tcga_code1[idx_sanger], tcga_code2[idx_tcga], 
                          dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
                          tcga_topk=0, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)

PCA_GSVA = plot_pca_batch(SANGER_RNA_GSVA[idx_sanger, ], 
                          TCGA_RNA_GSVA[idx_tcga, ], 
                          tcga_code1[idx_sanger], tcga_code2[idx_tcga], 
                          dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
                          tcga_topk=7, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)

# Raw, Subset of Samples
dir = mkdir("Batch Correction [Raw]")
PCA_Raw = plot_pca_batch(SANGER_RNA_TPM[idx_sanger, geneset_], 
                         TCGA_RNA_TPM[idx_tcga, geneset_], 
                         tcga_code1[idx_sanger], tcga_code2[idx_tcga], 
                         dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
                         tcga_topk=0, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)

PCA_Raw = plot_pca_batch(SANGER_RNA_TPM[idx_sanger, geneset_], 
                         TCGA_RNA_TPM[idx_tcga, geneset_], 
                         tcga_code1[idx_sanger], tcga_code2[idx_tcga], 
                         dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
                         tcga_topk=7, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)

# # Standardization, Subset of Samples
# SANGER_RNA_TPM_ = SANGER_RNA_TPM[, geneset_] %>% scale(center=T, scale=T) %>% as.data.frame
# TCGA_RNA_TPM_ = TCGA_RNA_TPM[, geneset_] %>% scale(center=T, scale=T) %>% as.data.frame
# TCGA_RNA_TPM_[is.na(TCGA_RNA_TPM_)] = 0
# 
# dir = mkdir("Batch Correction [Raw & Standardization]")
# plot_pca_batch(SANGER_RNA_TPM_, TCGA_RNA_TPM_[idx_tcga, ], tcga_code1, tcga_code2[idx_tcga],
#                dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
#                tcga_topk=0, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)
# 
# plot_pca_batch(SANGER_RNA_TPM_, TCGA_RNA_TPM_[idx_tcga, ], tcga_code1, tcga_code2[idx_tcga],
#                dir=dir, db_source="SANGER", db_target="TCGA", add_db=add_db, add_tcga=add_tcga,
#                tcga_topk=7, alpha=0.64, size=1, size_s=2.5, size_t=2.5, width=16, height=12, save=T)



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

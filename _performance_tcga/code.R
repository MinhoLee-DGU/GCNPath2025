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

Pred_GCNPath = read_pred_all(dir_list, "GCNPath", pattern)   # 22770 x 12
idx = match(Pred_GCNPath$Drug_CID, TCGA_Drug_Info$Drug_CID)

Pred_GCNPath = Pred_GCNPath %>% 
  mutate(Drug_Name=TCGA_Drug_Info$Drug_Name[idx], 
         GDSC=TCGA_Drug_Info$GDSC[idx])



##### 2-2. Responder vs Non-Responder [Boxplot]

boxplot_tcga = function(Pred, test_type="Normal", main=NULL, 
                        resp_class1=T, tumor_only=T, axis_tl=25, axis_tx=18, 
                        legend_tl=22.5, legend_tx=22.5, width=22.5, height=15, save=F) {
  
  normal_type = c("Solid Tissue Normal")
  resp_class = c("Complete Response", "Partial Response", 
                 "Stable Disease", "Clinical Progressive Disease")
  
  Pred = Pred %>% subset(Test_Type %in% test_type) %>% 
    mutate(Drug_Name=sprintf("%s [%s]", Drug_Name, ifelse(GDSC, "O", "X")))
  
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

  pl = pl %>% ggpar(legend="bottom")

  if (save) {
    save_fig(pl, main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
}

dir = mkdir("Responder and Non-Responder [Boxplot]")
main = sprintf("%s/GCNPath [Drug & Effect_Size]", dir)
Pred_GCNPath %>% boxplot_tcga(test_type=test_type, resp_class1=T, main=main, save=T)



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
                      axis_tl=22.5, axis_tx=18, width=9, height=15, save=F) {
  
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
  color = c("firebrick2", "royalblue2")
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
  pl = pl + stat_pwc(aes(group=Resp_Class), label="p.format", label.size=5, 
                     ref.group=labs[2], method.args=method_args,
                     tip.length=0, bracket.nudge.y=0.1, hide.ns=F)
  
  pl = pl %>% ggpar(legend="None")
  
  if (save) {
    save_fig(pl, main, width=width, height=height, units="cm", svg=T)
  } else print(pl)
  
  return(Pred_TN)
}

dir = mkdir("Tumor and Normal [Boxplot]")
Pred_GCNPath_TN = Pred_GCNPath %>% select_pat_tn

main = sprintf("%s/GCNPath [Tumor & Normal]", dir)
TN_Info = Pred_GCNPath_TN %>% boxplot_tn(test_type=test_type, resp_class1=T, main=main, save=T)

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
                            width=18.6, height=16, save=F) {
  
  Pred = Pred %>% mutate(GDSC = ifelse(GDSC, "O", "X")) %>% 
    mutate(Drug_Name = sprintf("%s [%s]", Drug_Name, GDSC))
  
  line_h = geom_hline(aes(yintercept=0), col="red", lty=2)
  line_v = geom_vline(aes(xintercept=-log(0.05, 10)), col="red", lty=2)
  color = scale_color_gradient2(low="beige", mid="yellow", high="red")
  repel = text_repel_def(Pred, label=Drug_Name, fontface="plain", size=inner_tx, force=4)
  
  ylab = "Effect Size"
  xlab = bquote("-Log"[10](Pval))
  # xlab = bquote(bold("-Log"["10"]~"(Pval)"))
  add_list = list(color, line_h, line_v, repel)
  
  Pred %>% plot_def(MLog10_Pval, Effect_Size, color=Freq, size=Freq, 
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
Sig_GCNPath %>% summarise(Num_Drug=sum(Effect_Size<0 & Pval<0.05, na.rm=T))

file = "Prediction [TCGA].csv"
write.csv(Pred_GCNPath, file=file, row.names=F)


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
  
  
  ### [Source Data] Fig. 5
  Pred_GCNPath_ = Pred_GCNPath %>% 
    subset(select=-c(Model, Dataset, Test_Type)) %>% 
    rename(Train_Seed=Seed, GDSC_Drug=GDSC) %>% 
    relocate(Train_Seed, .after=Prediction) %>% as.data.frame
  
  idx = match(Sig_GCNPath$Drug_Name, TCGA_Drug_Info$Drug_Name)
  
  Sig_GCNPath_ = Sig_GCNPath %>% 
    rename(GDSC_Drug=GDSC, Num_Total=Freq, 
           Minus_Log10_Pval=MLog10_Pval) %>% 
    mutate(Drug_CID=TCGA_Drug_Info$Drug_CID[idx], 
           Num_Pos=TCGA_Drug_Info$Resp_C1[idx], 
           Num_Neg=TCGA_Drug_Info$Non_C1[idx]) %>% 
    relocate(Drug_CID, GDSC_Drug, .after=Drug_Name) %>% 
    relocate(Num_Total, Num_Pos, Num_Neg, .after=Test_Type)
  
  TN_Info_ = TN_Info %>% 
    subset(select=-c(Model, Dataset, Test_Type)) %>% 
    rename(Pred_Normal=Normal, Pred_Tumor=Tumor, 
           Train_Seed=Seed, Response_Class=Resp_Class)
  
  Pred_TCGA = list(Pred_GCNPath_, Sig_GCNPath_, TN_Info_)
  Pred_TCGA %>% save_for_nc(num=5, suppl=F)
  
  
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

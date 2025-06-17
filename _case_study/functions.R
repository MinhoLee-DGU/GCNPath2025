
read_pred = function(dir, pattern, sep=",") {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_.*([0-9]+$)\\.csv", "Pred_CCLE_.*([0-9]+)\\.csv"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      try({
        Pred_TP = fread(df_name, header=T, sep=sep)
        # Pred_TP = read.table(df_name, header=T, sep=sep)
        df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
        nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
        model = strsplit(dir, "/")[[1]] %>% tail(1)
        
        # Pred_TP$Fold = nth
        Pred_TP$Seed = nth
        Pred_TP$Model = model
        Pred = Pred %>% rbind(Pred_TP)
      })
    }
    
    sprintf("# Predictions : %s", nrow(Pred)) %>% print
    return(Pred)
  }
}

subtype_def = function(Anno_Cells, category=NULL, 
                       col_category="CANCER_TYPE_DETAIL", 
                       col_cell="SANGER_MODEL_ID", names_rest="The Rest") {
  
  Anno_Cells$SubType_List = ifelse(
    Anno_Cells[[col_category]] %in% category, 
    Anno_Cells[[col_category]], 
    names_rest
  )
  
  SubType_List = split(Anno_Cells[[col_cell]], Anno_Cells$SubType_List)
  return(SubType_List)
}

filter_ic50 = function(Anno_Cells, filter=T) {
  if (filter) {
    Anno_Cells %>% subset(IC50)
  } else Anno_Cells
}

filter_non_can = function(Anno_Cells, filter=T) {
  if (filter) {
    Anno_Cells %>% subset(CANCER_TYPE!="Non-Cancerous")
  } else Anno_Cells
}

lookup_by_tissue = function(Anno_Cells, tissue, group="CANCER_TYPE", 
                            filt_ic50=F, filt_nc=T, num_min=3, return_num=T) {
  Anno_Cells = Anno_Cells %>%
    filter_ic50(filter=filt_ic50) %>% 
    filter_non_can(filter=filt_nc) %>% 
    subset(RNA & TISSUE %in% tissue) 
  
  Anno_Cells_ = Anno_Cells %>%
    group_by(across(all_of(group))) %>% summarise(Num=n()) %>%
    arrange(desc(Num)) %>% subset(Num>=num_min) %>% as.data.frame
    
  if (return_num) {
    return(Anno_Cells_)
  } else {
    Anno_Cells = Anno_Cells %>% filter(!!sym(group) %in% Anno_Cells_[[group]])
    return(Anno_Cells)
  }
}

get_anno_cd = function(Pred_Test, Anno_Cells, Anno_Drugs, 
                       filt_ic50=F, filt_nc=F, num_min=0) {
  
  Anno_Cells_Sub = Anno_Cells %>% 
    filter_ic50(filter=filt_ic50) %>% 
    filter_non_can(filter=filt_nc) %>% 
    subset(RNA) %>% 
    group_by(CANCER_TYPE_DETAIL) %>% 
    filter(n()>=num_min) %>% as.data.frame
  
  Pred_Test_Sub = Pred_Test %>% subset(Cell %in% Anno_Cells_Sub$SANGER_MODEL_ID)
  idx1 = match(Pred_Test_Sub$Cell, Anno_Cells$SANGER_MODEL_ID)
  idx2 = match(Pred_Test_Sub$Drug, Anno_Drugs$Drug_CID)
  
  Pred_Test_Sub = Pred_Test_Sub %>% 
    mutate(Cell_Name=Anno_Cells$MODEL_NAME[idx1],
           TCGA_Code=Anno_Cells$TCGA_CODE[idx1], 
           Tissue=Anno_Cells$TISSUE[idx1], 
           Cancer_Type=Anno_Cells$CANCER_TYPE[idx1], 
           Cancer_Detail=Anno_Cells$CANCER_TYPE_DETAIL[idx1], 
           Drug_Name=Anno_Drugs$Name[idx2],
           Drug_Pathway=Anno_Drugs$Target_Pathway[idx2])
  
  return(Pred_Test_Sub)
}

examine_subtype = function(Pred_Test, Anno_Cells, Anno_Drugs, cancer_type, 
                           filt_ic50=F, filt_nc=F, num_min=3) {
  
  min_max_diff = function(x) max(x)-min(x)
  
  Anno_Cells_Sub = Anno_Cells %>% 
    filter_ic50(filter=filt_ic50) %>% 
    filter_non_can(filter=filt_nc) %>% 
    subset(RNA & CANCER_TYPE %in% cancer_type) %>% 
    group_by(CANCER_TYPE_DETAIL) %>% 
    filter(n()>=num_min) %>% as.data.frame
  
  Pred_Test_Sub = Pred_Test %>% subset(Cell %in% Anno_Cells_Sub$SANGER_MODEL_ID)
  idx1 = match(Pred_Test_Sub$Cell, Anno_Cells$SANGER_MODEL_ID)
  idx2 = match(Pred_Test_Sub$Drug, Anno_Drugs$Drug_CID)
  
  Pred_Test_Sub = Pred_Test_Sub %>% 
    mutate(TCGA_Code=Anno_Cells$TCGA_CODE[idx1], 
           Cancer_Type=Anno_Cells$CANCER_TYPE[idx1], 
           Cancer_Detail=Anno_Cells$CANCER_TYPE_DETAIL[idx1], 
           Drug_Name=Anno_Drugs$Name[idx2],
           Drug_Pathway=Anno_Drugs$Target_Pathway[idx2])
  
  # Summarize by Drug_Pathway
  Pred_Median1 = Pred_Test_Sub %>% 
    group_by(Drug_Pathway) %>% 
    summarise(Median=median(Prediction)) %>% 
    arrange(Median) %>% as.data.frame
  
  Pred_Median_SubT1 = Pred_Test_Sub %>% 
    group_by(Cancer_Detail, Drug_Pathway) %>% 
    summarise(Median=median(Prediction), Num=n()) %>% 
    arrange(desc(Num)) %>% as.data.frame
  
  Pred_Median_Diff1 = Pred_Median_SubT1 %>% 
    group_by(Drug_Pathway) %>% 
    summarise(Median_Diff=min_max_diff(Median), 
              Median_Min=min(Median), Median_Max=max(Median)) %>% 
    arrange(desc(Median_Diff)) %>% as.data.frame
  
  # Summarize by Drug_Name
  Pred_Median2 = Pred_Test_Sub %>% 
    group_by(Drug_Name, Drug_Pathway) %>% 
    summarise(Median=median(Prediction)) %>% 
    arrange(Median) %>% as.data.frame
  
  Pred_Median_SubT2 = Pred_Test_Sub %>% 
    group_by(Cancer_Detail, Drug_Name) %>% 
    summarise(Median=median(Prediction), Num=n()) %>% 
    arrange(desc(Num)) %>% as.data.frame
  
  Pred_Median_Diff2 = Pred_Median_SubT2 %>% 
    group_by(Drug_Name) %>% 
    summarise(Median_Diff=min_max_diff(Median), 
              Median_Min=min(Median), Median_Max=max(Median)) %>% 
    arrange(desc(Median_Diff)) %>% as.data.frame
  
  Pred_Pathway = left_join(Pred_Median1, Pred_Median_Diff1, by="Drug_Pathway")
  Pred_Drug = left_join(Pred_Median2, Pred_Median_Diff2, by="Drug_Name")
  
  Pred_Drug = Pred_Drug %>% 
    arrange(desc(Median_Diff)) %>% 
    mutate(Top_Median = in_bottom(Median, 0.1),
           Top_Median_Diff = in_top(Median_Diff, 0.1))
  
  Pred_Pathway = Pred_Pathway %>% 
    arrange(desc(Median_Diff)) %>% 
    mutate(Top_Median = in_bottom(Median, 5), 
           Top_Median_Diff = in_top(Median_Diff, 5))
  
  Pred_Sum = list(Pred=Pred_Test_Sub, Pred_Drug=Pred_Drug, Pred_Pathway=Pred_Pathway)
  
  return(Pred_Sum)
}

in_rank = function(x, n, top=T) {
  if (n<0) stop("Parameter prob not adequete...")
  if (n>=length(x)) n = length(x)
  if (n>=0 & n<1) n = floor(length(x)*n) 
  cond = ifelse_def(top, x>=sort(x, decreasing=T)[n], x<=sort(x)[n])
  if (n==0) cond = rep(F, length(x))
  return(cond)
}

in_top = function(x, n) {
  in_rank(x, n, top=T)
}

in_bottom = function(x, n) {
  in_rank(x, n, top=F)
}

append_def = function(A, ...) {
  add_list = list(...)
  A = c(A, add_list)
  return(A)
}

plot_with_tpm = function(Pred, TPM, gene, marker=NULL, xlab=NULL, 
                         SubType_List=NULL, trend_line=T, anno_subtype=F, 
                         exclude_else=T, name_else="The Others", add=NULL, 
                         cells=NULL, col_cell="Cell_Name", width=18, height=18, 
                         axis_tl=30, axis_tx=22.5, legend_tl=20, legend_tx=18, ...) {
  
  if (is.null(add)) add = list()
  idx = match(Pred$Cell, rownames(TPM))
  Pred$Gene = TPM[idx, gene]
  Pred$Gene_Name = gene
  
  if (is.null(xlab)) xlab = gene
  ylab = bquote(Predicted~ln(IC[50]))
  # ylab = bquote(bold(Predicted~ln(IC["50"])))
  
  corr = Pred %>% with(cor(Gene, Prediction))
  sprintf("Cell : %s", nrow(Pred)) %>% print
  sprintf("Corr [Gene] : %.3f", corr) %>% print
  
  if (trend_line) {
    add = add %>% append_def(geom_smooth(method="lm"))
  }
  
  if (anno_subtype) {
    add = add %>% append_def(labs(color="Subtype"))
    if (!is.null(SubType_List)) {
      add = add %>% append_def(guides(shape=NULL))
      Subtype = stack(SubType_List) %>% setNames(c("Cell", "Subtype"))
      Pred = full_join(Pred, Subtype, by="Cell", relationship="many-to-many")
      
      if (exclude_else) {
        Pred = Pred %>% subset(!is.na(Subtype))
      } else {
        levels = c(levels(Pred$Subtype), name_else)
        Pred$Subtype = Pred$Subtype %>% as.character
        Pred$Subtype[is.na(Pred$Subtype)] = name_else
        Pred$Subtype = Pred$Subtype %>% factor(levels=levels)
      }
    }
    
    Pred = Pred %>% 
      mutate(GDSC_Cell=!Rest) %>% 
      subset(select=-Rest) %>% 
      relocate(Subtype, .after=Cell) %>% 
      relocate(Drug_Name, .after=Drug) %>%
      rename(Drug_CID=Drug, IC50_Missing=Missing) %>% 
      relocate(GDSC_Cell, .after=IC50_Missing) %>% as.data.frame
  }
  
  # Annotation of cell names
  anno_cells = T
  if (anno_cells) {
    idx = match(Pred$Cell, Anno_Cells$SANGER_MODEL_ID)
    Pred = Pred %>% mutate(Cell_Name=Anno_Cells$MODEL_NAME[idx]) %>% 
      relocate(Cell_Name, .after=Cell) %>% as.data.frame
  }
  
  if (!is.null(marker)) {
    add = add %>% append_def(labs(size=marker))
    idx = match(Pred$Cell, rownames(TPM))
    Pred$Marker = TPM[idx, marker]
    Pred$Marker_Name = marker
    corr_marker = Pred %>% with(cor(Marker, Prediction))
    sprintf("Corr [Marker] : %.3f", corr_marker) %>% print
  }
  
  if (!is.null(cells)) {
    Pred_ = Pred %>% filter(!!sym(col_cell) %in% cells)
    point_cells = Pred_ %>% geom_point(mapping=aes(Gene, Prediction), color="red", shape=24, size=3)
    text_cells = Pred_ %>% text_repel_def(!!sym(col_cell), size=8)
    add = add %>% append_def(point_cells)
    add = add %>% append_def(text_cells)
  }
  
  if (length(add)==0) add = NULL
  if (anno_subtype & !is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, color=Subtype, shape=Subtype, size=Marker, 
                      width=width, height=height, axis_tl=axis_tl, axis_tx=axis_tx, 
                      legend_tl=legend_tl, legend_tx=legend_tx, ...)
  } else if (anno_subtype & is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, color=Subtype, shape=Subtype, size=2.5, 
                      width=width, height=height, axis_tl=axis_tl, axis_tx=axis_tx, 
                      legend_tl=legend_tl, legend_tx=legend_tx, ...)
  } else if (!anno_subtype & !is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, size=Marker, 
                      width=width, height=height, axis_tl=axis_tl, axis_tx=axis_tx, 
                      legend_tl=legend_tl, legend_tx=legend_tx, ...)
  } else {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, 
                      width=width, height=height, axis_tl=axis_tl, axis_tx=axis_tx, 
                      legend_tl=legend_tl, legend_tx=legend_tx, ...)
  }
  return(Pred)
}

sort_subtype = function(TPM, Marker_List, Weight_List=NULL,
                        cells=NULL, quant=1/2, name_else="The Others") {
  
  Cell_List = list()
  if (!is.null(cells)) TPM = TPM[rownames(TPM) %in% cells, ]
  cells_else = rownames(TPM)
  
  for (i in 1:length(Marker_List)) {
    weight = Weight_List[[i]]
    markers = Marker_List[[i]]
    
    if (length(markers)>1) {
      mean_w = function(x, ...) weighted.mean(x, w=weight, ...)
      mean_ = ifelse_def(is.null(weight), mean, mean_w)
      tpm = TPM[, markers] %>% scale %>% apply(1, mean_, na.rm=T)
    } else tpm = setNames(TPM[, markers], rownames(TPM))
    
    idx = tpm>=quantile(tpm, quant, na.rm=T)
    Cell_List[[i]] = idx[idx] %>% names
    cells_else = cells_else %>% setdiff(Cell_List[[i]])
  }
  
  names(Cell_List) = names(Marker_List)
  Cell_List[[name_else]] = Cell_List[[name_else]] %>% c(cells_else)
  return(Cell_List)
}

pca_subtype = function(TPM, SubType_List, main=main, thr_sd=0.01, name_else="The Others",
                       exclude_else=F, width=20, height=16, return_val=T, save=F, ...) {
  
  Subtype = stack(SubType_List) %>% setNames(c("Cell", "Subtype"))
  if (exclude_else) TPM = TPM[rownames(TPM) %in% Subtype$Cell, ]
  cond_sd = sapply(TPM, sd)>=thr_sd
  TPM_PCA = TPM[, cond_sd] %>% prcomp(center=T, scale=T, rank=2)
  TPM_PCA = TPM_PCA$x %>% as.data.frame
  
  TPM_PCA$Cell = rownames(TPM_PCA)
  TPM_PCA = full_join(TPM_PCA, Subtype, by="Cell", relationship="many-to-many")
  
  if (!exclude_else) {
    levels = c(levels(TPM_PCA$Subtype), name_else)
    TPM_PCA$Subtype = TPM_PCA$Subtype %>% as.character
    TPM_PCA$Subtype[is.na(TPM_PCA$Subtype)] = name_else
    TPM_PCA$Subtype = TPM_PCA$Subtype %>% factor(levels=levels)
  }
  
  TPM_PCA = TPM_PCA %>% relocate(PC1, PC2, .after=everything()) %>% as.data.frame
  TPM_PCA %>% plot_def(PC1, PC2, main=main, color=Subtype, shape=Subtype,
                       size=3, stroke=1.5, alpha=0.8, margin=0.4, 
                       axis_tl=30, axis_tx=25, legend_tl=22.5, legend_tx=20,
                       width=width, height=height, save=save, ...)
  
  if (return_val) return(TPM_PCA)
}

boxplot_subtype = function(Pred, SubType_List, name_else="The Others", main=NULL, 
                           width=20, height=18, size_psig=6, exclude_else=T, use_ggpubr=T, save=F, ...) {
  
  suppressMessages(library(ggpubr))
  Subtype = stack(SubType_List) %>% setNames(c("Cell", "Subtype"))
  Pred = full_join(Pred, Subtype, by="Cell", relationship="many-to-many")
  
  if (exclude_else) {
    Pred = Pred %>% subset(!is.na(Subtype))
  } else {
    levels = c(levels(Pred$Subtype), name_else) %>% unique
    Pred$Subtype = Pred$Subtype %>% as.character
    Pred$Subtype[is.na(Pred$Subtype)] = name_else
    Pred$Subtype = Pred$Subtype %>% factor(levels=levels)
  }
  
  Pred = Pred %>% 
    mutate(GDSC_Cell=!Rest) %>% 
    subset(select=-Rest) %>% 
    relocate(Subtype, .after=Cell) %>% 
    relocate(Drug_Name, .after=Drug) %>%
    rename(Drug_CID=Drug, IC50_Missing=Missing) %>% 
    relocate(GDSC_Cell, .after=IC50_Missing) %>% as.data.frame
  
  # Annotation of cell names
  anno_cells = T
  if (anno_cells) {
    idx = match(Pred$Cell, Anno_Cells$SANGER_MODEL_ID)
    Pred = Pred %>% mutate(Cell_Name=Anno_Cells$MODEL_NAME[idx]) %>% 
      relocate(Cell_Name, .after=Cell) %>% as.data.frame
  }
  
  if (!use_ggpubr) {
    Pred %>% subset(!is.na(Prediction)) %>%
      boxplot_def(Subtype, Prediction, legend=F, force_bold=F,
                  width=widt, height=height, save=save, ...)
  } else {
    section = function(x, quant=0.5, na.rm=T) {
      if (na.rm) x = x %>% na.omit %>% as.numeric
      min(x)+(max(x)-min(x))*quant
    }
    
    ylab = bquote(ln(IC[50]))
    pos = position_dodge(width=0.8)
    margin1 = margin(5, 5, 5, 5, unit="pt")
    margin2 = margin(10, 10, 10, 10, unit="pt")
    margin3 = margin(0, 25, 0, 25, unit="pt")
    margin4 = margin(0, 25, 0, 10, unit="pt")
    
    # add.params = list(alpha=0.5, size=2)
    # ylab = bquote(Predicted~ln(IC[50]))
    
    font1 = font(object="ylab", size=30, margin=margin1)
    font2 = font(object="axis.text", size=22.5, margin=margin2, color="grey30")
    font3 = font(object="legend.title", size=25, margin=margin3)
    font4 = font(object="legend.text", size=25, margin=margin4)
    font = font1 + font2 + font3 + font4
    
    # pl = Pred %>% subset(!is.na(Prediction)) %>%
    #   ggboxplot(x="Subtype", y="Prediction", fill="Subtype", outlier.shape=NA,
    #             add="point", size=1, alpha=0.9, xlab=F, add.params=add.params, ...) +
    #   labs(y=ylab) + font + rotate_x_text(angle=30, hjust=1, vjust=1)
    
    col = c("LN_IC50", "Prediction")
    col_af = c("Response_Type", "Response_Value")
    color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
    
    labels = c("Actual", "Predicted")
    # labels = c(bquote(Actual~ln(IC[50])), bquote(Predicted~ln(IC[50])))
    
    Pred_ = Pred %>% as.data.frame %>% 
      reshape2::melt(measure.vars=col, variable.name=col_af[1], value.name=col_af[2]) %>%
      subset(!is.na(Response_Value)) %>% 
      mutate(Response_Type=ifelse(Response_Type=="LN_IC50", labels[1], labels[2])) %>% 
      mutate(Response_Type=Response_Type %>% factor(levels=labels))
    
    pl = Pred_ %>% 
      ggboxplot(x="Subtype", y="Response_Value", fill="Response_Type",
                outlier.shape=NA, size=1, alpha=0.9, xlab=F, ...) +
      geom_point(aes(group=Response_Type), size=2, alpha=0.5, position=pos) + 
      labs(y=ylab, fill="Response") +
      rotate_x_text(angle=30, hjust=1, vjust=1) + 
      scale_fill_manual(labels=labels, values=color) + font
    
    # 1. Comparison within group [Predicted vs Actual]
    pl = pl + geom_pwc(aes(group=Response_Type), 
                       label="p.adj.signif", p.adjust.method="fdr", label.size=size_psig,
                       tip.length=0, bracket.nudge.y=0.05, step.increase=0.12, hide.ns=F)
    
    # 2. Comparison between group [SCLC Type & NSCLC]
    pl = pl + geom_pwc(data=Pred_ %>% subset(Response_Type!="Actual"), 
                       label="p.adj.signif", p.adjust.method="fdr", label.size=size_psig,
                       tip.length=0.02, bracket.nudge.y=0.15+0.05, step.increase=0.12, hide.ns=T)
    
    pl = pl %>% ggpar(legend="bottom") + theme(legend.key.size=unit(1, "cm"))
    # pl = pl %>% ggpar(legend="none")
    
    if (save) {
      save_fig(pl, main, width=width, height=height, units="cm", svg=T)
    } else print(pl) ; return(Pred)
  }
  return(Pred)
}

boxplot_subtype_v2 = function(Pred, SubType_List, col_tn=NULL, col_cell="Sample", 
                              main=NULL, lvl_subtype=NULL, name_else="The Others", 
                              width=20, height=18, size_psig=6, exclude_else=T, save=F, ...) {
  
  suppressMessages(library(ggpubr))
  Subtype = stack(SubType_List) %>% setNames(c(col_cell, "Subtype"))
  Pred = full_join(Pred, Subtype, by=col_cell, relationship="many-to-many")
  if (is.null(lvl_subtype)) lvl_subtype = Pred$Subtype %>% unique %>% sort
  
  if (exclude_else) {
    Pred = Pred %>% subset(!is.na(Subtype))
  } else {
    lvl_subtype = c(lvl_subtype, name_else) %>% unique
    Pred$Subtype = Pred$Subtype %>% as.character
    Pred$Subtype[is.na(Pred$Subtype)] = name_else
  }
  
  Pred$Subtype = Pred$Subtype %>% factor(levels=lvl_subtype)
  Pred = Pred %>% relocate(Subtype, .after=!!sym(col_cell)) %>% as.data.frame
  
  ylab = bquote(Predicted~ln(IC[50]))
  pos = position_dodge(width=0.8)
  margin1 = margin(5, 5, 5, 5, unit="pt")
  margin2 = margin(10, 10, 10, 10, unit="pt")
  margin3 = margin(0, 25, 0, 25, unit="pt")
  margin4 = margin(0, 25, 0, 10, unit="pt")
  
  font1 = font(object="ylab", size=30, margin=margin1)
  font2 = font(object="axis.text", size=22.5, margin=margin2, color="grey30")
  font3 = font(object="legend.title", size=25, margin=margin3)
  font4 = font(object="legend.text", size=25, margin=margin4)
  font = font1 + font2 + font3 + font4
  
  if (is.null(col_tn)) {
    fill = RColorBrewer::brewer.pal(8, "Reds")[8]
    fill_labels = NULL
  } else {
    fill = col_tn
    labels = c("Tumor", "Normal")
    fill_labels = "Status"
    color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 3)]
  }
  
  pl = Pred %>% 
    ggboxplot(x="Subtype", y="Prediction", fill=fill,
              outlier.shape=NA, size=1, alpha=0.9, xlab=F, ...) +
    labs(y=ylab, fill=fill_labels) +
    rotate_x_text(angle=30, hjust=1, vjust=1) + font
  
  if (!is.null(col_tn)) {
    # 1. Comparison within group [Tumor vs Normal]
    Pred_ = Pred %>% filter(!!sym(col_tn)=="Tumor")
    
    pl = pl + scale_fill_manual(labels=labels, values=color) +
      geom_point(aes(group=!!sym(col_tn)), size=2, alpha=0.5, position=pos)
    
    pl = pl + geom_pwc(aes(group=!!sym(col_tn)), 
                       label="p.adj.signif", p.adjust.method="fdr", label.size=size_psig,
                       tip.length=0, bracket.nudge.y=0.05, step.increase=0.12, hide.ns=F)
  } else {
    Pred_ = Pred
    pl = pl + geom_point(size=2, alpha=0.5, position=pos)
  }
  
  # 2. Comparison between group [Subtypes]
  pl = pl + geom_pwc(data=Pred_, label="p.adj.signif", p.adjust.method="fdr", label.size=size_psig,
                     tip.length=0.02, bracket.nudge.y=0.15+0.05, step.increase=0.12, hide.ns=T)
  
  pl = pl %>% ggpar(legend="bottom") + theme(legend.key.size=unit(1, "cm"))
  # pl = pl %>% ggpar(legend="none")
  
  if (save) {
    save_fig(pl, main, width=width, height=height, units="cm", svg=T)
  } else print(pl) ; return(Pred)
}

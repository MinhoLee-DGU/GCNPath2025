#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(cogena))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))

source("../utils/functions.R")
loadings()


dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC.txt", dir)
IC50_GDSC = read.csv(file, sep="\t")

dir = "../processed/cell_data_biocarta"
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
SANGER_RNA_GSVA = read.csv(file, row.names=1)

dir = "../data/path_data"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.symbols.gmt", dir)
Path_List = gmt2list(file)
geneset = unique(unlist(Path_List))

file = "TPM_Sym.csv"
SANGER_RNA_TPM = fread_def(file)
# SANGER_RNA_TPM = SANGER_RNA_TPM[, colnames(SANGER_RNA_TPM) %in% geneset]

file = "Anno_Cells.csv"
Anno_Cells = read.csv(file)
Anno_Cells$TCGA_CODE[is.na(Anno_Cells$TCGA_CODE)] = "UNCLASSIFIED"

file = "Anno_Drugs.csv"
Anno_Drugs = read.csv(file)
Anno_Drugs$Drug_CID = Anno_Drugs$Drug_CID %>% as.character
Anno_Drugs$Target_Pathway[is.na(Anno_Drugs$Target_Pathway)] = "Unclassified"


IC50_GDSC_Miss = IC50_GDSC %>% complete(Cell, Drug) %>% 
  subset(is.na(LN_IC50)) %>% as.data.frame

idx = match(IC50_GDSC_Miss$Cell, IC50_GDSC$Cell)

IC50_GDSC_Miss = IC50_GDSC_Miss %>% 
  subset(select=-LN_IC50) %>% 
  mutate(Cell_BROAD=IC50_GDSC$Cell_BROAD[idx], 
         Cell_COSMIC=IC50_GDSC$Cell_COSMIC[idx])

dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC_Miss.txt", dir)
write.table(IC50_GDSC_Miss, file=file, row.names=F, sep="\t")



##### 2-1. Prediction of missing IC50 values (Interporation)

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

pattern = "pred_total_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Test = read_pred(dir, pattern=pattern)   # 3717510
Pred_Test$Drug = Pred_Test$Drug %>% as.character

pattern = "pred_miss_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Miss = read_pred(dir, pattern=pattern)   # 481530
Pred_Miss$Drug = Pred_Miss$Drug %>% as.character
Pred_Miss = Pred_Miss %>% mutate(LN_IC50=NA) %>% relocate(LN_IC50, .before=Prediction)

# GCNPath shows outstanding and stable performances...
Perf_Test = Pred_Test %>% group_by(Seed) %>% 
  summarise(RMSE=RMSE(LN_IC50, Prediction), 
            PCC=cor(LN_IC50, Prediction), 
            SCC=cor(LN_IC50, Prediction, method="spearman")) %>% as.data.frame

Perf_Test$RMSE %>% mean %>% round(3)   # 0.691
Perf_Test$PCC %>% mean %>% round(3)    # 0.967
Perf_Test$SCC %>% mean %>% round(3)    # 0.952

col = c("Cell", "Drug", "LN_IC50", "Prediction")
Pred_Test = Pred_Test[, col, with=F] %>% as.data.frame %>% mutate(Missing=F)
Pred_Miss = Pred_Miss[, col, with=F] %>% as.data.frame %>% mutate(Missing=T)

Pred_Test = Pred_Test %>% group_by(Cell, Drug) %>% 
  mutate(Prediction=mean(Prediction)) %>% 
  subset(select=c(Cell, Drug, LN_IC50, Prediction, Missing)) %>% distinct %>% as.data.frame

Pred_Miss = Pred_Miss %>% group_by(Cell, Drug) %>% 
  summarise(LN_IC50=NA, Prediction=mean(Prediction), Missing=T) %>% as.data.frame

Pred_Test = rbind(Pred_Test, Pred_Miss)
Pred_Test$Drug = Pred_Test$Drug %>% as.character

Pred_Test %>% nrow                     # 419904 [972*432]
Pred_Test$Cell %>% unique %>% length   # 972
Pred_Test$Drug %>% unique %>% length   # 432

Pred_Test %>% subset(Missing) %>% nrow    # 48153
Pred_Test %>% subset(!Missing) %>% nrow   # 371751
all(Pred_Test$Missing[is.na(Pred_Test$LN_IC50)])     # T
all(!Pred_Test$Missing[!is.na(Pred_Test$LN_IC50)])   # T



##### 2-2. Prediction of cell-lines without IC50 values

# Run the code "test_rest.sh" to predict the Pred_Rest
cells = rownames(SANGER_RNA_GSVA) %>% setdiff(unique(Pred_Test$Cell))   # 459
drugs = unique(Pred_Test$Drug)   # 432

dir = "../data/ic50_data"
file = sprintf("%s/IC50_GDSC_Rest.txt", dir)
Combn_Rest = expand.grid(Cell=cells, Drug=drugs) %>% as.data.frame
write.table(Combn_Rest, file=file, sep="\t", quote=F, row.names=F)

pattern = "pred_rest_seed([0-9]+).csv"
dir = "../results/IC50_GDSC/Normal/RGCN"
Pred_Rest = read_pred(dir, pattern=pattern)   # 1982880 [198288*10]
Pred_Rest$Drug = Pred_Rest$Drug %>% as.character

Pred_Rest = Pred_Rest %>% group_by(Cell, Drug) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 198288

Pred_Test$Rest = F
Pred_Rest = Pred_Rest %>% 
  mutate(LN_IC50=NA, Missing=F, Rest=T) %>% 
  relocate(LN_IC50, .before=Prediction) %>% as.data.frame

Pred_Test = Pred_Test %>% rbind(Pred_Rest)
Pred_Wide = Pred_Test %>% reshape2::acast(Cell~Drug, value.var="Prediction") %>% as.data.frame   # 1431 x 432


# Counting TCGA Code & Target Pathway [Cell & Drug]
Anno_Cells_SANGER = Anno_Cells %>% 
  subset(SANGER_MODEL_ID %in% rownames(SANGER_RNA_GSVA)) %>% 
  group_by(CANCER_TYPE_DETAIL) %>% summarise(Num=n()) %>% 
  arrange(desc(Num)) %>% as.data.frame

Anno_Drugs_SANGER = Anno_Drugs %>% 
  subset(Drug_CID %in% unique(IC50_GDSC$Drug)) %>% 
  group_by(Target_Pathway) %>% summarise(Num=n()) %>% 
  arrange(desc(Num)) %>% as.data.frame

# Complete cell annotation
# cells included in IC50_GDSC & SANGER_RNA?

IC50_GDSC_Num = IC50_GDSC %>% group_by(Cell) %>% 
  summarise(Num_IC50=n()) %>% as.data.frame
idx = match(Anno_Cells$SANGER_MODEL_ID, IC50_GDSC_Num$Cell)

Anno_Cells = Anno_Cells %>% 
  mutate(RNA=SANGER_MODEL_ID %in% rownames(SANGER_RNA_GSVA), 
         IC50=SANGER_MODEL_ID %in% unique(IC50_GDSC$Cell), 
         Num_IC50=IC50_GDSC_Num$Num_IC50[idx])

Anno_Cells$MODEL_NAME %>% is.na %>% sum   # 0



##### Examine cells by subtypes [SCLC]

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
                            filt_ic50=F, filt_nc=T, num_min=3) {
  Anno_Cells %>%
    filter_ic50(filter=filt_ic50) %>% 
    filter_non_can(filter=filt_nc) %>% 
    subset(RNA & TISSUE %in% tissue) %>%
    group_by(across(all_of(group))) %>% summarise(Num=n()) %>%
    arrange(desc(Num)) %>% subset(Num>=num_min) %>% as.data.frame
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
                         exclude_else=T, name_else="The Others", add=NULL, ...) {
  
  if (is.null(add)) add = list()
  idx = match(Pred$Cell, rownames(TPM))
  Pred$Gene = TPM[idx, gene]
  
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
    add = add %>% append_def(labs(color="SubType"))
    if (!is.null(SubType_List)) {
      add = add %>% append_def(guides(shape=NULL))
      SubType = stack(SubType_List) %>% setNames(c("Cell", "SubType"))
      Pred = full_join(Pred, SubType, by="Cell", relationship="many-to-many")
      
      if (exclude_else) {
        Pred = Pred %>% subset(!is.na(SubType))
      } else {
        levels = c(levels(Pred$SubType), name_else)
        Pred$SubType = Pred$SubType %>% as.character
        Pred$SubType[is.na(Pred$SubType)] = name_else
        Pred$SubType = Pred$SubType %>% factor(levels=levels)
      }
    }
    
    Pred = Pred %>% 
      mutate(GDSC_Cell=!Rest) %>% 
      subset(select=-Rest) %>% 
      relocate(SubType, .after=Cell) %>% 
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
    corr_marker = Pred %>% with(cor(Marker, Prediction))
    sprintf("Corr [Marker] : %.3f", corr_marker) %>% print
  }
  
  if (length(add)==0) add = NULL
  if (anno_subtype & !is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, color=SubType, shape=SubType, size=Marker, ...)
  } else if (anno_subtype & is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, color=SubType, shape=SubType, size=2.5, ...)
  } else if (!anno_subtype & !is.null(marker)) {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, size=Marker, ...)
  } else {
    Pred %>% plot_def(Gene, Prediction, xlab=xlab, ylab=ylab, alpha=0.8, 
                      margin=0.5, add=add, ...)
  }
  return(Pred)
}

# Lung
tissue = "Lung"
Anno_Cells %>% lookup_by_tissue(tissue)
Anno_Cells %>% lookup_by_tissue(tissue, group="CANCER_TYPE_DETAIL")
cancer_type = Anno_Cells %>% lookup_by_tissue(tissue) %>% pull(CANCER_TYPE)

except_type = c("Mesothelioma", "Other Solid Cancers")
cancer_type = cancer_type %>% setdiff(except_type)

Pred_Lung = examine_subtype(Pred_Test, Anno_Cells, Anno_Drugs, cancer_type)
Pred_Lung$Pred_Drug %>% head
Pred_Lung$Pred_Pathway %>% head


# SCLC
cancer = "Small Cell Lung Carcinoma"
cells_sclc = Anno_Cells %>% subset(CANCER_TYPE_DETAIL==cancer) %>% pull(SANGER_MODEL_ID)

# SCLC is classified as SCLC-A/N/P/I
# doi.org/10.3322/caac.21785
# [Review, 2023] Clinical insights into small cell lung cancer
# CD56 [NCAM1], LSD1 [KDM1A], PDL1 [CD274], PD1 [PDCD1], VISTA [VSIR]
# LSH1 gene in SCLC-N markers was not found...


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

marker_a = "ASCL1"
marker_n = "NEUROD1"
marker_p = "POU2F3"
markers = c(marker_a, marker_n, marker_p, "PD-1", "PD-L1")
markers_ori = c(marker_a, marker_n, marker_p, "CD274", "PDCD1")

dir = mkdir("Case Study [SCLC]")
main = sprintf("%s/TPM Histogram [%s]", dir, markers)
median_tpm = SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, markers_ori] %>% sapply(median)

xlab = sprintf("Expression of %s", markers)
# xlab = sprintf('bquote(bold(log["2"]~("TPM+1")~"[%s]"))', markers)
# xlab = xlab %>% lapply(function(x) eval(parse(text=x)))

# median_tpm %>% round(2)
# ASCL1   NEUROD1   POU2F3   CD274   PDCD1 
# 9.69    2.30      0.11     1.00    0.00 

SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, marker_a] %>% 
  hist_def(main=main[1], xlab=xlab[1], margin=0.4, margin_pl=0.8,
           axis_tl=27, axis_tx=22.5, dist=F, text_info=F, save=T)
SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, marker_n] %>% 
  hist_def(main=main[2], xlab=xlab[2], margin=0.4, margin_pl=0.8,
           axis_tl=27, axis_tx=22.5, dist=F, text_info=F, save=T)
SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, marker_p] %>% 
  hist_def(main=main[3], xlab=xlab[3], margin=0.4, margin_pl=0.8,
           axis_tl=27, axis_tx=22.5, dist=F, text_info=F, save=T)
SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, "CD274"] %>% 
  hist_def(main=main[4], xlab=xlab[4], margin=0.4, margin_pl=0.8,
           axis_tl=27, axis_tx=22.5, dist=F, text_info=F, save=T)   # PD1 [SCLC-I]
SANGER_RNA_TPM[rownames(SANGER_RNA_TPM) %in% cells_sclc, "PDCD1"] %>% 
  hist_def(main=main[5], xlab=xlab[5], margin=0.4, margin_pl=0.8,
           axis_tl=27, axis_tx=22.5, dist=F, text_info=F, save=T)   # PDL1 [SCLC-I]

labels = sprintf("SCLC-%s", c("A", "N", "P"))
Marker_List = list(marker_a, marker_n, marker_p)
names(Marker_List) = labels

quant = 1/2
SCLC_Type = SANGER_RNA_TPM %>%
  sort_subtype(Marker_List, cells=cells_sclc, name_else="SCLC-I", quant=quant)

main = sprintf("%s/Venn diagram of SCLC subtypes", dir)
SCLC_Type %>% venn_def(main=main, set_size=7.5, text_size=8.4, 
                       width=24, height=18, legend=F, save=T)

main = sprintf("%s/Venn diagram of SCLC subtypes [No Labels]", dir)
SCLC_Type %>% venn_def(main=main, labels=F, text_size=8.4, 
                       width=24, height=18, legend=F, save=T)


pca_subtype = function(TPM, SubType_List, main=main, thr_sd=0.01, name_else="The Others",
                       exclude_else=F, width=20, height=16, return_val=T, save=F, ...) {
  
  SubType = stack(SubType_List) %>% setNames(c("Cell", "SubType"))
  if (exclude_else) TPM = TPM[rownames(TPM) %in% SubType$Cell, ]
  cond_sd = sapply(TPM, sd)>=thr_sd
  TPM_PCA = TPM[, cond_sd] %>% prcomp(center=T, scale=T, rank=2)
  TPM_PCA = TPM_PCA$x %>% as.data.frame
  
  TPM_PCA$Cell = rownames(TPM_PCA)
  TPM_PCA = full_join(TPM_PCA, SubType, by="Cell", relationship="many-to-many")
  
  if (!exclude_else) {
    levels = c(levels(TPM_PCA$SubType), name_else)
    TPM_PCA$SubType = TPM_PCA$SubType %>% as.character
    TPM_PCA$SubType[is.na(TPM_PCA$SubType)] = name_else
    TPM_PCA$SubType = TPM_PCA$SubType %>% factor(levels=levels)
  }
  
  TPM_PCA = TPM_PCA %>% relocate(PC1, PC2, .after=everything()) %>% as.data.frame
  TPM_PCA %>% plot_def(PC1, PC2, main=main, color=SubType, shape=SubType,
                       size=3, stroke=1.5, alpha=0.8, margin=0.4, 
                       axis_tl=30, axis_tx=25, legend_tl=22.5, legend_tx=20,
                       width=width, height=height, save=save, ...)
  
  if (return_val) return(TPM_PCA)
}

boxplot_subtype = function(Pred, SubType_List, name_else="The Others",
                           main=NULL, pos_anova=0.95, width=20, height=18,
                           size_psig=6, exclude_else=T, use_ggpubr=T, save=F, ...) {

  suppressMessages(library(ggpubr))
  SubType = stack(SubType_List) %>% setNames(c("Cell", "SubType"))
  Pred = full_join(Pred, SubType, by="Cell", relationship="many-to-many")

  if (exclude_else) {
    Pred = Pred %>% subset(!is.na(SubType))
  } else {
    levels = c(levels(Pred$SubType), name_else)
    Pred$SubType = Pred$SubType %>% as.character
    Pred$SubType[is.na(Pred$SubType)] = name_else
    Pred$SubType = Pred$SubType %>% factor(levels=levels)
  }

  Pred = Pred %>% 
    mutate(GDSC_Cell=!Rest) %>% 
    subset(select=-Rest) %>% 
    relocate(SubType, .after=Cell) %>% 
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
      boxplot_def(SubType, Prediction, legend=F, force_bold=F,
                  width=widt, height=height, save=save, ...)
  } else {
    section = function(x, quant=0.5, na.rm=T) {
      if (na.rm) x = x %>% na.omit %>% as.numeric
      min(x)+(max(x)-min(x))*quant
    }

    # ylab = bquote(Predicted~ln(IC[50]))
    # pos_anova = unlist(Pred$Prediction) %>% section(quant=pos_anova, na.rm=T)
    # stat_compare_means(method="anova", label.y=pos_anova)
    
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
    #   ggboxplot(x="SubType", y="Prediction", fill="SubType", outlier.shape=NA,
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
      ggboxplot(x="SubType", y="Response_Value", fill="Response_Type",
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

# Paired wilcox.test
# When implementing paired t-test, it is not necessary to specify the sample ID column
# Just arrange sample IDs from two groups in the same order...
# https://rpkgs.datanovia.com/rstatix/reference/wilcox_test.html
# https://www.datanovia.com/en/lessons/how-to-do-a-t-test-in-r-calculation-and-reporting/#paired-t-test

# Some cell-lines do not have actual ln(IC50) values
# This causes the error of paired wilcox.test
# So we decided to implement unpaired test within each subtype

pc2_type = c("LN_IC50", "Log2TPM", "GSVA")
shape = c(22, 23, 24, 25, 21)
color = c("brown1", "gold", "seagreen3", "royalblue1", "mediumorchid")

add = list(scale_color_manual(values=color), 
           scale_shape_manual(values=shape), 
           theme(legend.key.size=unit(1, "cm")))

dir = "Case Study [SCLC]"
main = sprintf("%s/PC2 Plot [%s]", dir, pc2_type)
cells_lung = rownames(SANGER_RNA_TPM) %in% unique(Pred_Lung$Pred$Cell)   # 202

PC2_IC50 = Pred_Wide[cells_lung, ] %>% 
  pca_subtype(SCLC_Type, main=main[1], name_else="NSCLC", add=add, save=T)
PC2_TPM = SANGER_RNA_TPM[cells_lung, ] %>% 
  pca_subtype(SCLC_Type, main=main[2], name_else="NSCLC", add=add, save=T)
PC2_GSVA = SANGER_RNA_GSVA[cells_lung, ] %>% 
  pca_subtype(SCLC_Type, main=main[3], name_else="NSCLC", add=add, save=T)

drug_bcl2 = c("Venetoclax", "Navitoclax", "Obatoclax Mesylate", "Sabutoclax", "ABT737", "TW 37")
drug_aurka = c("Alisertib", "Tozasertib", "ZM447439", "CD532", "HG-5-113-01", "GSK1070916")
drug_igf1r = c("Linsitinib", "BMS-536924", "BMS-754807", "GSK1904529A", "NVP-ADW742")
drug_parp = c("Niraparib", "Olaparib", "Rucaparib", "Talazoparib", "Veliparib")


targets = c("BCL2", "AURKA", "IGF1R", "PARP")
dir = sprintf("Case Study [SCLC]/Target %s", targets)
for (d in dir) mkdir(d)

Pred_Lung_BCL2 = list()
for (drug in drug_bcl2) {
  main = sprintf("%s/SCLC Prediction [Target %s, %s]", dir[1], targets[1], drug)
  Pred_Lung_BCL2[[drug]] = Pred_Lung$Pred %>% subset(Drug_Name==drug) %>% 
    boxplot_subtype(SCLC_Type, main=main, exclude_else=F, name_else="NSCLC", save=T)
}

Pred_Lung_AURKA = list()
for (drug in drug_aurka) {
  main = sprintf("%s/SCLC Prediction [Target %s, %s]", dir[2], targets[2], drug)
  Pred_Lung_AURKA[[drug]] = Pred_Lung$Pred %>% subset(Drug_Name==drug) %>% 
    boxplot_subtype(SCLC_Type, main=main, exclude_else=F, name_else="NSCLC", save=T)
}

Pred_Lung_IGF1R = list()
for (drug in drug_igf1r) {
  main = sprintf("%s/SCLC Prediction [Target %s, %s]", dir[3], targets[3], drug)
  Pred_Lung_IGF1R[[drug]] = Pred_Lung$Pred %>% subset(Drug_Name==drug) %>% 
    boxplot_subtype(SCLC_Type, main=main, exclude_else=F, name_else="NSCLC", save=T)
}

Pred_Lung_PARP = list()
for (drug in drug_parp) {
  main = sprintf("%s/SCLC Prediction [Target %s, %s]", dir[4], targets[4], drug)
  Pred_Lung_PARP[[drug]] = Pred_Lung$Pred %>% subset(Drug_Name==drug) %>% 
    boxplot_subtype(SCLC_Type, main=main, exclude_else=F, name_else="NSCLC", save=T)
}

# Bcl-2 [SCLC-A]
# Venetoclax, Navitoclax, Obatoclax, Sabutoclax
dir_ex = mkdir(sprintf("%s [Example]", dir[1]))
main1 = sprintf("%s/Venetoclax-BCL2-ASCL1", dir_ex)
main2 = sprintf("%s/Navitoclax-BCL2-ASCL1", dir_ex)
main3 = sprintf("%s/Venetoclax-BCL2-Apoptosis", dir_ex)
main4 = sprintf("%s/Navitoclax-BCL2-Apoptosis", dir_ex)

shape = c(22, 23, 24, 25, 21)
color = c("brown1", "gold", "seagreen3", "royalblue1", "mediumorchid")

add = list(scale_color_manual(values=color), 
           scale_shape_manual(values=shape), 
           theme(legend.key.size=unit(0.9, "cm")))

xlab = "Expression of BCL2"

Pred_Ve_BCL2 = Pred_Lung$Pred %>% 
  subset(Drug_Name=="Venetoclax") %>% 
  plot_with_tpm(SANGER_RNA_TPM, gene="BCL2", marker="ASCL1", xlab=xlab, add=add,
                SubType_List=SCLC_Type, anno_subtype=T, exclude_else=F,
                main=main1, name_else="NSCLC", width=20, height=16, 
                axis_tl=30, axis_tx=22.5, legend_tl=20, legend_tx=18, save=T)
# [BCL2] -0.759, [ASCL1] -0.603

Pred_Na_BCL2 = Pred_Lung$Pred %>% 
  subset(Drug_Name=="Navitoclax") %>% 
  plot_with_tpm(SANGER_RNA_TPM, gene="BCL2", marker="ASCL1", xlab=xlab, add=add,
                SubType_List=SCLC_Type, anno_subtype=T, exclude_else=F,
                main=main2, name_else="NSCLC", width=20, height=16, 
                axis_tl=30, axis_tx=22.5, legend_tl=20, legend_tx=18, save=T)
# [BCL2] -0.735, [ASCL1] -0.665

path = "BIOCARTA_DEATH_PATHWAY"
xlab = "BIOCARTA Apoptosis"

Pred_Ve_Apop = Pred_Lung$Pred %>% 
  subset(Drug_Name=="Venetoclax") %>% 
  plot_with_tpm(SANGER_RNA_GSVA, gene=path, xlab=xlab, add=add, 
                SubType_List=SCLC_Type, anno_subtype=T, exclude_else=F,
                main=main3, name_else="NSCLC", width=20, height=15, 
                axis_tl=30, axis_tx=22.5, legend_tl=20, legend_tx=20, save=T)
# [Apoptosis] 0.408

Pred_Na_Apop = Pred_Lung$Pred %>% 
  subset(Drug_Name=="Navitoclax") %>% 
  plot_with_tpm(SANGER_RNA_GSVA, gene=path, xlab=xlab, add=add, 
                SubType_List=SCLC_Type, anno_subtype=T, exclude_else=F,
                main=main4, name_else="NSCLC", width=20, height=15, 
                axis_tl=30, axis_tx=22.5, legend_tl=20, legend_tx=20, save=T)
# [Apoptosis] 0.429


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
  
  
  ### [Source Data] Fig. 6
  
  Pred_Ve = Pred_Ve_BCL2 %>% 
    rename(BCL2=Gene, ASCL1=Marker, 
           Drug_Target_Pathway=Drug_Pathway) %>% 
    mutate(BIOCARTA_Apoptosis=Pred_Ve_Apop$Gene)
  
  Pred_Na = Pred_Na_BCL2 %>% 
    rename(BCL2=Gene, ASCL1=Marker, 
           Drug_Target_Pathway=Drug_Pathway) %>% 
    mutate(BIOCARTA_Apoptosis=Pred_Na_Apop$Gene)
  
  Pred_Ve_Na_ = list(Pred_Ve, Pred_Na)
  
  num_fig = c("a-c", "d-f")
  Pred_Ve_Na_ %>% save_for_nc(num=6, suppl=F, num_fig=num_fig)
  
  
  ### [Source Data] Supplementary Fig. 34
  SCLC_Type_ = SCLC_Type %>% stack %>% setNames(c("Cell", "SubType"))
  cells_lung = rownames(SANGER_RNA_TPM) %in% unique(Pred_Lung$Pred$Cell)
  cells_lung = rownames(SANGER_RNA_TPM)[cells_lung]   # 202
  
  cells_nsclc = cells_lung[!(cells_lung %in% SCLC_Type_$Cell)]
  NSCLC_Type_ = data.frame(Cell=cells_nsclc, SubType="NSCLC")
  Lung_Type = rbind(SCLC_Type_, NSCLC_Type_)
  
  idx = match(Lung_Type$Cell, Anno_Cells$SANGER_MODEL_ID)
  Lung_Type$Cell_Name = Anno_Cells$MODEL_NAME[idx]
  
  markers_sclc = c(marker_a, marker_n, marker_p, "CD274", "PDCD1")
  markers_sclc_ = c(marker_a, marker_n, marker_p, "CD274 [PD-L1]", "PDCD1 [PD-1]")
  TPM_Marker = SANGER_RNA_TPM[cells_lung, markers_sclc]
  colnames(TPM_Marker) = markers_sclc_
  
  # 40 [=33+7] SCLC cell-lines have multiple subtype membership...
  SCLC_Type_$Cell %>% table %>% table   # 1 [32], 2 [33], 3 [7]
  SCLC_Feature_ = list(Lung_Type, PC2_IC50, PC2_TPM, PC2_GSVA, TPM_Marker)
  
  num_fig = c(letters[1:4], "e-i")
  SCLC_Feature_ %>% save_for_nc(num=34, suppl=T, num_fig=num_fig, rowNames=c(F, F, F, F, T))
  
  
  ### [Source Data] Supplementary Fig. 35
  Pred_Lung_BCL2 %>% save_for_nc(num=35, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 36
  Pred_Lung_AURKA %>% save_for_nc(num=36, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 37
  Pred_Lung_IGF1R %>% save_for_nc(num=37, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 38
  Pred_Lung_PARP %>% save_for_nc(num=38, suppl=T)
}

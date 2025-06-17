#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(cogena))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))

source("../functions.R")
source("functions.R")
loadings()


dir = "../../processed_data/cell_data/BIOCARTA"
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
SANGER_RNA_GSVA = read.csv(file, row.names=1)

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Sym.csv", dir)
SANGER_RNA_TPM = fread_def(file)

# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/Anno_Cells.csv", dir)
# Anno_Cells = read.csv(file)
# Anno_Cells$TCGA_CODE[is.na(Anno_Cells$TCGA_CODE)] = "UNCLASSIFIED"

# dir = "../../processed_data/drug_data/GDSC"
# file = sprintf("%s/Anno_Drugs.csv", dir)
# Anno_Drugs = read.csv(file)
# Anno_Drugs$Drug_CID = Anno_Drugs$Drug_CID %>% as.character
# Anno_Drugs$Target_Pathway[is.na(Anno_Drugs$Target_Pathway)] = "Unclassified"

dir = "../../project/_Liu_Qian_2024/"
file = sprintf("%s/1-s2.0-S0092867423013351-mmc6.xlsx", dir)
SubType = read.xlsx(file, sheet="Table S6A")

colnames(SubType) = c("Tumor", "Subtype")
SubType = SubType %>% 
  mutate(Normal=gsub("T", "N", Tumor)) %>% 
  relocate(Normal, .after=Tumor) %>% as.data.frame

SubType = SubType %>% reshape2::melt(id.vars=c("Subtype"))
colnames(SubType) = c("Subtype", "Status", "Sample")
SubType$Subtype = gsub("nmf", "NMF", SubType$Subtype)




##### 2-1. Prediction [ComBat X]

pattern = "pred_liu24_invivo_seed([0-9]+).csv"

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGSA_SG = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_TGSA_SG$Model = "TGDRP_SG"

dir = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_DRPreter_SG$Model = "DRPreter_SG"

dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_GCNPath$Model = "GCNPath"


Pred_List = list(Pred_TGSA_SG, Pred_DRPreter_SG, Pred_GCNPath)
Pred_List = Reduce(rbind, Pred_List)   # 38520 [214*6*3*10]

Pred_List = Pred_List %>% group_by(Sample, Drug, Model) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 3852 [214*6*3]

idx = match(Pred_List$Sample, SubType$Sample)
Pred_List$Status = SubType$Status[idx]

# drugs = Pred_List$Drug %>% unique
# models = Pred_List$Model %>% unique
# dir = mkdir("Boxplot per Subtypes")
# 
# for (model in models) {
#   for (drug in drugs) {
#     main = sprintf("%s/Boxplot per SCLC Subtypes [%s, %s]", dir, model, drug)
#     Pred_List %>% subset(Model==model & Drug==drug) %>% 
#       boxplot_subtype(SubType, col_tn="Status", main=main, save=T)
#   }
# }


Pred_List_TN = Pred_List %>% 
  mutate(Group=gsub("[T, N]", "", Sample)) %>% 
  reshape2::dcast(Drug+Model+Group~Status, value.var="Prediction") %>% 
  mutate(Diff_Pred=Tumor-Normal) %>% as.data.frame

Pred_List_TN = Pred_List_TN %>% 
  mutate(Sample=paste0("T", Group)) %>% 
  left_join(SubType, by="Sample")

labels = c("TGDRP_SG", "DRPreter_SG", "GCNPath")
Pred_List_TN = Pred_List_TN %>% 
  mutate(Model = Model %>% factor(levels=labels)) %>% 
  relocate(Sample, .before=everything()) %>% 
  relocate(Subtype, .after=Group) %>% 
  subset(select=-c(Group, Status)) %>% 
  rename(Pred_Tumor=Tumor, Pred_Normal=Normal)

drugs = Pred_List$Drug %>% unique
color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
hline = geom_hline(yintercept=0, color="red", linetype=2)

margin1 = margin(0, 0.25, 0, 0.25, unit="cm")
margin2 = margin(0, 0.125, 0, 0.125, unit="cm")

theme_lg = theme(legend.position="bottom",
                 legend.key.size=unit(0.8, 'cm'),
                 legend.title=element_text(margin=margin1), 
                 legend.text=element_text(margin=margin2))

add = list(hline, theme_lg, scale_fill_manual(labels=labels, values=color))

ylab = "△Pred [Tumor-Normal]"
dir = mkdir("Boxplot per Subtypes [Tumor-Normal, ComBat X]")

for (drug in drugs) {
  main = sprintf("%s/Boxplot per SCLC Subtypes [%s]", dir, drug)
  Pred_List_TN %>% subset(Drug==drug) %>% 
    boxplot_def(Subtype, Diff_Pred, fill=Model, main=main, 
                add=add, ylab=ylab, alpha=0.9, width=17.1, 
                axis_tl=18, axis_tx=18, legend_tl=15, legend_tx=15, save=T)
}



##### 2-2. Prediction [ComBat O]

pattern = "pred_liu24_invivo_seed([0-9]+)_combat.csv"

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGSA_SG_CB = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_TGSA_SG_CB$Model = "TGDRP_SG"

dir = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG_CB = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_DRPreter_SG_CB$Model = "DRPreter_SG"

dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath_CB = read_pred(dir, pattern=pattern)   # 12840 [214*6*10]
Pred_GCNPath_CB$Model = "GCNPath"


Pred_List_CB = list(Pred_TGSA_SG_CB, Pred_DRPreter_SG_CB, Pred_GCNPath_CB)
Pred_List_CB = Reduce(rbind, Pred_List_CB)   # 38520 [214*6*3*10]

Pred_List_CB = Pred_List_CB %>% group_by(Sample, Drug, Model) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 3852 [214*6*3]

idx = match(Pred_List_CB$Sample, SubType$Sample)
Pred_List_CB$Status = SubType$Status[idx]


Pred_List_TN_CB = Pred_List_CB %>% 
  mutate(Group=gsub("[T, N]", "", Sample)) %>% 
  reshape2::dcast(Drug+Model+Group~Status, value.var="Prediction") %>% 
  mutate(Diff_Pred=Tumor-Normal) %>% as.data.frame

Pred_List_TN_CB = Pred_List_TN_CB %>% 
  mutate(Sample=paste0("T", Group)) %>% 
  left_join(SubType, by="Sample")

labels = c("TGDRP_SG", "DRPreter_SG", "GCNPath")
Pred_List_TN_CB = Pred_List_TN_CB %>% 
  mutate(Model = Model %>% factor(levels=labels)) %>% 
  relocate(Sample, .before=everything()) %>% 
  relocate(Subtype, .after=Group) %>% 
  subset(select=-c(Group, Status)) %>% 
  rename(Pred_Tumor=Tumor, Pred_Normal=Normal)

margin1 = margin(0, 0.25, 0, 0.25, unit="cm")
margin2 = margin(0, 0.125, 0, 0.125, unit="cm")

theme_lg = theme(legend.position="bottom",
                 legend.key.size=unit(0.8, 'cm'),
                 legend.title=element_text(margin=margin1), 
                 legend.text=element_text(margin=margin2))

add = list(hline, theme_lg, scale_fill_manual(labels=labels, values=color))

ylab = "△Pred [Tumor-Normal]"
dir = mkdir("Boxplot per Subtypes [Tumor-Normal, ComBat O]")

for (drug in drugs) {
  main = sprintf("%s/Boxplot per SCLC Subtypes [%s]", dir, drug)
  Pred_List_TN_CB %>% subset(Drug==drug) %>% 
    boxplot_def(Subtype, Diff_Pred, fill=Model, main=main, 
                add=add, ylab=ylab, alpha=0.9, width=17.1, 
                axis_tl=18, axis_tx=18, legend_tl=15, legend_tx=15, save=T)
}


# subtypes = SubType$Subtype %>% unique %>% sort
# dir = mkdir("Boxplot per Subtypes [Tumor-Normal, ComBat O (by SubType)]")
# 
# for (subtype in subtypes) {
#   main = sprintf("%s/Boxplot per SCLC Subtypes [%s]", dir, subtype)
#   Pred_List_TN_CB %>% subset(Subtype==subtype) %>% 
#     boxplot_def(Drug, Diff_Pred, fill=Model, main=main, 
#                 add=add, ylab=ylab, alpha=0.9, width=20, 
#                 axis_tl=20, axis_tx=18, legend_tl=18, legend_tx=15, save=T)
# }




##### 2-3. Prediction [Proteome, ComBat]

pattern_prot = "pred_liu24_invivo_prot_seed([0-9]+).csv"
dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath_Prot = read_pred(dir, pattern=pattern_prot)   # 12840 [214*6*10]

Pred_GCNPath = Pred_GCNPath %>% group_by(Sample, Drug, Model) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 1284 [214*6]

Pred_GCNPath_Prot = Pred_GCNPath_Prot %>% group_by(Sample, Drug, Model) %>% 
  summarise(Prediction=mean(Prediction)) %>% as.data.frame   # 1284 [214*6]

Pred_GCNPath$Model = "GCNPath_TPM"
Pred_GCNPath_Prot$Model = "GCNPath_Proteome"

file = "Table S1H.csv"
Corr_TPM_Pt = read.csv(file)
all(Corr_TPM_Pt$adjusted.P.value<0.05)   # T

col = c("Sample", "Corr")
Corr_TPM_Pt = Corr_TPM_Pt[, 1:2] %>% setNames(col)

# Prediction [TPM & Proteome]
Pred_GCNPath_ = rbind(Pred_GCNPath, Pred_GCNPath_Prot) %>% 
  reshape2::dcast(Sample+Drug~Model, value.var="Prediction")

idx = match(Pred_GCNPath_$Sample, Corr_TPM_Pt$Sample)
Pred_GCNPath_$Corr_TPM_Prot = Corr_TPM_Pt$Corr[idx]

Pred_GCNPath_ = Pred_GCNPath_ %>% 
  dplyr::rename(Pred_TPM=GCNPath_TPM, Pred_Proteome=GCNPath_Proteome)

xlab = "Prediction with TPM"
ylab = "Prediction with Proteome"
legend = bquote(Corr["(TPM, Proteome)"])

add = list(scale_color_gradient(low="yellow", high="firebrick1"))
main = "Prediction [TPM & Proteome]"

Pred_GCNPath_ %>% 
  plot_def(Pred_TPM, Pred_Proteome, color=Corr_TPM_Prot, 
           shape=Drug, main=main, add=add, alpha=0.5, 
           xlab=xlab, ylab=ylab, legend=legend, xy_line=T, 
           axis_tl=24, axis_tx=24, width=20, save=T)

Pred_GCNPath_ %>% with(cor(Pred_Proteome, Pred_TPM))   # 0.8833926

dir = "../GCNPath/processed/cell_data_biocarta"
file1 = sprintf("%s/Liu24_invivo_RNA_GSVA.csv", dir)
file2 = sprintf("%s/Liu24_invivo_Prot_GSVA.csv", dir)

Liu24_RNA_GSVA = read.csv(file1, row.names=1)
Liu24_Prot_GSVA = read.csv(file2, row.names=1)

row = intersect(rownames(Liu24_RNA_GSVA), rownames(Liu24_Prot_GSVA))
col = intersect(colnames(Liu24_RNA_GSVA), colnames(Liu24_Prot_GSVA))
Liu24_RNA_GSVA = Liu24_RNA_GSVA[row, col]     # 214 x 292
Liu24_Prot_GSVA = Liu24_Prot_GSVA[row, col]   # 214 x 292

dir = "../../project/_Liu_Qian_2024"
file1 = sprintf("%s/Table S1D_RNA.csv", dir)
file2 = sprintf("%s/Table S1E_Prot.csv", dir)

Liu24_RNA = read.csv(file1, row.names=1)
Liu24_Prot = read.csv(file2, row.names=1)

Liu24_RNA = Liu24_RNA %>% t %>% as.data.frame
Liu24_Prot = Liu24_Prot %>% t %>% as.data.frame

row = intersect(rownames(Liu24_RNA), rownames(Liu24_Prot))
col = intersect(colnames(Liu24_RNA), colnames(Liu24_Prot))
Liu24_RNA = Liu24_RNA[row, col]     # 214 x 9428
Liu24_Prot = Liu24_Prot[row, col]   # 214 x 9428

col_raw = c()
col_gsva = c()
for (i in rownames(Liu24_RNA)) {
  corr_raw_ = cor.test(as.numeric(Liu24_RNA[i, ]), as.numeric(Liu24_Prot[i, ]), method="spearman")
  corr_gsva_ = cor.test(as.numeric(Liu24_RNA_GSVA[i, ]), as.numeric(Liu24_Prot_GSVA[i, ]), method="spearman")
  col_raw = col_raw %>% c(corr_raw_$estimate)
  col_gsva = col_gsva %>% c(corr_gsva_$estimate)
}

names(col_raw) = rownames(Liu24_RNA)
names(col_gsva) = rownames(Liu24_RNA)

col_raw %>% is.na %>% sum    # 0
col_gsva %>% is.na %>% sum   # 0
Corr_TPM_Pt = data.frame(PCC_Raw_Data=col_raw, PCC_GSVA=col_gsva) %>% 
  mutate(Sample=rownames(Liu24_RNA)) %>% 
  mutate(Status=ifelse(grepl("T", Sample), "Tumor", "Normal")) %>% 
  relocate(Sample, Status, .before=everything())

rownames(Corr_TPM_Pt) = NULL
Corr_TPM_Pt$PCC_Raw_Data %>% range %>% round(3)   # -0.121  0.580
Corr_TPM_Pt$PCC_GSVA %>% range %>% round(3)       # 0.000 0.825

ceil_def = function(x, digit=0) {
  ceiling(x * 10**digit) / 10**digit
}

floor_def = function(x, digit=0) {
  floor(x * 10**digit) / 10**digit
}

range_raw = Corr_TPM_Pt$PCC_Raw_Data %>% range
range_gsva = Corr_TPM_Pt$PCC_GSVA %>% range
min_ = min(range_raw[1], range_gsva[1]) %>% floor_def(2)
max_ = max(range_raw[2], range_gsva[2]) %>% ceil_def(2)
lim = c(min_, max_)

xlab = bquote(SCC["(TPM, Proteome)"]~"[Raw Data]")
ylab = bquote(SCC["(TPM, Proteome)"]~"[GSVA]")

main = "SCC in TPM & Proteome"
Corr_TPM_Pt %>% plot_def(PCC_Raw_Data, PCC_GSVA, color=Status, shape=Status, 
                         main=main, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, 
                         size=2, alpha=0.8, width=15, height=11.4, xy_line=T, save=T)

sample_raw = Corr_TPM_Pt %>% slice_max(PCC_Raw_Data) %>% pull(Sample)   # T170393 [0.579656]
sample_gsva = Corr_TPM_Pt %>% slice_max(PCC_GSVA) %>% pull(Sample)      # T161286 [0.8247455]

Ex_Raw = data.frame(
  Gene=colnames(Liu24_RNA), 
  TPM=Liu24_RNA[sample_raw, ] %>% as.numeric, 
  Proteome=Liu24_Prot[sample_raw, ] %>% as.numeric
)

Ex_GSVA = data.frame(
  Pathway=colnames(Liu24_RNA_GSVA),
  TPM=Liu24_RNA_GSVA[sample_gsva, ] %>% as.numeric, 
  Proteome=Liu24_Prot_GSVA[sample_gsva, ] %>% as.numeric
)

xlab = bquote(Log[2](TPM+1))
ylab = "Proteome"

main = sprintf("SCC in TPM & Proteome [Raw, %s]", sample_raw)
# Ex_Raw %>% na.omit %>% dim   # 8927
Ex_Raw %>% plot_def(TPM, Proteome, main=main, 
                    xlab=xlab, ylab=ylab, alpha=0.5, width=12, height=12, save=T)

xlab = bquote(GSVA~Score~"["~log[2]~"(TPM+1)"~"]")
ylab = "GSVA Score [Proteome]"

main = sprintf("SCC in TPM & Proteome [GSVA, %s]", sample_gsva)
# Ex_GSVA %>% na.omit %>% dim   # 292
Ex_GSVA %>% plot_def(TPM, Proteome, main=main, 
                     xlab=xlab, ylab=ylab, alpha=0.5, width=12, height=12, save=T, xy_line=T)

# idx = match(Corr_TPM_Pt_Ori$Sample, names(col_raw))
# identical(Corr_TPM_Pt_Ori$Corr, col_raw[idx])   # F 
# plot(Corr_TPM_Pt_Ori$Corr, col_raw[idx])



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
  
  
  ### [Source Data] Fig. 8
  dir = "../GCNPath/_case_study"
  file = sprintf("%s/SCLC_in_vitro.csv", dir)
  Pred_In_Vitro = read.csv(file)
  
  models = c("TGDRP_SG", "DRPreter_SG", "GCNPath")
  drugs = c("Etoposide", "Cisplatin", "Anlotinib", "Alisertib", "Barasertib", "AMG-900")
  
  Pred_List_TN_ = Pred_List_TN %>% 
    mutate(Model=Model %>% factor(levels=models), 
           Drug=Drug %>% factor(levels=drugs)) %>% arrange(Model, Drug, Sample)
  
  Pred_GCNPath_ = Pred_GCNPath_ %>% arrange(desc(Corr_TPM_Prot))
  Corr_TPM_Pt_ = Corr_TPM_Pt %>% arrange(desc(PCC_GSVA))
  
  Temp = list(Pred_In_Vitro, Pred_List_TN_, 
              Pred_GCNPath_, Corr_TPM_Pt_, Ex_Raw, Ex_GSVA)
  
  num_fig = c("a", "b", "c1", "c2", "c3", "c4")
  Temp %>% save_for_nc(num=8, suppl=F, num_fig=num_fig)
  
  
  ### [Source Data] Supplementary Fig. 37
  Pred_List_TN_CB_ = Pred_List_TN_CB %>% 
    mutate(Model=Model %>% factor(levels=models), 
           Drug=Drug %>% factor(levels=drugs)) %>% arrange(Model, Drug, Sample)
  
  Pred_List_TN_CB_ %>% save_for_nc(num=37, suppl=T)
}

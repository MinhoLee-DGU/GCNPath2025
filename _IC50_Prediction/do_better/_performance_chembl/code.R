#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

source("../functions.R")
loadings()


# cf. Gather cell list trained in all models

# [GDSC] MUT+CNV (GraphDRP, tCNNS)
file = "../GraphDRP/cell_feat.csv"
cells_graphdrp = fread_def(file) %>% rownames   # 990

# [GDSC] EXP (HiDRA)
file = "../HiDRA/_data/GDSC_RNA_Array.csv"
cells_hidra = fread_def(file) %>% rownames   # 968

# [CCLE] EXP (DRPreter)
file = "../DRPreter/_data/EXP.csv"
cells_drpreter = fread_def(file) %>% rownames   # 1389

# [CCLE] MUT+EXP+CNV (TGDRP, TGSA)
file1 = "../TGSA/_data/MUT.csv"
file2 = "../TGSA/_data/EXP.csv"
file3 = "../TGSA/_data/CNV.csv"

cells_tgsa_1 = fread_def(file1, header=T) %>% rownames   # 1399
cells_tgsa_2 = fread_def(file2, header=T) %>% rownames   # 1399
cells_tgsa_3 = fread_def(file3, header=T) %>% rownames   # 1399
cells_tgsa = Reduce(intersect, list(cells_tgsa_1, cells_tgsa_2, cells_tgsa_3))   # 1399

# [SANGER] EXP (GCNPath, DRPreter_SG, PaccMann_SG)
file = "../GCNPath/processed/cell_data_kegg/SANGER_RNA_GSVA.csv"
cells_gcnpath = fread_def(file) %>% rownames   # 1431

# [SANGER] MUT+EXP+CNV (TGDRP_SG, TGSA_SG)
file1 = "../TGSA_SANGER/_data/MUT.csv"
file2 = "../TGSA_SANGER/_data/EXP.csv"
file3 = "../TGSA_SANGER/_data/CNV.csv"

cells_tgsa_sg_1 = fread_def(file1, header=T) %>% rownames   # 1183
cells_tgsa_sg_2 = fread_def(file2, header=T) %>% rownames   # 1183
cells_tgsa_sg_3 = fread_def(file3, header=T) %>% rownames   # 1183
cells_tgsa_sg = Reduce(intersect, list(cells_tgsa_sg_1, cells_tgsa_sg_2, cells_tgsa_sg_3))   # 1183

# [Intersect] Cells commonly utilized in training...
file = "../../processed_data/ic50_data/GDSC/IC50_GDSC.csv"
IC50_GDSC = read.csv(file)   # 373681

cells_cosmic = Reduce(intersect, list(cells_graphdrp, cells_hidra))    # 962
cells_broad = Reduce(intersect, list(cells_drpreter, cells_tgsa))      # 1337
cells_sanger = Reduce(intersect, list(cells_gcnpath, cells_tgsa_sg))   # 1183

cells_train = IC50_GDSC %>% 
  subset(COSMIC_ID %in% cells_cosmic) %>% 
  subset(BROAD_ID %in% cells_broad) %>% 
  subset(SANGER_MODEL_ID %in% cells_sanger) %>% 
  pull(SANGER_MODEL_ID) %>% unique   # 670



##### 2. Read prediction files

read_pred = function(dir, pattern, sep=",", model_name=NULL, pred_mode="regression") {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_[^CCLE]", "Pred_CCLE"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      try({
        Pred_TP = fread(df_name, header=T, sep=sep)
        if (pred_mode=="regression") {
          Pred_TP = Pred_TP[Standard_Relation=="=", ]
        } else if (pred_mode=="classification") {
          Pred_TP = Pred_TP[!(Standard_Relation %in% c("=", "~")), ]
        } else {
          Pred_TP = Pred_TP[Standard_Relation!="~", ]
        }
        
        df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
        nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
        model = strsplit(dir, "/")[[1]][2]
        
        Pred_TP$Seed = nth
        Pred_TP$Model = model
        Pred = Pred %>% rbind(Pred_TP)
      })
    }
    
    sprintf("# Number of Predictions : %s", nrow(Pred)) %>% print
    if (!is.null(model_name)) Pred$Model = model_name
    print("# Number of Predictions in each seed : ")
    Pred$Seed %>% table %>% print
    return(Pred)
  }
}

filter_invalid = function(Pred, cells=NULL, drugs=NULL, assays=NULL) {
  n_before = Pred %>% nrow
  Pred = Pred[!(Cell_ChEMBL_ID %in% cells) & !(Molecule_ChEMBL_ID %in% drugs) & !(Assay_ChEMBL_ID %in% assays), ]
  sprintf("# Invalid data excluded : %s > %s", n_before, nrow(Pred)) %>% print
  return(Pred)
}

calc_perf = function(Pred, sum_seed=F, test_db=F, cells_no=NULL, drugs_no=NULL, assays_no=NULL) {
  
  col_by = "Model"
  if (!sum_seed) col_by = col_by %>% c("Seed")
  if (test_db) col_by = col_by %>% c("Test_DB")
  
  Pred = Pred %>% filter_invalid(cells_no, drugs_no, assays_no)
  pred_inf = Pred$Prediction %>% is.infinite
  
  if (any(pred_inf)) {
    seed_inf = Pred[pred_inf, "Seed"] %>% unique %>% unlist %>% sort %>% paste(collapse=" & ")
    Pred = Pred[!pred_inf, ]
    sprintf("Prediction Inf found... [n=%s in Seed %s]", sum(pred_inf), seed_inf) %>% print
  }
  
  Perf = Pred %>% 
    group_by(across(all_of(col_by))) %>% 
    summarize(N_Test = n(),
              RMSE = RMSE(LN_IC50, Prediction),
              MAE = MAE(LN_IC50, Prediction), 
              R2 = R2(LN_IC50, Prediction),
              PCC = cor(LN_IC50, Prediction), 
              SCC = cor(LN_IC50, Prediction, method="spearman"))
  
  Perf = Perf %>% as.data.frame
  if (test_db) Perf = Perf %>% arrange(Model, Test_DB, Seed) %>% as.data.frame
  return(Perf)
}

pattern = "pred_chembl_([0-9]+).csv"
pattern_tcnns = "pred_chembl_seed([0-9]+)_([0-9]+).csv"

dir = "../tCNNS/Results/IC50_GDSC/Normal"
Pred_tCNNS = read_pred(dir, pattern=pattern_tcnns, model_name="tCNNS")   # 2636420

dir = "../HiDRA/Results/IC50_GDSC/Normal"
Pred_HiDRA = read_pred(dir, pattern=pattern, model_name="HiDRA")   # 2576790

# # Decided to exclude PaccMann due to significantly low number of available cell data...
# dir = "../PaccMann/Results/IC50_GDSC/Normal"
# Pred_PaccMann = read_pred(dir, pattern=pattern, model_name="PaccMann")

dir = "../PaccMann_SANGER/Results/IC50_GDSC/Normal"
Pred_PaccMann_SG = read_pred(dir, pattern=pattern, model_name="PaccMann_SG")   # 2899080

dir = "../GraphDRP/Results/IC50_GDSC/Normal"
Pred_GraphDRP = read_pred(dir, pattern=pattern, model_name="GraphDRP")   # 2636560

dir = "../TGSA/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP = read_pred(dir, pattern=pattern, model_name="TGDRP")   # 2765620

dir = "../TGSA/Results/IC50_GDSC/Normal/TGSA"
Pred_TGSA = read_pred(dir, pattern=pattern, model_name="TGSA")   # 2463930

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP_SG = read_pred(dir, pattern=pattern, model_name="TGDRP_SG")   # 2668330

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGSA"
Pred_TGSA_SG = read_pred(dir, pattern=pattern, model_name="TGSA_SG")   # 2635410

dir = "../DRPreter/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter = read_pred(dir, pattern=pattern, model_name="DRPreter")   # 2775610

dir = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG = read_pred(dir, pattern=pattern, model_name="DRPreter_SG")   # 2899070

dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath = read_pred(dir, pattern=pattern, model_name="GCNPath")   # 2899070

col = colnames(Pred_GCNPath)
Pred_ChEMBL = Reduce(rbind, list(Pred_tCNNS[, col, with=F], Pred_HiDRA,
                                 Pred_PaccMann_SG, Pred_GraphDRP[, col, with=F],
                                 Pred_TGDRP, Pred_TGSA, Pred_TGDRP_SG, Pred_TGSA_SG,
                                 Pred_DRPreter, Pred_DRPreter_SG, Pred_GCNPath))   # 29855890

model_list = c("tCNNS", "HiDRA", "PaccMann_SG", "GraphDRP", 
               "TGDRP", "TGSA", "TGDRP_SG", "TGSA_SG",
               "DRPreter", "DRPreter_SG", "GCNPath")

Pred_ChEMBL$Model = Pred_ChEMBL$Model %>% factor(levels=model_list)
# pred_names = sprintf("Pred_%s", model_list)
# rm(list=ls()[ls() %in% pred_names])



##### Save cell & drug pairs intersecting in all models

filter_pair = function(Pred, model_list=NULL, tag_ref=NULL, 
                       cells_no=NULL, drugs_no=NULL, assays_no=NULL, seed_ex=2030) {
  
  # Pick up cell x drug predicted in all models...
  
  if (!is.null(model_list)) {
    n_model = length(model_list)
  } else n_model = Pred$Model %>% unique %>% length
  
  Pred = Pred %>% filter_invalid(cells_no, drugs_no, assays_no)
  Pred[, Tag := paste0(Cell_ChEMBL_ID, "@", Molecule_ChEMBL_ID)]
  if (!is.null(tag_ref)) Pred = Pred[Tag %in% tag_ref, ]
  
  if (!is.null(seed_ex)) {
    Pred_ = Pred[Seed==seed_ex]
  } else Pred_ = Pred
  Tag_Num = Pred_[, .(Num=uniqueN(Model)), by=Tag]
  
  tag_int = Tag_Num %>% subset(Num==n_model) %>% pull(Tag)
  sprintf("# Cell x Drug predicted in all models : %s [Models=%s]", length(tag_int), n_model) %>% print
  Pred = Pred[Tag %in% tag_int, ]
  Pred[, Tag:=NULL]
  return(Pred)
}

count_num = function(Pred, seed_ex=2030) {
  if (!is.null(seed_ex)) Pred = Pred[Seed==seed_ex]
  Stat_Num = Pred %>% group_by(Model) %>%
    summarise(Cell=ulen(Cell_ChEMBL_ID), 
              Drug=ulen(Molecule_ChEMBL_ID), 
              IC50=n()) %>% as.data.frame
  return(Stat_Num)
}

Pred_ChEMBL = Pred_ChEMBL %>% filter_pair(model_list=model_list)   # 29855890 > 26432340
Pred_ChEMBL = Pred_ChEMBL %>% subset(SANGER_MODEL_ID %in% cells_train & !Drug_GDSC)   # 26432340 > 26141060
Pred_Num = Pred_ChEMBL %>% count_num(seed_ex=2030)   # 237646



##### 3. Performance [Regression]

Perf_ChEMBL = Pred_ChEMBL %>% calc_perf(sum_seed=F)
col = c("RMSE", "MAE", "R2", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Avg = Perf_ChEMBL %>% group_by(Model) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame


select_best_seed = function(Pred, Perf_Best, test_db=F) {
  Pred_Best = list()
  for (i in 1:nrow(Perf_Best)) {
    model = Perf_Best$Model[i]
    seed = Perf_Best$Seed[i]
    Pred_Best[[i]] = Pred[Model==model & Seed==seed, ]
    
    if (test_db) {
      db = Perf_Best$Test_DB[i]
      Pred_Best[[i]] = Pred_Best[[i]][Test_DB==db, ]
    }
  }
  Pred_Best = data.table::rbindlist(Pred_Best)
  return(Pred_Best)
}

Perf_ChEMBL_Best = Perf_ChEMBL %>% group_by(Model) %>% filter(RMSE==min(RMSE)) %>% as.data.frame
Pred_ChEMBL_Best = Pred_ChEMBL %>% select_best_seed(Perf_Best=Perf_ChEMBL_Best)   # 2614106

draw_pred = F
dir = mkdir("Performance")
xlab = bquote(Actual~ln(IC[50]))
ylab = bquote(Predicted~ln(IC[50]))

if (draw_pred) {
  for (model in model_list) {
    main = sprintf("%s/Prediction [%s]", dir, model)
    Pred_ChEMBL_Best %>% subset(Model==model) %>%
      plot_def(LN_IC50, Prediction, main=main,
               xlab=xlab, ylab=ylab, alpha=0.5,
               axis_tl=27, axis_tx=22.5, width=15, height=15, dpi=1200,
               xy_line=T, force_bold=F, raster=T, save=T, save_svg=T)
  }
}

file = "Performance.xlsx"
sheet = c("Perf", "Perf_Avg", "Perf_Best")
df_list = list(Perf_ChEMBL, Perf_ChEMBL_Avg, Perf_ChEMBL_Best)
write.xlsx(df_list, file=file, rowNames=F, sheetName=sheet)


is_inf_def = function(x) {
  inf_p = sum(x>0 & is.infinite(x))
  inf_m = sum(x<0 & is.infinite(x))
  sprintf("# Total : %s", length(x)) %>% print
  if (inf_p!=0) sprintf("# Infinite (+) : %s", inf_p) %>% print
  if (inf_m!=0) sprintf("# Infinite (-) : %s", inf_m) %>% print
  sprintf("# Infinite (Total) : %s", sum(is.infinite(x))) %>% print
}

Pred_ChEMBL %>% subset(Model=="tCNNS") %>% pull(Prediction) %>% is_inf_def
# Infinite (+) : 221745
# Infinite (-) : 162343

Pred_ChEMBL_Best %>% subset(Model=="tCNNS") %>% pull(Prediction) %>% is_inf_def
# Infinite (+) : 29234
# Infinite (-) : 2



##### 4. Read prediction files [Batch Effect of RNA Data]

pattern_gdsc = "pred_chembl_gdsc_seed([0-9]+).csv"
pattern_ccle = "pred_chembl_ccle_seed([0-9]+).csv"
pattern_sanger = "pred_chembl_sanger_seed([0-9]+).csv"

dir = "../HiDRA/Results/IC50_GDSC/Normal"
Pred_HiDRA_C = read_pred(dir, pattern=pattern_ccle, model_name="HiDRA")     # 2776400
Pred_HiDRA_S = read_pred(dir, pattern=pattern_sanger, model_name="HiDRA")   # 2899080

dir = "../PaccMann_SANGER/Results/IC50_GDSC/Normal"
Pred_PaccMann_SG_G = read_pred(dir, pattern=pattern_gdsc, model_name="PaccMann_SG")   # 2576790
Pred_PaccMann_SG_C = read_pred(dir, pattern=pattern_ccle, model_name="PaccMann_SG")   # 2775620

dir = "../TGSA/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP_G = read_pred(dir, pattern=pattern_gdsc, model_name="TGDRP")     # 2404570
Pred_TGDRP_S = read_pred(dir, pattern=pattern_sanger, model_name="TGDRP")   # 2727030

dir = "../TGSA/Results/IC50_GDSC/Normal/TGSA"
Pred_TGSA_G = read_pred(dir, pattern=pattern_gdsc, model_name="TGSA")     # 2404310
Pred_TGSA_S = read_pred(dir, pattern=pattern_sanger, model_name="TGSA")   # 2463930

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Pred_TGDRP_SG_G = read_pred(dir, pattern=pattern_gdsc, model_name="TGDRP_SG")   # 2575710
Pred_TGDRP_SG_C = read_pred(dir, pattern=pattern_ccle, model_name="TGDRP_SG")   # 2497220

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGSA"
Pred_TGSA_SG_G = read_pred(dir, pattern=pattern_gdsc, model_name="TGSA_SG")   # 2575640
Pred_TGSA_SG_C = read_pred(dir, pattern=pattern_ccle, model_name="TGSA_SG")   # 2464300

dir = "../DRPreter/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_G = read_pred(dir, pattern=pattern_gdsc, model_name="DRPreter")     # 2576780
Pred_DRPreter_S = read_pred(dir, pattern=pattern_sanger, model_name="DRPreter")   # 2899070

dir = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Pred_DRPreter_SG_G = read_pred(dir, pattern=pattern_gdsc, model_name="DRPreter_SG")   # 2576780
Pred_DRPreter_SG_C = read_pred(dir, pattern=pattern_ccle, model_name="DRPreter_SG")   # 2775610

dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Pred_GCNPath_G = read_pred(dir, pattern=pattern_gdsc, model_name="GCNPath")   # 2585110
Pred_GCNPath_C = read_pred(dir, pattern=pattern_ccle, model_name="GCNPath")   # 2776390

Pred_ChEMBL_G = Reduce(rbind, list(Pred_PaccMann_SG_G, Pred_TGDRP_G, Pred_TGSA_G, 
                                   Pred_TGDRP_SG_G, Pred_TGSA_SG_G,
                                   Pred_DRPreter_G, Pred_DRPreter_SG_G, Pred_GCNPath_G))   # 20275690

Pred_ChEMBL_C = Reduce(rbind, list(Pred_HiDRA_C, Pred_PaccMann_SG_C, 
                                   Pred_TGDRP_SG_C, Pred_TGSA_SG_C,
                                   Pred_DRPreter_SG_C, Pred_GCNPath_C))   # 16065540

Pred_ChEMBL_S = Reduce(rbind, list(Pred_HiDRA_S, Pred_TGDRP_S, Pred_TGSA_S, Pred_DRPreter_S))   # 10989110


model_list_exp = c("HiDRA", "PaccMann_SG", 
                   "TGDRP", "TGSA", "TGDRP_SG", "TGSA_SG",
                   "DRPreter", "DRPreter_SG", "GCNPath")

pred_names_g = sprintf("Pred_%s_G", model_list_exp)
pred_names_c = sprintf("Pred_%s_C", model_list_exp)
pred_names_s = sprintf("Pred_%s_S", model_list_exp)
# pred_names_exp = c(pred_names_g, pred_names_c, pred_names_s)
# rm(list=ls()[ls() %in% pred_names_exp])

Pred_ChEMBL_G$Model = Pred_ChEMBL_G$Model %>% factor(levels=model_list)
Pred_ChEMBL_C$Model = Pred_ChEMBL_C$Model %>% factor(levels=model_list)
Pred_ChEMBL_S$Model = Pred_ChEMBL_S$Model %>% factor(levels=model_list)

Pred_ChEMBL_G$Test_DB = "GDSC"
Pred_ChEMBL_C$Test_DB = "CCLE"
Pred_ChEMBL_S$Test_DB = "SANGER"

# Get only Cell x Drug combinations from original test
tag_ref = Pred_ChEMBL %>% with(paste0(Cell_ChEMBL_ID, "@", Molecule_ChEMBL_ID)) %>% unique   # 237646

Pred_ChEMBL_G = Pred_ChEMBL_G %>% filter_pair(model_list=NULL, tag_ref=tag_ref)   # 20275690 > 19011680
Pred_ChEMBL_C = Pred_ChEMBL_C %>% filter_pair(model_list=NULL, tag_ref=tag_ref)   # 16065540 > 14258760
Pred_ChEMBL_S = Pred_ChEMBL_S %>% filter_pair(model_list=NULL, tag_ref=tag_ref)   # 10989110 > 9505840

count_num_db = function(Pred, seed_ex=2030) {
  ulen = function(x) length(unique(x))
  if (!is.null(seed_ex)) Pred = Pred[Seed==seed_ex]
  
  Stat_Num = Pred %>% 
    group_by(Model, Test_DB) %>%
    summarise(Cell=ulen(Cell_ChEMBL_ID), 
              Drug=ulen(Molecule_ChEMBL_ID), 
              IC50=n()) %>% as.data.frame
  
  return(Stat_Num)
}

Pred_ChEMBL_Ext = Reduce(rbind, list(Pred_ChEMBL_G, Pred_ChEMBL_C, Pred_ChEMBL_S))   # 42776280
Pred_Num_Ext = Pred_ChEMBL_Ext %>% count_num_db(seed_ex=2030)

levels = c("GDSC", "CCLE", "SANGER")
Pred_ChEMBL_Ext$Test_DB = Pred_ChEMBL_Ext$Test_DB %>% factor(levels=levels)
Pred_Num_Ext = Pred_Num_Ext %>% 
  mutate(Test_DB=Test_DB %>% factor(levels=levels)) %>% 
  arrange(Model, Test_DB) %>% as.data.frame

# Fortunately, All cell x drug combinations from original test exist in external test
Pred_Num$IC50 %>% unique       # 237646
Pred_Num_Ext$IC50 %>% unique   # 237646



##### 5. Performance [Regression, Batch Effect of RNA Data]

Perf_ChEMBL_Ext = Pred_ChEMBL_Ext %>% 
  calc_perf(sum_seed=F, test_db=T)   # 180 [9 x 2 x 10]

col = c("RMSE", "MAE", "R2", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Ext_Avg = Perf_ChEMBL_Ext %>% group_by(Model, Test_DB) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame   # 18

Perf_ChEMBL_Ext_Best = Perf_ChEMBL_Ext %>% 
  group_by(Model, Test_DB) %>% filter(RMSE==min(RMSE)) %>% as.data.frame   # 18
Pred_ChEMBL_Ext_Best = Pred_ChEMBL_Ext %>% 
  select_best_seed(Perf_Best=Perf_ChEMBL_Ext_Best, test_db=T)   # 42776280 > 4277628


# Boxplot
model_g = c("HiDRA")
model_c = c("TGDRP", "TGSA", "DRPreter")
model_s = c("PaccMann_SG", "TGDRP_SG", "TGSA_SG", "DRPreter_SG", "GCNPath")

model = c(model_g, model_c, model_s)
test_db = c(rep("GDSC", length(model_g)), rep("CCLE", length(model_c)), rep("SANGER", length(model_s)))
Train_DB = data.frame(Model=model, Test_DB=test_db)

Perf_ChEMBL_Total = Perf_ChEMBL %>% subset(Model %in% Train_DB$Model)
idx = match(Perf_ChEMBL_Total$Model, Train_DB$Model)

Perf_ChEMBL_Total = Perf_ChEMBL_Total %>% 
  mutate(Test_DB=Train_DB$Test_DB[idx], Train_DB=T) %>% 
  relocate(Test_DB, Train_DB, .after=Seed) %>% as.data.frame

Perf_ChEMBL_Ext = Perf_ChEMBL_Ext %>% 
  mutate(Train_DB=F) %>% relocate(Train_DB, .after=Test_DB) %>% as.data.frame

levels = c("GDSC", "CCLE", "SANGER")
Perf_ChEMBL_Total = rbind(Perf_ChEMBL_Total, Perf_ChEMBL_Ext)
Perf_ChEMBL_Total$Test_DB = Perf_ChEMBL_Total$Test_DB %>% factor(levels=levels)


append_def = function(A, ...) {
  add_list = list(...)
  A = c(A, add_list)
  return(A)
}

plot_perf = function(Perf_List, score="RMSE", dir=".", models_no=NULL, 
                     col_highlight="firebrick2", test_db=T, scale_y=F, 
                     axis_tl=36, axis_tx=24, legend_tl=27, legend_tx=27, size_point=5,
                     alpha=0.9, margin=0.6, pos_dodge=0.8, width=30, height=20, save=F) {
  
  # [score] RMSE, PCC, R2
  suppressMessages(library(ggtext))
  pos = position_dodge(pos_dodge)
  if (!is.null(models_no)) models_no = paste(models_no, collapse=" & ")
  
  highlight = function(x, pattern, color="black", family="") {
    suppressMessages(library(glue))
    suppressMessages(library(ggtext))
    ifelse(grepl(pattern, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
  }
  
  highlight_x = function(x) highlight(x, "GCNPath", color=col_highlight, family="bold")
  add_common = list(scale_x_discrete(labels=highlight_x), theme(axis.text.x=element_markdown()))
  
  if (test_db) {
    Perf_Sum = Perf_List %>% group_by(Model, Test_DB, Train_DB) %>% 
      summarise(Median=median(object(score))) %>% as.data.frame
    Perf_Sum$Median[!Perf_Sum$Train_DB] = NA
    
    pos_mark = position_dodge2(width=pos_dodge, preserve="single")
    pl_mark1 = geom_point(data=Perf_Sum, aes(Model, Median, fill=Test_DB), 
                          color="grey30", shape=17, size=size_point, position=pos_mark, show.legend=F)
    pl_mark2 = geom_point(data=Perf_Sum, aes(Model, Median, fill=Test_DB),
                          color="springgreen", shape=17, size=size_point-1, position=pos_mark, show.legend=F)
    
    ylim = range(Perf_List[, score])
    ylim = c(ylim[1]-0.1, ylim[2]+0.1)
    color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
    scale_y_axis = coord_trans(y="sqrt", ylim=ylim)
    
    pos_legend = theme(legend.position="bottom")
    margin1 = margin(0, 1.0, 0, 1.0, unit="cm")
    margin2 = margin(0, 0.5, 0, 0.25, unit="cm")
    margin_lg = theme(legend.title=element_text(size=legend_tl, margin=margin1), 
                      legend.text=element_text(size=legend_tx, margin=margin2))
    
    add = list(scale_fill_manual(values=color), pl_mark1, pl_mark2, pos_legend, margin_lg)
    add = add %>% c(add_common)
    if (scale_y) add = add %>% append_def(scale_y_axis)
    
    if (is.null(models_no)) {
      main = sprintf("%s/Performances [%s, RNA Batch Effect]", dir, score)
    } else main = sprintf("%s/Performances [%s, RNA Batch Effect, Except %s]", dir, score, models_no)
    
    Perf_List %>% boxplot_def(Model, object(score), fill=Test_DB,
                              main=main, ylab=score, alpha=alpha, pos=pos, legend="Dataset",
                              axis_tl=axis_tl, axis_tx=axis_tx, vjust=1, hjust=1, add=add,
                              legend_tl=legend_tl, legend_tx=legend_tx, margin=margin,
                              width=width, height=height, save=save, save_svg=T)
    
  } else {
    if (is.null(models_no)) {
      main = sprintf("%s/Performances [%s]", dir, score)
    } else main = sprintf("%s/Performances [%s, Except %s]", dir, score, models_no)
    
    Perf_List %>% boxplot_def(Model, object(score), legend=F,
                              main=main, ylab=score, alpha=alpha, pos=pos, add=add_common,
                              axis_tl=axis_tl, axis_tx=axis_tx, vjust=1, hjust=1, margin=margin,
                              width=width, height=height, save=save, save_svg=T)
  }
}

plot_perf_bar = function(Perf_List, score="RMSE", col_highlight="firebrick2",  
                         dir=NULL, test_db=T, point=F, pos_dodge=0.8, 
                         axis_tl=36, axis_tx=24, legend_tl=22.5, legend_tx=22.5, 
                         size_point=4.5, width=18, height=15, reorder_score=F, models_no=NULL, save=T) {
  
  # width=30, height=18
  suppressMessages(library(ggtext))
  suppressMessages(library(ggpubr))
  
  highlight = function(x, pattern, color="black", family="") {
    suppressMessages(library(glue))
    suppressMessages(library(ggtext))
    ifelse(grepl(pattern, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
  }
  
  pos = position_dodge(pos_dodge)
  if (!is.null(models_no)) models_no = paste(models_no, collapse=" & ")
  
  font1 = font("ylab", size=axis_tl, margin=margin(r=10, unit="pt"))
  font2 = font("x.text", color="grey30", size=axis_tx, margin=margin(t=10, unit="pt"))
  font3 = font("y.text", color="grey30", size=axis_tx, margin=margin(r=10, unit="pt"))
  font4 = font("legend.title", size=legend_tl)
  font5 = font("legend.text", size=legend_tx)
  
  font = font1 + font2 + font3 + font4 + font5
  color = RColorBrewer::brewer.pal(8, "Reds")[c(8, 6, 2)]
  main = sprintf("%s/Performance [%s]", dir, score)
  # levels(Perf_List$Model)[levels(Perf_List$Model)=="GCNPath"] = "<b>GCNPath</b>"
  
  if (!is.null(models_no)) {
    Perf_List = Perf_List %>% subset(!(Model %in% models_no))
    Perf_List$Model = Perf_List$Model %>% base::droplevels()
    main = sprintf("%s/Performance [%s, (Except)]", dir, score)
  }
  
  if (reorder_score) {
    lvl_models = Perf_List %>% group_by(Model) %>%
      summarise(Mean=mean(object(score), na.rm=T)) %>%
      arrange(desc(Mean)) %>% pull(Model)
    if (!(score %in% c("RMSE", "MAE"))) {
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
  
  if (score %in% c("RMSE", "MAE")) {
    ymin = min(yval, na.rm=T)-0.002
    ymax = max(yval, na.rm=T)+0.002
  } else {
    ymin = max(min(yval, na.rm=T)-0.002, -1)
    ymax = min(max(yval, na.rm=T)+0.002, 1)
  }
  
  if (test_db) {
    Perf_Sum = Perf_List %>% group_by(Model, Test_DB, Train_DB) %>% 
      summarise(Median=median(object(score))) %>% as.data.frame
    Perf_Sum$Median[!Perf_Sum$Train_DB] = NA
    
    pos_mark = position_dodge2(width=pos_dodge, preserve="single")
    pl_mark1 = geom_point(data=Perf_Sum, aes(Model, Median, fill=Test_DB), 
                          color="grey30", shape=17, size=size_point, position=pos_mark, show.legend=F)
    pl_mark2 = geom_point(data=Perf_Sum, aes(Model, Median, fill=Test_DB),
                          color="springgreen", shape=17, size=size_point-1, position=pos_mark, show.legend=F)
    
    pl = Perf_List %>%
      ggbarplot(x="Model", y=score, color="black", fill="Test_DB",
                xlab=F, add="mean_se", position=pos) +
      scale_fill_manual(values=color) + pl_mark1 + pl_mark2
  } else {
    pl = Perf_List %>%
      ggbarplot(x="Model", y=score, color="black", fill="#539ed6",
                xlab=F, add="mean_se", position=pos)
    # [V0] fill=color[1]
    # [V1] fill="#ff4d59"
    # [V2] fill="#539ed6"
  }
  
  pl = pl + font + rotate_x_text(30)
  pl = pl %>% ggpar(ylim=c(ymin, ymax))
  pl = pl + scale_x_discrete(labels=highlight_x) + theme(axis.text.x=element_markdown())
  if (point) pl = pl + geom_point(alpha=0.5)
  
  if (save) {
    pl %>% save_fig_ggpubr(main=main, width=width, height=height, svg=T)
  } else print(pl)
}

score = c("RMSE", "PCC", "SCC")
dir = mkdir("Performances [Total]")

Perf_ChEMBL %>% plot_perf(score=score[1], dir=dir, save=T, test_db=F,
                          axis_tl=27, axis_tx=16.5, width=20, height=15)
Perf_ChEMBL %>% plot_perf(score=score[2], dir=dir, save=T, test_db=F, 
                          axis_tl=27, axis_tx=16.5, width=20, height=15)
Perf_ChEMBL %>% plot_perf(score=score[3], dir=dir, save=T, test_db=F, 
                          axis_tl=27, axis_tx=16.5, width=20, height=15)

# dir = mkdir("Performances [Total, Barplot]")
# Perf_ChEMBL %>% plot_perf_bar(score=score[1], dir=dir, save=T, test_db=F, width=20, height=15)
# Perf_ChEMBL %>% plot_perf_bar(score=score[2], dir=dir, save=T, test_db=F, width=20, height=15)
# Perf_ChEMBL %>% plot_perf_bar(score=score[3], dir=dir, save=T, test_db=F, width=20, height=15)


models_no = c("tCNNS", "PaccMann_SG")
models_no_ = paste(models_no, collapse=" & ")
dir = mkdir(sprintf("Performances [Total, Except %s]", models_no_))

Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[1], dir=dir, save=T, test_db=F, 
            models_no=models_no, axis_tl=27, axis_tx=18, width=20, height=15)
Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[2], dir=dir, save=T, test_db=F, 
            models_no=models_no, axis_tl=27, axis_tx=18, width=20, height=15)
Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[3], dir=dir, save=T, test_db=F, 
            models_no=models_no, axis_tl=27, axis_tx=18, width=20, height=15)

# dir = mkdir(sprintf("Performances [Total, Except %s, Barplot]", models_no_))
# 
# Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
#   plot_perf_bar(score=score[1], dir=dir, save=T, test_db=F, 
#                 models_no=models_no, width=20, height=15)
# Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
#   plot_perf_bar(score=score[2], dir=dir, save=T, test_db=F, 
#                 models_no=models_no, width=20, height=15)
# Perf_ChEMBL %>% subset(!(Model %in% models_no)) %>%
#   plot_perf_bar(score=score[3], dir=dir, save=T, test_db=F, 
#                 models_no=models_no, width=20, height=15)


dir = mkdir("Performances [RNA Batch Effect]")
Perf_ChEMBL_Total %>% plot_perf(score=score[1], dir=dir, width=27, height=20, save=T, test_db=T)
Perf_ChEMBL_Total %>% plot_perf(score=score[2], dir=dir, width=27, height=20, save=T, test_db=T)
Perf_ChEMBL_Total %>% plot_perf(score=score[3], dir=dir, width=27, height=20, save=T, test_db=T)

# dir = mkdir("Performances [RNA Batch Effect, Barplot]")
# Perf_ChEMBL_Total %>% plot_perf_bar(score=score[1], axis_tl=24, axis_tx=20, dir=dir, save=T, test_db=T)
# Perf_ChEMBL_Total %>% plot_perf_bar(score=score[2], axis_tl=24, axis_tx=20, dir=dir, save=T, test_db=T)
# Perf_ChEMBL_Total %>% plot_perf_bar(score=score[3], axis_tl=24, axis_tx=20, dir=dir, save=T, test_db=T)

# models_no = c("HiDRA", "PaccMann_SG")
models_no = "PaccMann_SG"
models_no_ = paste(models_no, collapse=" & ")
dir = mkdir(sprintf("Performances [RNA Batch Effect, Except %s]", models_no_))

Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[1], dir=dir, save=T, test_db=T, 
            width=32, axis_tx=30, models_no=models_no, size_point=6)
Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[2], dir=dir, save=T, test_db=T, 
            width=32, axis_tx=30, models_no=models_no, size_point=6)
Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>%
  plot_perf(score=score[3], dir=dir, save=T, test_db=T, 
            width=32, axis_tx=30, models_no=models_no, size_point=6)

# dir = mkdir(sprintf("Performances [RNA Batch Effect, Except %s, Barplot]", models_no_))
# 
# Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>% 
#   plot_perf_bar(score=score[1], dir=dir, save=T, test_db=T, 
#                 width=24, axis_tl=24, axis_tx=24, models_no=models_no, size_point=6)
# Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>% 
#   plot_perf_bar(score=score[2], dir=dir, save=T, test_db=T, 
#                 width=24, axis_tl=24, axis_tx=24, models_no=models_no, size_point=6)
# Perf_ChEMBL_Total %>% subset(!(Model %in% models_no)) %>% 
#   plot_perf_bar(score=score[3], dir=dir, save=T, test_db=T, 
#                 width=24, axis_tl=24, axis_tx=24, models_no=models_no, size_point=6)

col = c("RMSE", "MAE", "R2", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Total_Avg = Perf_ChEMBL_Total %>% group_by(Model, Test_DB) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame   # 18

Perf_ChEMBL_Total_Best = Perf_ChEMBL_Total %>% 
  group_by(Model, Test_DB) %>% filter(RMSE==min(RMSE)) %>% 
  arrange(Model, Test_DB) %>% as.data.frame   # 18

file = "Performance [Total].xlsx"
sheet = c("Perf", "Perf_Avg", "Perf_Best")
df_list = list(Perf_ChEMBL_Total, Perf_ChEMBL_Total_Avg, Perf_ChEMBL_Total_Best)
write.xlsx(df_list, file=file, rowNames=F, sheetName=sheet)



##### 6. Performance [Regression, Batch Effect of RNA Data]

wilcox_test_sub = function(Perf_List, model_h1, model_h0) {
  
  Perf_Model_H1 = Perf_List %>% subset(Model==model_h1)
  Perf_Model_H0 = Perf_List %>% subset(Model==model_h0)
  
  utest_rmse = wilcox.test(Perf_Model_H1$RMSE, Perf_Model_H0$RMSE, alternative="less")
  utest_pcc = wilcox.test(Perf_Model_H1$PCC, Perf_Model_H0$PCC, alternative="greater")
  utest_scc = wilcox.test(Perf_Model_H1$SCC, Perf_Model_H0$SCC, alternative="greater")
  
  utest_res = c(model_h0, model_h1, 
                utest_rmse$alternative, utest_pcc$alternative, utest_scc$alternative, 
                utest_rmse$statistic, utest_pcc$statistic, utest_scc$statistic,
                utest_rmse$p.value, utest_pcc$p.value, utest_scc$p.value)
  
  return(utest_res)
}

wilcox_test = function(Perf_List, fdr_adjust=T) {
  
  Test_Res = data.frame()
  model_list = Perf_List$Model %>% unique
  Model_Pair = expand.grid(Model_H1=model_list, Model_H0=model_list)
  Model_Pair = Model_Pair %>% subset(Model_H1!=Model_H0)
  
  for (i in 1:nrow(Model_Pair)) {
    model_h1 = Model_Pair$Model_H1[i] %>% as.character
    model_h0 = Model_Pair$Model_H0[i] %>% as.character
    test_res = wilcox_test_sub(Perf_List, model_h1, model_h0)
    Test_Res = Test_Res %>% rbind(test_res)
  }
  
  colnames(Test_Res) = c("Model_H0", "Model_H1", 
                         "Hypothesis_RMSE", "Hypothesis_PCC", "Hypothesis_SCC", 
                         "W_RMSE", "W_PCC", "W_SCC", "Pval_RMSE", "Pval_PCC", "Pval_SCC")
  
  Test_Res[, 6:11] = Test_Res[, 6:11] %>% sapply(as.numeric)
  fdr_adj = function(x) p.adjust(x, "fdr")
  mlog10 = function(x) -log10(x)
  signif = function(x) x<0.05
  
  metrics = c("RMSE", "SCC", "PCC")
  col_fdr = sprintf("FDR_%s", metrics)
  col_pval = sprintf("Pval_%s", metrics)
  
  if (fdr_adjust) {
    Test_Res = Test_Res %>% 
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
                             round=2, width=24, height=18, asterisk=T, save=F) {
  
  # axis_tx=20, legend_tl=18, legend_tx=16.5, size=4.8,
  # size=4.8, round=2, width=36, height=30, 
  
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
    
    # legend1 = bquote(bold(atop(-log["10"](FDR), p(RMSE["y"]<RMSE["x"]))))
    # legend2 = bquote(bold(atop(-log["10"](FDR), p(PCC["y"]>PCC["x"]))))
    # legend3 = bquote(bold(atop(-log["10"](FDR), p(SCC["y"]>SCC["x"]))))
    
    Perf_UTest = Perf_UTest %>% 
      mutate(Pval_RMSE_Log=-log(FDR_RMSE, 10), 
             Pval_PCC_Log=-log(FDR_PCC, 10), 
             Pval_SCC_Log=-log(FDR_SCC, 10))
  } else {
    legend1 = bquote(atop("-Log"["10"](Pval), p(RMSE["y"]<RMSE["x"])))
    legend2 = bquote(atop("-Log"["10"](Pval), p(PCC["y"]>PCC["x"])))
    legend3 = bquote(atop("-Log"["10"](Pval), p(SCC["y"]>SCC["x"])))
    
    # legend1 = bquote(bold(atop(-log["10"](Pval), p(RMSE["y"]<RMSE["x"]))))
    # legend2 = bquote(bold(atop(-log["10"](Pval), p(PCC["y"]>PCC["x"]))))
    # legend3 = bquote(bold(atop(-log["10"](Pval), p(SCC["y"]>SCC["x"]))))
    
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
                          add=add, width=width+1.0, height=height, mean_summary=F, save=save)
  
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

rename_model = function(Perf_UTest, models, lvl_model=NULL) {
  Perf_UTest = Perf_UTest %>% 
    mutate(Model_H0=ifelse(Model_H0 %in% models, sprintf("%s*", Model_H0), Model_H0), 
           Model_H1=ifelse(Model_H1 %in% models, sprintf("%s*", Model_H1), Model_H1))
  
  if (!is.null(lvl_model)) {
    lvl_model = ifelse(lvl_model %in% models, sprintf("%s*", lvl_model), lvl_model)
    Perf_UTest = Perf_UTest %>% 
      mutate(Model_H0 = Model_H0 %>% factor(levels=lvl_model), 
             Model_H1 = Model_H1 %>% factor(levels=lvl_model))
  }
  return(Perf_UTest)
}

Perf_UTest = Perf_ChEMBL %>% wilcox_test(fdr_adjust=T)
Perf_UTest_G = Perf_ChEMBL_Total %>% subset(Test_DB=="GDSC") %>% wilcox_test(fdr_adjust=T)
Perf_UTest_C = Perf_ChEMBL_Total %>% subset(Test_DB=="CCLE") %>% wilcox_test(fdr_adjust=T)
Perf_UTest_S = Perf_ChEMBL_Total %>% subset(Test_DB=="SANGER") %>% wilcox_test(fdr_adjust=T)

Perf_UTest_G_ = Perf_UTest_G %>% rename_model(model_g, lvl_model=model_list_exp)
Perf_UTest_C_ = Perf_UTest_C %>% rename_model(model_c, lvl_model=model_list_exp)
Perf_UTest_S_ = Perf_UTest_S %>% rename_model(model_s, lvl_model=model_list_exp)


metrics = c("RMSE", "PCC", "SCC")
dir1 = mkdir("Performance [Wilcox Test, Grid]/All Models")
dir2 = mkdir("Performance [Wilcox Test, Grid]/GDSC")
dir3 = mkdir("Performance [Wilcox Test, Grid]/CCLE")
dir4 = mkdir("Performance [Wilcox Test, Grid]/SANGER")

main1 = sprintf("%s/%s Wilcox of Models", dir1, metrics)
main2 = sprintf("%s/%s Wilcox of Models [GDSC]", dir2, metrics)
main3 = sprintf("%s/%s Wilcox of Models [CCLE]", dir3, metrics)
main4 = sprintf("%s/%s Wilcox of Models [SANGER]", dir4, metrics)

Perf_UTest %>% wilcox_model_grid(main=main1, save=T, size=4, width=25.8, height=20, lvl_model=model_list)
Perf_UTest_G_ %>% wilcox_model_grid(main=main2, save=T)
Perf_UTest_C_ %>% wilcox_model_grid(main=main3, save=T)
Perf_UTest_S_ %>% wilcox_model_grid(main=main4, save=T)


# Tag for Re-Visualization of duplicate IC50 values [GDSC & ChEMBL]

Cell_Drug_Pair = data.frame(
  Cell_ChEMBL_ID = gsub("(.*)@(.*)", "\\1", tag_ref), 
  Drug_ChEMBL_ID = gsub("(.*)@(.*)", "\\2", tag_ref)
)

file = "Cell_Drug_Pair.csv"
fwrite(Cell_Drug_Pair, file=file)

dir = "Performance"
file = sprintf("%s/Prediction.csv", dir)
fwrite(Pred_ChEMBL, file=file)

# Cf. PaccMann
dir = "../PaccMann/data/gene_expression"
file = sprintf("%s/gdsc-rnaseq_gene-expression.csv", dir)
cells_paccmann = fread_def(file) %>% rownames
Pred_ChEMBL[Model=="GCNPath" & Seed==2021 & Cell_Line_Name %in% cells_paccmann, ] %>% nrow   # 60539


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
  
  process_utest = function(Perf_UTest, score="RMSE") {
    
    score_ = sprintf("_%s", score)
    col_rm = c("RMSE", "PCC", "SCC") %>% setdiff(score)
    
    Perf_UTest_ = Perf_UTest %>% 
      select(!contains(col_rm)) %>% 
      rename_with(~gsub(score_, "", .x)) %>% 
      rename(Minus_Log10_FDR=MLog10_FDR) %>% 
      relocate(Minus_Log10_FDR, .after=FDR) %>% as.data.frame
    
    return(Perf_UTest_)
  }
  
  
  ### [Source Data] Fig. 4
  Perf_ChEMBL_ = Perf_ChEMBL
  Perf_ChEMBL_ = Perf_ChEMBL_ %>% 
    rename(Train_Seed=Seed, Num_Test=N_Test)
  
  Perf_ChEMBL_Total_ = Perf_ChEMBL_Total
  Perf_ChEMBL_Total_ = Perf_ChEMBL_Total_ %>% 
    rename(Train_Seed=Seed, Train_DB_Cell=Train_DB, Test_DB_Cell=Test_DB, Num_Test=N_Test)
  
  Perf_ChEMBL_List_ = list(Performance=Perf_ChEMBL_, 
                           Performance_RNA=Perf_ChEMBL_Total_)
  
  num_fig = c("a-c", "d-f")
  Perf_ChEMBL_List_ %>% save_for_nc(num=4, suppl=F, num_fig=num_fig)
  
  
  ### [Source Data] Supplementary Fig. 24
  # [Warning] Too large volumn (2.22GB, 23,896,950 rows)
  # Write the file in CSV format with data.table::fwrite (without sheet names)
  col = c("Assay_ChEMBL_ID", "Cell_ChEMBL_ID", "SANGER_MODEL_ID", 
          "Molecule_ChEMBL_ID", "Drug_CID", "LN_IC50", "Prediction", "Seed", "Model")
  
  Pred_ChEMBL_ = Pred_ChEMBL[, col, with=F]
  file = "Prediction [ChEMBL].csv"
  fwrite(Pred_ChEMBL_, file=file, row.names=F)
  
  
  ### [Source Data] Supplementary Fig. 25
  UTest_ChEMBL_ = list(Perf_UTest, Perf_UTest_G, Perf_UTest_C, Perf_UTest_S)
  UTest_ChEMBL_RMSE_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="RMSE"))
  UTest_ChEMBL_RMSE_ %>% save_for_nc(num=25, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 26
  UTest_ChEMBL_PCC_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="PCC"))
  UTest_ChEMBL_PCC_ %>% save_for_nc(num=26, suppl=T)
  
  
  ### [Source Data] Supplementary Fig. 27
  UTest_ChEMBL_SCC_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="SCC"))
  UTest_ChEMBL_SCC_ %>% save_for_nc(num=27, suppl=T)
}

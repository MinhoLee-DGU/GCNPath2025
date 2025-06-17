#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

source("../functions.R")
source("functions.R")
loadings()


# cf. Gather cell list trained in all models

# [GDSC] RF
file = "../RF/_data/SANGER_RNA_TPM.csv"
cells_rf = fread_def(file, col_numeric=T) %>% rownames   # 1432

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

cells_cosmic = Reduce(intersect, list(cells_graphdrp, cells_hidra, cells_rf))   # 962
cells_broad = Reduce(intersect, list(cells_drpreter, cells_tgsa))               # 1337
cells_sanger = Reduce(intersect, list(cells_gcnpath, cells_tgsa_sg))            # 1183

cells_train = IC50_GDSC %>% 
  subset(COSMIC_ID %in% cells_cosmic) %>% 
  subset(BROAD_ID %in% cells_broad) %>% 
  subset(SANGER_MODEL_ID %in% cells_sanger) %>% 
  pull(SANGER_MODEL_ID) %>% unique   # 670



##### 2. Read prediction files

pattern = "pred_chembl_seed([0-9]+).csv"
pattern_tcnns = "pred_chembl_seed([0-9]+)_([0-9]+).csv"

dir = "../RF/Results/IC50_GDSC/Normal"
Pred_RF = read_pred(dir, pattern=pattern, model_name="RF")   # 2375800 [Timing is already considered...]

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
Pred_ChEMBL = Reduce(rbind, list(Pred_RF, Pred_tCNNS[, col, with=F], Pred_HiDRA,
                                 Pred_PaccMann_SG, Pred_GraphDRP[, col, with=F],
                                 Pred_TGDRP, Pred_TGSA, Pred_TGDRP_SG, Pred_TGSA_SG,
                                 Pred_DRPreter, Pred_DRPreter_SG, Pred_GCNPath))   # 32231690

model_list = c("RF", "tCNNS", "HiDRA", "PaccMann_SG", "GraphDRP", 
               "TGDRP", "TGSA", "TGDRP_SG", "TGSA_SG",
               "DRPreter", "DRPreter_SG", "GCNPath")

Pred_ChEMBL$Model = Pred_ChEMBL$Model %>% factor(levels=model_list)
# pred_names = sprintf("Pred_%s", model_list)
# rm(list=ls()[ls() %in% pred_names])



##### Save cell & drug pairs intersecting in all models

Pred_ChEMBL = Pred_ChEMBL %>% filter_pair(model_list=model_list)
Pred_ChEMBL = Pred_ChEMBL %>% subset(SANGER_MODEL_ID %in% cells_train & !Drug_GDSC)
Pred_Num = Pred_ChEMBL %>% count_num(seed_ex=2030)   # 237646 IC50 for each model


# cf. Duplicated synonyms of SIDM00210 [A2780/A2780S]
# In SANGER Passports, two cell lines were considered as identified ones...
# For simplicity, we chose A2780...

col = c("Cell_ChEMBL_ID", "SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID", "Cell_Line_Name")
Anno_Cells_ChEMBL = Pred_ChEMBL[, col, with=F] %>% distinct   # 404 x 3
Anno_Cells_ChEMBL$Cell_ChEMBL_ID %>% unique %>% length        # 405
Anno_Cells_ChEMBL$SANGER_MODEL_ID %>% unique %>% length       # 404
Anno_Cells_ChEMBL %>% group_by(SANGER_MODEL_ID) %>% filter(n()>1)
# CHEMBL3308421  SIDM00210       ACH-000657    906804 A2780         
# CHEMBL4295440  SIDM00210       ACH-000657    906804 A2780S

# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/Anno_Cells.csv", dir)
# Anno_Cells = read.csv(file)
# 
# col = c("SANGER_MODEL_ID", "MODEL_NAME", "SYNONYMS")
# Anno_Cells[, col] %>% subset(MODEL_NAME=="A2780")
# # SANGER_MODEL_ID MODEL_NAME SYNONYMS
# #       SIDM00210      A2780   A2780S
# 
# dir = "../../raw_data/SANGER_Passports"
# file = sprintf("%s/model_list_20230517.csv", dir)
# Anno_Cells_SANGER = read.csv(file)
# 
# col = c("model_id", "model_name", "synonyms")
# Anno_Cells_SANGER[, col] %>% subset(model_name=="A2780")
# #  model_id model_name synonyms
# # SIDM00210      A2780   A2780S

# For simplicity, we chose A2780...
# IC50_GDSC %>% subset(grepl("A2780", CELL_LINE_NAME, ignore.case=T)) %>% pull(CELL_LINE_NAME) %>% unique   # A2780
Pred_ChEMBL = Pred_ChEMBL %>% subset(Cell_ChEMBL_ID!="CHEMBL4295440")   # 26141060 > 26133800
Pred_Num = Pred_ChEMBL %>% count_num(seed_ex=2030)   # 237580 IC50 for each model



##### 3. Performance [Regression]

Perf_ChEMBL = Pred_ChEMBL %>% calc_perf(sum_seed=F)
col = c("RMSE", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Avg = Perf_ChEMBL %>% group_by(Model) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame


Perf_ChEMBL_Best = Perf_ChEMBL %>% group_by(Model) %>% filter(RMSE==min(RMSE)) %>% as.data.frame
Pred_ChEMBL_Best = Pred_ChEMBL %>% select_best_seed(Perf_Best=Perf_ChEMBL_Best)   # 2850960

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
# Total : 2375800
# Infinite (+) : 221708
# Infinite (-) : 162319
# Infinite (Total) : 384027

Pred_ChEMBL_Best %>% subset(Model=="tCNNS") %>% pull(Prediction) %>% is_inf_def
# Total : 237580
# Infinite (+) : 29231
# Infinite (-) : 2
# Infinite (Total) : 29233



##### 4. Read prediction files [Batch Effect of RNA Data]

pattern_gdsc = "pred_chembl_gdsc_seed([0-9]+).csv"
pattern_ccle = "pred_chembl_ccle_seed([0-9]+).csv"
pattern_sanger = "pred_chembl_sanger_seed([0-9]+).csv"

dir = "../RF/Results/IC50_GDSC/Normal"
Pred_RF_C = read_pred(dir, pattern=pattern_ccle, model_name="RF")   # 2375800
Pred_RF_G = read_pred(dir, pattern=pattern_gdsc, model_name="RF")   # 2375800

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

Pred_ChEMBL_G = Reduce(rbind, list(Pred_RF_G, Pred_PaccMann_SG_G, Pred_TGDRP_G, Pred_TGSA_G, 
                                   Pred_TGDRP_SG_G, Pred_TGSA_SG_G,
                                   Pred_DRPreter_G, Pred_DRPreter_SG_G, Pred_GCNPath_G))   # 22651490

Pred_ChEMBL_C = Reduce(rbind, list(Pred_RF_C, Pred_HiDRA_C, Pred_PaccMann_SG_C, 
                                   Pred_TGDRP_SG_C, Pred_TGSA_SG_C,
                                   Pred_DRPreter_SG_C, Pred_GCNPath_C))   # 18442120

Pred_ChEMBL_S = Reduce(rbind, list(Pred_HiDRA_S, Pred_TGDRP_S, Pred_TGSA_S, Pred_DRPreter_S))   # 10989110


model_list_exp = c("RF", "HiDRA", "PaccMann_SG", 
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
tag_ref = Pred_ChEMBL %>% with(paste0(Cell_ChEMBL_ID, "@", Molecule_ChEMBL_ID)) %>% unique   # 237580

Pred_ChEMBL_G = Pred_ChEMBL_G %>% filter_pair(model_list=NULL, tag_ref=tag_ref)
Pred_ChEMBL_C = Pred_ChEMBL_C %>% filter_pair(model_list=NULL, tag_ref=tag_ref)
Pred_ChEMBL_S = Pred_ChEMBL_S %>% filter_pair(model_list=NULL, tag_ref=tag_ref)

Pred_ChEMBL_Ext = Reduce(rbind, list(Pred_ChEMBL_G, Pred_ChEMBL_C, Pred_ChEMBL_S))
Pred_Num_Ext = Pred_ChEMBL_Ext %>% count_num_db(seed_ex=2030)

levels = c("GDSC", "CCLE", "SANGER")
Pred_ChEMBL_Ext$Test_DB = Pred_ChEMBL_Ext$Test_DB %>% factor(levels=levels)
Pred_Num_Ext = Pred_Num_Ext %>% 
  mutate(Test_DB=Test_DB %>% factor(levels=levels)) %>% 
  arrange(Model, Test_DB) %>% as.data.frame

# Fortunately, All cell x drug combinations from original test exist in external test
Pred_Num$IC50 %>% unique       # 237580 IC50 for each model
Pred_Num_Ext$IC50 %>% unique   # 237580 IC50 for each model



##### 5. Performance [Regression, Batch Effect of RNA Data]

Perf_ChEMBL_Ext = Pred_ChEMBL_Ext %>% 
  calc_perf(sum_seed=F, test_db=T)   # 180 [9 x 2 x 10]

col = c("RMSE", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Ext_Avg = Perf_ChEMBL_Ext %>% group_by(Model, Test_DB) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame   # 18

Perf_ChEMBL_Ext_Best = Perf_ChEMBL_Ext %>% 
  group_by(Model, Test_DB) %>% filter(RMSE==min(RMSE)) %>% as.data.frame   # 18
Pred_ChEMBL_Ext_Best = Pred_ChEMBL_Ext %>% 
  select_best_seed(Perf_Best=Perf_ChEMBL_Ext_Best, test_db=T)   # 47516000 > 4751600


# Boxplot
model_g = c("HiDRA")
model_c = c("TGDRP", "TGSA", "DRPreter")
model_s = c("RF", "PaccMann_SG", "TGDRP_SG", "TGSA_SG", "DRPreter_SG", "GCNPath")

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


score = c("RMSE", "PCC", "SCC")
dir = mkdir("Performances [Total]")

Perf_ChEMBL %>% plot_perf(score=score[1], dir=dir, save=T, test_db=F,
                          axis_tl=27, axis_tx=16.5, width=20, height=15)
Perf_ChEMBL %>% plot_perf(score=score[2], dir=dir, save=T, test_db=F, 
                          axis_tl=27, axis_tx=16.5, width=20, height=15)
Perf_ChEMBL %>% plot_perf(score=score[3], dir=dir, save=T, test_db=F, 
                          axis_tl=27, axis_tx=16.5, width=20, height=15)

dir = mkdir("Performances [Total, RMSE<5]")
Perf_ChEMBL %>% 
  plot_perf(score=score[1], dir=dir, test_db=F, ymax=5,
            axis_tl=27, axis_tx=16.5, width=20, height=15, save=T)

dir = mkdir("Performances [RNA Batch Effect]")
Perf_ChEMBL_Total %>% plot_perf(score=score[1], dir=dir, width=27, height=20, save=T, test_db=T)
Perf_ChEMBL_Total %>% plot_perf(score=score[2], dir=dir, width=27, height=20, save=T, test_db=T)
Perf_ChEMBL_Total %>% plot_perf(score=score[3], dir=dir, width=27, height=20, save=T, test_db=T)

dir = mkdir("Performances [RNA Batch Effect, RMSE<5]")
Perf_ChEMBL_Total %>% 
  plot_perf(score=score[1], dir=dir, width=27, height=20, ymax=5, save=T, test_db=T)

col = c("RMSE", "PCC", "SCC")
avg_plus_sd = function(x) sprintf("%.3f±%.3f", mean(x), sd(x))
Perf_ChEMBL_Total_Avg = Perf_ChEMBL_Total %>% group_by(Model, Test_DB) %>%
  reframe(N_Test=max(N_Test), across(all_of(col), avg_plus_sd)) %>% distinct %>% as.data.frame

Perf_ChEMBL_Total_Best = Perf_ChEMBL_Total %>% 
  group_by(Model, Test_DB) %>% filter(RMSE==min(RMSE)) %>% 
  arrange(Model, Test_DB) %>% as.data.frame

file = "Performance [Total].xlsx"
sheet = c("Perf", "Perf_Avg", "Perf_Best")
df_list = list(Perf_ChEMBL_Total, Perf_ChEMBL_Total_Avg, Perf_ChEMBL_Total_Best)
write.xlsx(df_list, file=file, rowNames=F, sheetName=sheet)




##### 6. Performance [Regression, Batch Effect of RNA Data]

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

col = c("Cell_ChEMBL_ID", "SANGER_MODEL_ID", "BROAD_ID", "COSMIC_ID")
Anno_Cells_ChEMBL = Pred_ChEMBL[, col, with=F] %>% distinct   # 404 x 3

file = "Cell_Drug_Pair.csv"
fwrite(Cell_Drug_Pair, file=file)

file = "Anno_Cells_ChEMBL.csv"
fwrite(Anno_Cells_ChEMBL, file=file)

file = "Prediction.csv"
fwrite(Pred_ChEMBL, file=file)

# # Cf. PaccMann
# dir = "../PaccMann/data/gene_expression"
# file = sprintf("%s/gdsc-rnaseq_gene-expression.csv", dir)
# cells_paccmann = fread_def(file) %>% rownames
# Pred_ChEMBL[Model=="GCNPath" & Seed==2021 & Cell_Line_Name %in% cells_paccmann, ] %>% nrow   # 60539




# Estimate Inference Time
IC50_ChEMBL_Time = Pred_ChEMBL[Model=="GCNPath" & Seed==2021, 1:16, with=F]

file = "IC50_ChEMBL_Time.txt"
fwrite(IC50_ChEMBL_Time, file=file, sep="\t")

# The conda environment setting of 
# tCNNS, HiDRA were not compatible to NVIDIA RTX3090...
# Those models were trained and tested with CPUs only...
# Inference time were not estimated to these models...

read_time = function(dir, pattern, sep=",", model_name=NULL) {
  Pred = data.frame()
  df_name_list = list.files(path=dir, pattern=pattern, full.names=T)
  # [pattern] "Pred_[^CCLE]", "Pred_CCLE"
  
  if (length(df_name_list)!=0) {
    for (df_name in df_name_list) {
      try({
        Pred_TP = fread(df_name, header=T, sep=sep)
        df_name_ = strsplit(df_name, "/")[[1]] %>% tail(1)
        nth = gsub(pattern, "\\1", df_name_) %>% as.numeric
        
        Pred_TP$Seed = nth
        Pred = Pred %>% rbind(Pred_TP)
      })
    }
    
    sprintf("# Number of Predictions : %s", nrow(Pred)) %>% print
    if (!is.null(model_name)) Pred$Model = model_name
    return(Pred)
  }
}


# Inference Time with GPU
pattern = "log_time_chembl_seed([0-9]+).csv"

dir = "../RF/Results/IC50_GDSC/Normal"
Time_RF = read_time(dir, pattern=pattern, model_name="RF")

dir = "../GraphDRP/Results/IC50_GDSC/Normal"
Time_GraphDRP = read_time(dir, pattern=pattern, model_name="GraphDRP")

dir = "../PaccMann_SANGER/Results/IC50_GDSC/Normal"
Time_PaccMann_SG = read_time(dir, pattern=pattern, model_name="PaccMann_SG")

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGDRP"
Time_TGDRP_SG = read_time(dir, pattern=pattern, model_name="TGDRP_SG")

dir = "../TGSA_SANGER/Results/IC50_GDSC/Normal/TGSA"
Time_TGSA_SG = read_time(dir, pattern=pattern, model_name="TGSA_SG")

dir = "../DRPreter_SANGER/Results/IC50_GDSC/Normal/DRPreter"
Time_DRPreter_SG = read_time(dir, pattern=pattern, model_name="DRPreter_SG")

dir = "../GCNPath/results/IC50_GDSC/Normal/RGCN"
Time_GCNPath = read_time(dir, pattern=pattern, model_name="GCNPath")


# Inference Time with CPU
pattern_tcnns = "log_time_chembl_seed([0-9]+)_([0-9]+).csv"
dir = "../tCNNS/Results/IC50_GDSC/Normal"
Time_tCNNS = read_time(dir, pattern=pattern_tcnns, model_name="tCNNS")
Time_tCNNS = Time_tCNNS %>% group_by(Seed) %>% summarise(Time=sum(Time)) %>% 
  mutate(Step="Test", Model="tCNNS") %>% relocate(Step, Time, Seed, Model)

dir = "../HiDRA/Results/IC50_GDSC/Normal"
Time_HiDRA = read_time(dir, pattern=pattern, model_name="HiDRA")


# Time units from ms into s (milli-second > second)
Time_List = list(Time_GraphDRP, Time_PaccMann_SG, Time_TGDRP_SG, 
                 Time_TGSA_SG, Time_DRPreter_SG, Time_GCNPath)

Time_List = Reduce(rbind, Time_List) %>% as.data.frame
Time_List$Time = Time_List$Time / 1000


# Barplot of Inference Time
margin1 = margin(r=0.25, l=0.25, unit="cm")
margin2 = margin(0.25, 0.25, 0.25, 0.25, unit="cm")

ylim = c(0, max(Time_List$Time)*1.05)
font1 = font(object="ylab", size=20, margin=margin1)
font2 = font(object="axis.text", size=15, margin=margin2, color="grey30")
font = font1 + font2

pl = Time_List %>% 
  ggbarplot(x="Model", y="Time", color="black", fill="steelblue3",
            label=T, lab.vjust=-1.2, lab.nb.digits=2, xlab=F, ylab="Time [s]", 
            ylim=ylim, position=position_dodge(0.64), 
            add=c("point", "mean_se"), add.params=list(alpha=0.5)) + 
  theme(text=element_text(size=6)) + font + rotate_x_text(36)

main = "Inference Time [NVIDIA RTX3090]"
pl = pl %>% ggpar(legend=NULL)
pl %>% save_fig_ggpubr(main=main, width=15, height=12, svg=T)


# Time units from ms into s (milli-second > second)
Time_List_CPU = list(Time_RF, Time_tCNNS, Time_HiDRA)
Time_List_CPU = Reduce(rbind, Time_List_CPU) %>% as.data.frame
Time_List_CPU$Time = Time_List_CPU$Time / 1000
ylim = c(0, max(Time_List_CPU$Time)*1.05)

pl = Time_List_CPU %>% 
  ggbarplot(x="Model", y="Time", color="black", fill="steelblue3",
            label=T, lab.nb.digits=2, xlab=F, ylab="Time [s]", 
            ylim=ylim, position=position_dodge(0.64), lab.vjust=c(-1.2, -4.8, -1.2), 
            add=c("point", "mean_se"), add.params=list(alpha=0.5)) + 
  theme(text=element_text(size=6)) + font + rotate_x_text(36)

main = "Inference Time [Intel(R) Xeon(R) Gold 5220 CPU]"
pl = pl %>% ggpar(legend=NULL)
pl %>% save_fig_ggpubr(main=main, width=8, height=12, svg=T)




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
  Perf_ChEMBL_ = Perf_ChEMBL %>% 
    rename(Train_Seed=Seed, Num_Test=N_Test)
  
  Perf_ChEMBL_Total_ = Perf_ChEMBL_Total %>% 
    rename(Train_Seed=Seed, Train_DB_Cell=Train_DB, 
           Test_DB_Cell=Test_DB, Num_Test=N_Test)
  
  Perf_ChEMBL_List_ = list(Performance=Perf_ChEMBL_, 
                           Performance_RNA=Perf_ChEMBL_Total_)
  
  num_fig = c("a-c", "d-f")
  Perf_ChEMBL_List_ %>% save_for_nc(num=4, suppl=F, num_fig=num_fig)
  
  
  # ### [Source Data] Supplementary Fig. 21
  # # [Warning] Too large volumn (2.x GB, ~26,000,000 rows)
  # # Write the file in CSV format with data.table::fwrite (without sheet names)
  # col = c("Assay_ChEMBL_ID", "Cell_ChEMBL_ID", "SANGER_MODEL_ID", 
  #         "Molecule_ChEMBL_ID", "Drug_CID", "LN_IC50", "Prediction", "Seed", "Model")
  # 
  # Pred_ChEMBL_ = Pred_ChEMBL[, col, with=F]
  # file = "Prediction [ChEMBL].csv"
  # fwrite(Pred_ChEMBL_, file=file, row.names=F)
  
  
  ### [Source Data] Supplementary Fig. 29
  UTest_ChEMBL_ = list(Perf_UTest, Perf_UTest_G, Perf_UTest_C, Perf_UTest_S)
  UTest_ChEMBL_RMSE_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="RMSE"))
  UTest_ChEMBL_RMSE_ %>% save_for_nc(num=29, suppl=T)
  
  ### [Source Data] Supplementary Fig. 30
  UTest_ChEMBL_PCC_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="PCC"))
  UTest_ChEMBL_PCC_ %>% save_for_nc(num=30, suppl=T)
  
  ### [Source Data] Supplementary Fig. 31
  UTest_ChEMBL_SCC_ = UTest_ChEMBL_ %>% lapply(function(df) df %>% process_utest(score="SCC"))
  UTest_ChEMBL_SCC_ %>% save_for_nc(num=31, suppl=T)
  
  ### [Source Data] Supplementary Fig. 33
  Time_ChEMBL_ = rbind(Time_List, Time_List_CPU) %>% 
    subset(select=-Step) %>% 
    rename(Train_Seed=Seed, Time_Second=Time) %>% 
    relocate(Model, Train_Seed, .before=everything())
  Time_ChEMBL_ %>% save_for_nc(num=33, suppl=T)
}

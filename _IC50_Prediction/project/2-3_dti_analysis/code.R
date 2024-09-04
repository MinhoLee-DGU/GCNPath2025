#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(coop))
suppressMessages(library(cogena))
suppressMessages(library(reshape2))

source("../functions.R")
loadings()

dist_norm = function(df, norm=2) {
  df_names = colnames(df)
  df_dist = do.call(rbind, lapply(df_names, function(i) {
    data.frame(Sample1=i, Sample2=df_names, Dist=colMeans((df-df[, i])^norm, na.rm=T))}))

  # df_dist = df_dist[!duplicated(t(apply(df_dist[1:2], 1, sort))), ]
  df_dist$Dist = df_dist$Dist ** (1/norm)
  return(df_dist)
}

plot_dfs_knn = function(Drug_IC50_Info, x, y, knn=5, mode="max", return_val=T, ...) {
  
  knn_ = sprintf("KNN%s", knn)
  non_knn_ = sprintf("Non-KNN%s", knn)
  
  x = deparse(substitute(x))
  y = deparse(substitute(y))
  Drug_IC50_Info$X_Interest = Drug_IC50_Info[[x]]
  Drug_IC50_Info$Y_Interest = Drug_IC50_Info[[y]]
  
  ifelse_def = function(cond, x, y) {
    ifelse(cond, return(x), return(y))
  }
  
  # Prepare calculating KNN
  col = c("Drug1", "Drug2", "X_Interest", "Y_Interest")
  Drug_IC50_Info = Drug_IC50_Info[, col]
  Drug_IC50_Info_Rev = Drug_IC50_Info %>% dplyr::rename(Drug2=Drug1, Drug1=Drug2)
  
  Drug_IC50_Info = Drug_IC50_Info %>% 
    rbind(Drug_IC50_Info_Rev[, c(2, 1, 3, 4)]) %>% 
    subset(Drug1!=Drug2) %>% distinct %>% as.data.frame
  sprintf("Total number of rows : %s", nrow(Drug_IC50_Info)) %>% print
  
  select_def = ifelse_def(mode=="max", tail, head)
  Drug_IC50_Info = Drug_IC50_Info %>% group_by(Drug1) %>% 
    mutate(KNN=ifelse(X_Interest %in% select_def(sort(X_Interest), n=knn), knn_, non_knn_))
  sprintf("# Total %s pairs are in KNN%s...", sum(Drug_IC50_Info$KNN==knn_), knn) %>% print
  
  col = setNames(c("firebrick1", "black"), c(knn_, non_knn_))
  add = list(scale_colour_manual(values=col))
  
  Drug_IC50_Info = Drug_IC50_Info %>% as.data.frame %>% 
    mutate(KNN=factor(KNN, levels=c(knn_, non_knn_)))
  Drug_IC50_Info %>% plot_def(X_Interest, Y_Interest, color=KNN, add=add, legend="Drug-Pair", ...)
  
  if (return_val) {
    Drug_IC50_Info = Drug_IC50_Info %>% dplyr::rename(Drug_CID1=Drug1, Drug_CID2=Drug2)
    return(Drug_IC50_Info)
  }
}

plot_dfs_multi = function(Drug_IC50_Info, x, knn=5, mode="max", xlab="Tanimoto", 
                          dir=NULL, subtitle=NULL, width=20, height=16, alpha=0.2, 
                          axis_tl=25, axis_tx=20, legend_tl=20, legend_tx=18, return_val=T, save=T) {
  
  x = deparse(substitute(x))
  Drug_IC50_Info$X_Interest = Drug_IC50_Info[[x]]
  
  option = c("Dist_IC50", "Corr_IC50", "Sim_IC50")
  main = sprintf("%s/%s & %s", dir, x, option)
  if (!is.null(subtitle)) main = sprintf("%s [%s]", main, subtitle)
  
  ylab = c(bquote(ln(IC["50"])~Distance~"["*NRMSE*"]"), 
           bquote(ln(IC["50"])~Similarity~"["*PCC*"]"), 
           bquote(ln(IC["50"])~Similarity~"["*RBF*"]"))
  
  DFS_NRMSE = Drug_IC50_Info %>% 
    plot_dfs_knn(X_Interest, Dist_IC50, knn=knn, mode=mode, 
                 main=main[1], xlab=xlab, ylab=ylab[[1]], 
                 size=1.5, alpha=alpha, width=width, height=height, 
                 axis_tl=axis_tl, axis_tx=axis_tx, 
                 legend_tl=legend_tl, legend_tx=legend_tx,
                 raster=T, return_val=return_val, save=save, save_svg=save)
  
  DFS_PCC = Drug_IC50_Info %>% 
    plot_dfs_knn(X_Interest, PCC_IC50, knn=knn, mode=mode,
                 main=main[2], xlab=xlab, ylab=ylab[[2]], 
                 size=1.5, alpha=alpha, width=width, height=height,
                 axis_tl=axis_tl, axis_tx=axis_tx, 
                 legend_tl=legend_tl, legend_tx=legend_tx,
                 raster=T, return_val=return_val, save=save, save_svg=save)
  
  # DFS_Sim = Drug_IC50_Info %>% 
  #   plot_dfs_knn(X_Interest, Sim_IC50, knn=knn, mode=mode,
  #                main=main[3], xlab=xlab, ylab=ylab[[3]], 
  #                size=1.5, alpha=alpha, width=width, height=height,
  #                axis_tl=axis_tl, axis_tx=axis_tx, 
  #                legend_tl=legend_tl, legend_tx=legend_tx,
  #                raster=T, return_val=return_val, save=save, save_svg=save)
  
  colnames(DFS_NRMSE)[3:4] = c("Drug_Similarity", "NRMSE_LN_IC50")
  colnames(DFS_PCC)[3:4] = c("Drug_Similarity", "PCC_LN_IC50")
  # colnames(DFS_Sim)[3:4] = c("Drug_Similarity", "Sim_LN_IC50")
  
  # DFS = list(DFS_NRMSE=DFS_NRMSE, DFS_PCC=DFS_PCC, DFS_Sim=DFS_Sim)
  DFS = list(DFS_NRMSE=DFS_NRMSE, DFS_PCC=DFS_PCC)
  if (return_val) return(DFS)
}

# dist_norm = function(x, y, norm=2, except_inf=T) {
#   if (length(x)!=length(y)) {
#     if (except_inf) diff = diff[!is.infinite(diff)]
#     dist = mean(diff, na.rm=T)**(1/norm)
#   }
#   return(dist)
# }


dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/GDSC_Drug_Morgan.csv", dir)
GDSC_Morgan = read.csv(file, row.names=1)
# Just run the code "gdsc_last/project/2-3_dti_analysis/process_drug.py" to create the file

file = sprintf("%s/Anno_Drugs.csv", dir)
Anno_Drugs = read.csv(file)

dir = "../../processed_data/ic50_data/GDSC"
file = sprintf("%s/IC50_GDSC.csv", dir)
IC50_GDSC = read.csv(file)

dir = "../../raw_data/MSigDB"
file = sprintf("%s/c2.cp.biocarta.v2023.1.Hs.symbols.gmt", dir)
Path_List = gmt2list(file)
Path_List %>% unlist %>% unique %>% length   # 1509

dir = "../../processed_data/net_data/STRING"
file = sprintf("%s/STRING_Filt_Sym.csv", dir)
STRING = fread(file)
STRING = STRING %>% as.data.frame

dir = "../../processed_data/drug_data/STiTCH"
file = sprintf("%s/GDSC_DTI_Wide.csv", dir)
GDSC_DTI_Wide = read.csv(file, row.names=1)

genes_string = union(STRING$Node1, STRING$Node2)   # 16812
sum(colnames(GDSC_DTI_Wide) %in% genes_string)     # 2130 [from 2188]


# Calculate Dist, PCC and RBF Similarities between LN_IC50s of drug-pairs
IC50_GDSC_Wide = acast(IC50_GDSC, SANGER_MODEL_ID~DRUG_CID, value.var="LN_IC50")
PCC_Drug_IC50 = IC50_GDSC_Wide %>% cor(method="pearson", use="pairwise.complete.obs") %>% as.data.frame
SCC_Drug_IC50 = IC50_GDSC_Wide %>% cor(method="spearman", use="pairwise.complete.obs") %>% as.data.frame

PCC_Drug_IC50 = PCC_Drug_IC50 %>% as.matrix %>% reshape2::melt() %>% as.data.frame
SCC_Drug_IC50 = SCC_Drug_IC50 %>% as.matrix %>% reshape2::melt() %>% as.data.frame

col1 = c("Drug1", "Drug2", "PCC_IC50")
col2 = c("Drug1", "Drug2", "SCC_IC50")

colnames(PCC_Drug_IC50) = col1
colnames(SCC_Drug_IC50) = col2

# Dist_Drug_IC50 = IC50_GDSC_Wide %>% stat_self(method=dist_norm, by_row=F)
Dist_Drug_IC50 = IC50_GDSC_Wide %>% dist_norm(self_dist=T)
col = c("Drug1", "Drug2", "Dist_IC50")
colnames(Dist_Drug_IC50) = col

PCC_Drug_IC50 = PCC_Drug_IC50 %>% mutate(Drug1=as.character(Drug1), Drug2=as.character(Drug2))
SCC_Drug_IC50 = SCC_Drug_IC50 %>% mutate(Drug1=as.character(Drug1), Drug2=as.character(Drug2))
Dist_Drug_IC50 = Dist_Drug_IC50 %>% mutate(Drug1=as.character(Drug1), Drug2=as.character(Drug2))

by = c("Drug1"="Drug1", "Drug2"="Drug2")
Drug_IC50_Info = inner_join(Dist_Drug_IC50, PCC_Drug_IC50, by=by)
Drug_IC50_Info = inner_join(Drug_IC50_Info, SCC_Drug_IC50, by=by)
# plot(Drug_IC50_Info$PCC_IC50, Drug_IC50_Info$SCC_IC50)   # PCC~SCC Mostly Same...

sigma = Drug_IC50_Info$Dist_IC50 %>% mean   # 3.151784
Drug_IC50_Info = Drug_IC50_Info %>% mutate(Sim_IC50=exp(-Dist_IC50/sigma))
# Drug_IC50_Info$Sim_IC50 %>% hist



##### 2. Drug Functional Similarity

### 2-0. Structural Similarity [Tanimoto, Morgan]

cal_tanimoto = function(ECFP, cores=1, rev=T) {
  
  suppressMessages(library(Rcpi))
  suppressMessages(library(rJava))
  nC2 = combn(nrow(ECFP), 2) %>% t
  
  Tanimoto = data.frame(Drug1 = rownames(ECFP)[nC2[, 1]], 
                        Drug2 = rownames(ECFP)[nC2[, 2]])
  
  if (cores>=2) {
    cluster = multicores(cores)
    on.exit(stopCluster(cluster))
    tanimoto = foreach(i=1:nrow(nC2), .combine=c) %dopar% {
      tc = calcDrugFPSim(ECFP[nC2[i, 1], ], ECFP[nC2[i, 2], ], "complete", "tanimoto")
      return(tc)
    }
  } else {
    tanimoto = c()
    for (i in 1:nrow(nC2)) {
      tc = calcDrugFPSim(ECFP[nC2[i, 1], ], ECFP[nC2[i, 2], ], "complete", "tanimoto")
      tanimoto = tanimoto %>% c(tc)
    }
  }
  
  Tanimoto = Tanimoto %>% cbind(tanimoto)
  colnames(Tanimoto)[3] = "Tanimoto"
  
  if (rev) {
    Tanimoto_Rev = Tanimoto %>% dplyr::rename(Drug2=Drug1, Drug1=Drug2)
    Tanimoto = rbind(Tanimoto, Tanimoto_Rev) %>% distinct(.keep_all=T)
  }
  
  return(Tanimoto)
}

GDSC_Morgan_Tan = cal_tanimoto(GDSC_Morgan, cores=20, rev=T)
# Drug1, Drug2, Tanimoto [3, D(D-1)]

by = c("Drug1"="Drug1", "Drug2"="Drug2")
Drug_IC50_Info = full_join(Drug_IC50_Info, GDSC_Morgan_Tan, by=by)
Drug_IC50_Info$Tanimoto[is.na(Drug_IC50_Info$Tanimoto)] = 1


dir = mkdir("../../processed_data/drug_data/DTI_Analysis/Tanimoto")
xlab = "Structural Similarity [Tanimoto, Morgan]"
GDSC_DFS_Tan_Morgan_ = Drug_IC50_Info %>% 
  plot_dfs_multi(Tanimoto, xlab=xlab, dir=dir, save=T)

file = sprintf("%s/GDSC_Morgan_Tanimoto.csv", dir)
write.csv(GDSC_Morgan_Tan, file=file, row.names=F)



### 2-0. Structural Similarity [Tanimoto, PubChem]

suppressMessages(library(Rcpi))
suppressMessages(library(rcdk))

dir = "../../processed_data/drug_data/GDSC"
file = sprintf("%s/SMILES_GDSC.csv", dir)
SMILES_GDSC = read.csv(file)

GDSC_PubChem = data.frame()
mol_gdsc = SMILES_GDSC$SMILES_CAN %>% parse.smiles

for (i in 1:length(mol_gdsc)) {
  pubchem_fp = extractDrugPubChemComplete(mol_gdsc[[i]])
  GDSC_PubChem = GDSC_PubChem %>% rbind(pubchem_fp)
}

rownames(GDSC_PubChem) = SMILES_GDSC$Drug_CID
colnames(GDSC_PubChem) = sprintf("FP%s", 0:(ncol(GDSC_PubChem)-1))
# GDSC_PubChem %>% rowSums %>% hist   # 400 x 881

GDSC_PubChem_Tan = cal_tanimoto(GDSC_PubChem, cores=20)
colnames(GDSC_PubChem_Tan)[3] = "Tanimoto_PC"
# Drug1, Drug2, Tanimoto [3, D(D-1)]

by = c("Drug1"="Drug1", "Drug2"="Drug2")
Drug_IC50_Info = full_join(Drug_IC50_Info, GDSC_PubChem_Tan, by=by)
Drug_IC50_Info$Tanimoto_PC[is.na(Drug_IC50_Info$Tanimoto_PC)] = 1


dir = "../../processed_data/drug_data/DTI_Analysis/Tanimoto"
xlab = "Structural Similarity [Tanimoto, PubChem]"
GDSC_DFS_Tan_PubChem_ = Drug_IC50_Info %>% 
  plot_dfs_multi(Tanimoto_PC, xlab=xlab, dir=dir, axis_tl=22.5, save=T)

xlab = "Structural Similarity [Tanimoto, Morgan]"
ylab = "Structural Similarity [Tanimoto, PubChem]"
main = sprintf("%s/Tanimoto & Tanimoto_PC", dir)
Drug_IC50_Info %>% plot_def(Tanimoto, Tanimoto_PC, main=main[1], 
                            xlab=xlab, ylab=ylab, width=20, height=20,
                            axis_tl=20, axis_tx=18, xy_line=T, save=T)

file = sprintf("%s/GDSC_PubChem_Tanimoto.csv", dir)
write.csv(GDSC_PubChem_Tan, file=file, row.names=F)



### 2-1. RWR + ssGSEA or singscore

fill_row_na = function(df, row, value=0) {
  
  df = df %>% as.data.frame
  row_na = row[!(row %in% rownames(df))]
  
  df_row_na = rep(value, length(row_na)*ncol(df)) %>% 
    matrix(nrow=length(row_na), ncol=ncol(df)) %>% as.data.frame
  rownames(df_row_na) = row_na
  colnames(df_row_na) = colnames(df)
  
  df = df %>% rbind(df_row_na)
  return(df)
}

fill_col_na = function(df, col, value=0) {
  
  df = df %>% as.data.frame
  col_na = col[!(col %in% colnames(df))]
  
  df_col_na = rep(value, nrow(df)*length(col_na)) %>% 
    matrix(nrow=nrow(df), ncol=length(col_na)) %>% as.data.frame
  rownames(df_col_na) = rownames(df)
  colnames(df_col_na) = col_na
  
  df = df %>% cbind(df_col_na)
  return(df)
}

gsva_def = function(Omics, Path_List, method="gsva", 
                    filt_genes=T, do_qnorm=F, cores=10, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  suppressMessages(library(stats))
  gene_list = Path_List %>% unlist %>% unique
  
  n_before = length(gene_list)
  n_after = sum(colnames(Omics) %in% gene_list)
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    n_before = ncol(Omics)
    Omics = Omics[, colnames(Omics) %in% gene_list]
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  Omics = Omics %>% t %>% as.matrix
  
  if (do_qnorm) {
    suppressMessages(library(preprocessCore))
    Omics = Omics %>% normalize.quantiles(keep.names=T)
  }
  
  Omics_Path = gsva(Omics, Path_List, method=method, parallel.sz=cores, ...)
  Omics_Path = Omics_Path %>% t %>% as.data.frame
  return(Omics_Path)
}

singscore_def = function(Omics, Path_List, path_names=NULL, filt_genes=T, 
                         direction=T, center=T, stableGenes=NULL, ...) {
  
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  
  suppressMessages(library(GSEABase))
  suppressMessages(library(singscore))
  gene_list = Path_List %>% unlist %>% unique
  
  n_before = length(gene_list)
  n_after = sum(colnames(Omics) %in% gene_list)
  int_ratio = round(100 * n_after/n_before, 2)
  sprintf("# The omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before) %>% print
  
  if (filt_genes) {
    n_before = ncol(Omics)
    Omics = Omics[, colnames(Omics) %in% gene_list]
    sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after) %>% print
  } else print("# Genes are not filtered from omics data...")
  
  if (is.null(path_names)) path_names = names(Path_List)
  Path_List = Path_List %>% lapply(GeneSet)
  Omics = Omics %>% t %>% rankGenes(stableGenes=stableGenes)
  for (i in 1:length(path_names)) Path_List[[i]]@setName = path_names[[i]]
  
  Omics_Path = multiScore(Omics, upSetColc=Path_List, knownDirection=direction, centerScore=center, ...)
  Omics_Path = Omics_Path[[1]] %>% t %>% as.data.frame
  return(Omics_Path)
}

rwr_path = function(Omics, Net, Path_List, col_node1="Node1", col_node2="Node2", 
                    mode="singscore", normalise="laplacian", normalise_res="quantile", restart=0.7, cores=15) {
  
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  
  suppressMessages(library(dnet))
  suppressMessages(library(igraph))
  # To install dnet, 2 Bioconductor packages must be installed first (supraHex, Rgraphviz)
  
  genes_net = union(Net[, col_node1], Net[, col_node2])
  genes_out = setdiff(colnames(Omics), genes_net)
  sprintf("# Genes not in Network were filtered out... [%s]", length(genes_out)) %>% print
  
  Omics = Omics[, colnames(Omics) %in% genes_net]
  Omics = Omics %>% fill_col_na(genes_net)
  Omics = Omics[, genes_net] %>% t %>% as.data.frame
  
  genes = rownames(Omics)
  samples = colnames(Omics)
  Net = Net[, c(col_node1, col_node2)]
  colnames(Net) = c("from", "to")
  
  Net = graph_from_data_frame(Net, vertices=genes, directed=F)
  Omics = dRWR(Net, setSeeds=Omics, restart=restart, verbose=F, 
               normalise=normalise, normalise.affinity.matrix=normalise_res)
  
  norm = function(x) x / sum(x)
  Omics = Omics %>% as.matrix %>% apply(2, norm) %>% t %>% as.data.frame

  rownames(Omics) = samples
  colnames(Omics) = genes
  sprintf("# %s-bit fingerprints are processed...", length(genes)) %>% print

  if (mode %in% c("ssgsea", "gsva")) {
    Omics_Path = Omics %>% gsva_def(method=mode, Path_List=Path_List, cores=cores)
  } else {
    Omics_Path = Omics %>% singscore_def(Path_List=Path_List)
  }

  return(Omics_Path)
}

analyze_dpi = function(DPI, Drug_IC50_Info, subtitle="RWR_DTI", method="pearson", dir=".") {
  
  DPI_Corr = DPI %>% t %>% cor(method=method) %>% as.data.frame
  DPI_Corr = DPI_Corr %>% as.matrix %>% reshape2::melt() %>% as.data.frame
  
  colnames(DPI_Corr) = c("Drug1", "Drug2", "Corr_DPI")
  DPI_Corr = DPI_Corr %>% mutate(Drug1=as.character(Drug1), Drug2=as.character(Drug2))
  
  by = c("Drug1"="Drug1", "Drug2"="Drug2")
  col = c("Drug1", "Drug2", "Dist_IC50", "PCC_IC50", "SCC_IC50", "Sim_IC50")
  DFS = inner_join(DPI_Corr, Drug_IC50_Info[, col], by=by)
  
  xlab = "Target Similarity [PCC]"
  DFS %>% plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle, save=T)
  return(DFS)
}

analyze_dpi_rwr = function(DTI, Net, Path_List, Drug_IC50_Info, 
                        mode="ssgsea", restart=0.7, cores=15, dir=".") {
  
  mode_ = list(ssgsea="ssGSEA", singscore="SING")
  mode_ = mode_[[mode]]
  
  subtitle = sprintf("RWR+%s, restart=%s", mode_, restart)
  DPI_RWR = DTI %>% rwr_path(Net, Path_List, mode=mode, restart=restart, cores=cores)
  DFS_RWR = analyze_dpi(DPI_RWR, Drug_IC50_Info, subtitle=subtitle, dir=dir)
  
  file = sprintf("%s/GDSC_DPI_%s_%s.csv", dir, mode_, restart)
  write.csv(DPI_RWR, file=file, row.names=T)
  
  DFS = list(DPI=DPI_RWR, DFS=DFS_RWR)
  return(DFS)
}

restart = c(0.3, 0.5, 0.7, 0.9)
dir = mkdir("../../processed_data/drug_data/DTI_Analysis/RWR_ssGSEA")
# 16,812-bit fingerprints are processed... [16,812 STRING genes]
# Omics data contain 98.28% pathway genes... [1483 / 1509 Pathway genes]

GDSC_DFS_ssGSEA_0.3 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="ssgsea", restart=restart[1], dir=dir)
GDSC_DFS_ssGSEA_0.5 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="ssgsea", restart=restart[2], dir=dir)
GDSC_DFS_ssGSEA_0.7 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="ssgsea", restart=restart[3], dir=dir)
GDSC_DFS_ssGSEA_0.9 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="ssgsea", restart=restart[4], dir=dir)

dir = mkdir("../../processed_data/drug_data/DTI_Analysis/RWR_SING")
# 16,812-bit fingerprints are processed... [16,812 STRING genes]
# Omics data contain 98.28% pathway genes... [1483 / 1509 Pathway genes]

GDSC_DFS_SING_0.3 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="singscore", restart=restart[1], dir=dir)
GDSC_DFS_SING_0.5 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="singscore", restart=restart[2], dir=dir)
GDSC_DFS_SING_0.7 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="singscore", restart=restart[3], dir=dir)
GDSC_DFS_SING_0.9 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwr(STRING, Path_List, Drug_IC50_Info, mode="singscore", restart=restart[4], dir=dir)

plot_only = T
restart_list = c(0.3, 0.5, 0.7, 0.9)
xlab = "Target Similarity [PCC]"

if (plot_only) {
  dir = "../../processed_data/drug_data/DTI_Analysis/RWR_ssGSEA"
  subtitle = sprintf("RWR+ssGSEA, restart=%s", restart_list)
  GDSC_DFS_ssGSEA_0.3_ = GDSC_DFS_ssGSEA_0.3$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[1], save=T)
  GDSC_DFS_ssGSEA_0.5_ = GDSC_DFS_ssGSEA_0.5$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[2], save=T)
  GDSC_DFS_ssGSEA_0.7_ = GDSC_DFS_ssGSEA_0.7$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[3], save=T)
  GDSC_DFS_ssGSEA_0.9_ = GDSC_DFS_ssGSEA_0.9$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[4], save=T)
  
  dir = "../../processed_data/drug_data/DTI_Analysis/RWR_SING"
  subtitle = sprintf("RWR+singscore, restart=%s", restart_list)
  GDSC_DFS_SING_0.3_ = GDSC_DFS_SING_0.3$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[1], save=T)
  GDSC_DFS_SING_0.5_ = GDSC_DFS_SING_0.5$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[2], save=T)
  GDSC_DFS_SING_0.7_ = GDSC_DFS_SING_0.7$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[3], save=T)
  GDSC_DFS_SING_0.9_ = GDSC_DFS_SING_0.9$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[4], save=T)
}



### 2-2. RWR similarity using Target-Pathway hegerogenous network [Drug Target 1 & Non-Target 0]
# [2020] MAP, Predicting miRNA-based disease-disease relationships through network diffusion on multi-omics biological data
# doi.org/10.1038/s41598-020-65633-6

rwr_path_hetero = function(Omics, Net, Path_List, col_node1="Node1", col_node2="Node2", 
                           normalise="laplacian", normalise_res="quantile", restart=0.7, cores=15) {
  
  suppressMessages(library(dnet))
  suppressMessages(library(igraph))
  # Omics [C x G] > Path [C x Path]
  # To install dnet, 2 Bioconductor packages must be installed first (supraHex, Rgraphviz)
  
  genes_net = union(Net[, col_node1], Net[, col_node2])
  genes_out = setdiff(colnames(Omics), genes_net)
  sprintf("# Genes not in Network were filtered out... [%s]", length(genes_out)) %>% print
  
  Omics = Omics[, colnames(Omics) %in% genes_net]
  Omics = Omics %>% fill_col_na(genes_net)
  Omics = Omics %>% t %>% as.data.frame
  
  path_genes = Path_List %>% unlist %>% as.character
  path_names = rep(names(Path_List), sapply(Path_List, length))
  Net_Path = data.frame(from=path_genes, to=path_names) %>% subset(from %in% genes_net)
  
  Net = Net[, c(col_node1, col_node2)]
  colnames(Net) = c("from", "to")
  Net = Net %>% rbind(Net_Path)
  Omics = Omics %>% fill_row_na(names(Path_List)) %>% as.data.frame
  
  genes = rownames(Omics)
  samples = colnames(Omics)
  Net = graph_from_data_frame(Net, vertices=genes, directed=F)
  
  Omics = dRWR(Net, setSeeds=Omics, restart=restart, verbose=F,
               normalise=normalise, normalise.affinity.matrix=normalise_res)
  
  Omics = Omics %>% as.matrix %>% t %>% as.data.frame
  rownames(Omics) = samples
  colnames(Omics) = genes
  sprintf("# %s-bit fingerprints are processed...", length(genes)) %>% print
  
  return(Omics)
}

analyze_dpi_rwrh = function(DTI, Net, Path_List, Drug_IC50_Info, 
                        restart=0.7, cores=15, dir=".") {
  
  subtitle = sprintf("RWR_Hetero, restart=%s", restart)
  DPI_RWR = DTI %>% rwr_path_hetero(Net, Path_List, restart=restart, cores=cores)
  DFS_RWR = analyze_dpi(DPI_RWR, Drug_IC50_Info, subtitle=subtitle, dir=dir)
  
  file = sprintf("%s/GDSC_DPI_Hetero_%s.csv", dir, restart)
  write.csv(DPI_RWR, file=file, row.names=T)
  
  DFS = list(DPI=DPI_RWR, DFS=DFS_RWR)
  return(DFS)
}

restart = c(0.3, 0.5, 0.7, 0.9)
dir = mkdir("../../processed_data/drug_data/DTI_Analysis/RWR_Path_Hetero")
# 17104-bit fingerprints are processed...  [16,812 STRING genes + 292 Pathway]

GDSC_DFS_RWRH_0.3 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwrh(STRING, Path_List, Drug_IC50_Info, restart=restart[1], dir=dir)
GDSC_DFS_RWRH_0.5 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwrh(STRING, Path_List, Drug_IC50_Info, restart=restart[2], dir=dir)
GDSC_DFS_RWRH_0.7 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwrh(STRING, Path_List, Drug_IC50_Info, restart=restart[3], dir=dir)
GDSC_DFS_RWRH_0.9 = GDSC_DTI_Wide %>% 
  analyze_dpi_rwrh(STRING, Path_List, Drug_IC50_Info, restart=restart[4], dir=dir)

plot_only = T
restart_list = c(0.3, 0.5, 0.7, 0.9)
xlab = "Target Similarity [PCC]"

if (plot_only) {
  dir = "../../processed_data/drug_data/DTI_Analysis/RWR_Path_Hetero"
  subtitle = sprintf("RWR_Hetero, restart=%s", restart_list)
  GDSC_DFS_RWRH_0.3_ = GDSC_DFS_RWRH_0.3$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[1], save=T)
  GDSC_DFS_RWRH_0.5_ = GDSC_DFS_RWRH_0.5$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[2], save=T)
  GDSC_DFS_RWRH_0.7_ = GDSC_DFS_RWRH_0.7$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[3], save=T)
  GDSC_DFS_RWRH_0.9_ = GDSC_DFS_RWRH_0.9$DFS %>% 
    plot_dfs_multi(Corr_DPI, xlab=xlab, dir=dir, subtitle=subtitle[4], save=T)
}



### 2-3. Pathway similarity using Target-Pathway hegerogenous network [PathSim]
# [2021] drugSim, Investigation of pharmacological mechanism of natural product using pathway fingerprints similarity based on drug-target-pathway heterogenous network
# doi.org/10.1186/s13321-021-00549-5
# https://github.com/huihui1126/drugSim-pathway

into_adj = function(adj_list) {
  # Adjecancy List to Adjecancy Matrix
  v_from = adj_list %>% unlist %>% as.character
  v_to = rep(names(adj_list), sapply(adj_list, length))
  Links = data.frame(from = v_from, to = v_to)
  Links = Links %>% acast(from~to, value.var="from", fun.aggregate=length)
  return(Links)
}

mult_adj = function(adj1, adj2) {
  # Multiply 2 Matrix with intersecting y_adj1 & x_adj2
  rows = rownames(adj1)
  cols = colnames(adj2)
  v = intersect(colnames(adj1), rownames(adj2))
  adj = as.matrix(adj1[, v]) %*% as.matrix(adj2[v, ]) %>% as.data.frame
  
  rownames(adj) = rows
  colnames(adj) = cols
  return(adj)
}

mult_adj_round = function(...) {
  # Multiply Multiple Matrix
  args = list(...)
  adj = Reduce(mult_adj, args)
  adj = mult_adj(adj, t(adj))
  return(adj)
}

calc_pathsim = function(DTI_Adj=NULL, TPI_Adj=NULL, PPI_Adj=NULL, return_adj=T, na_zero=T) {
  
  if (is.null(PPI_Adj)) {
    PathSim_Adj = mult_adj_round(DTI_Adj, TPI_Adj)
    # [D x D] Drug > Target > Pathway > Target > Drug
  } else {
    PathSim_Adj = mult_adj_round(DTI_Adj, PPI_Adj, TPI_Adj)
    # [D x D] Drug > Target > Target > Pathway > Target > Target > Drug
  }
  
  PathSim = PathSim_Adj
  diag_adj = PathSim_Adj %>% as.matrix %>% diag %>% as.numeric
  
  for (i in 1:nrow(PathSim_Adj)) {
    for (j in 1:ncol(PathSim_Adj)) {
      sym = 2*PathSim_Adj[i, j] / (diag_adj[i] + diag_adj[j])
      # Normalize by the sum of the number of metapath to themselves
      sym = ifelse(na_zero & is.na(sym), 0, sym)
      PathSim[i, j] = sym
    }
  }
  
  if (!return_adj) {
    return(PathSim)
  } else {
    PathSim = list(Adj=PathSim_Adj, DFS=PathSim)
    return(PathSim)
  }
}

analyze_dfs = function(DFS, Drug_IC50_Info, subtitle="PathSim", dir=".", neighbor=0, return_knn=F) {
  
  DFS = DFS %>% as.matrix %>% reshape2::melt() %>% as.data.frame
  colnames(DFS) = c("Drug1", "Drug2", "PathSim")
  DFS = DFS %>% mutate(Drug1=as.character(Drug1), Drug2=as.character(Drug2))
  
  by = c("Drug1"="Drug1", "Drug2"="Drug2")
  col = c("Drug1", "Drug2", "Dist_IC50", "PCC_IC50", "SCC_IC50", "Sim_IC50")
  DFS = inner_join(DFS, Drug_IC50_Info[, col], by=by)
  
  xlab = sprintf("Target Similarity [PathSim, %s-hop]", neighbor)
  DFS_KNN = DFS %>% plot_dfs_multi(PathSim, xlab=xlab, dir=dir, subtitle=subtitle, return_val=T, save=T)
  
  if (return_knn) {
    return(DFS_KNN)
  } else return(DFS)
}

Path_Adj = Path_List %>% into_adj   # 1509 x 292
STRING_Adj = STRING %>% acast(Node1~Node2, value.var="Node2", fun.aggregate=length)   # 16812 x 16812
# STRING_Adj2 = STRING_Adj %*% STRING_Adj   # Too slow... not recommended

PathSim_0hop = calc_pathsim(DTI_Adj=GDSC_DTI_Wide, TPI_Adj=Path_Adj)
PathSim_1hop = calc_pathsim(DTI_Adj=GDSC_DTI_Wide, TPI_Adj=Path_Adj, PPI_Adj=STRING_Adj)
# PathSim_2hop = calc_pathsym(DTI_Adj=GDSC_DTI_Wide, TPI_Adj=Path_Adj, PPI_Adj=STRING_Adj2)


subtitle = sprintf("PathSim (%s-Hop)", 0:2)
dir = mkdir("../../processed_data/drug_data/DTI_Analysis/PathSim")
GDSC_DFS_PathSim0 = analyze_dfs(PathSim_0hop$DFS, Drug_IC50_Info, dir=dir, subtitle=subtitle[1], neighbor=0)
GDSC_DFS_PathSim1 = analyze_dfs(PathSim_1hop$DFS, Drug_IC50_Info, dir=dir, subtitle=subtitle[2], neighbor=1)
# GDSC_DFS_PathSim2 = analyze_dfs(PathSim_2hop$DFS, Drug_IC50_Info, dir=dir, subtitle=subtitle[3], neighbor=2)

GDSC_DFS_PathSim0_ = analyze_dfs(PathSim_0hop$DFS, Drug_IC50_Info, dir=dir, 
                                 subtitle=subtitle[1], neighbor=0, return_knn=T)
GDSC_DFS_PathSim1_ = analyze_dfs(PathSim_1hop$DFS, Drug_IC50_Info, dir=dir, 
                                 subtitle=subtitle[2], neighbor=1, return_knn=T)

GDSC_DFS_PathSim0 = GDSC_DFS_PathSim0 %>% 
  subset(select=c(Drug1, Drug2, PathSim))
GDSC_DFS_PathSim1 = GDSC_DFS_PathSim1 %>% 
  subset(select=c(Drug1, Drug2, PathSim))

file = sprintf("%s/PathSim_Adj_%shop.csv", dir, 0:1)
write.csv(PathSim_0hop$Adj, file=file[1])
write.csv(PathSim_1hop$Adj, file=file[2])



### 2-4. Target GO semantic similarity

# run "sbatch target_go_sim.sh"



### 2-5. Target Seq Similarity
# Calculate normalized sequence similarities between targets
# Smith-Waterman [Local] & Needleman-Wunsch [Global]
# [2011] SITAR, Combining Drug and Gene Similarity Measures for Drug-Target Elucidation
# doi.org/10.1089/cmb.2010.0213

# run "sbatch target_seq_sim.sh"



### Visualize target similarities calculated in 2-4 & 2-5

into_pairwise = function(DFS) {
  DFS_Rev = DFS
  colnames(DFS_Rev)[1:2] = c("Drug2", "Drug1")
  DFS = DFS %>% rbind(DFS_Rev) %>% distinct
  sprintf("# MoA Similarity : %s", nrow(DFS)) %>% print
  return(DFS)
}

dir = "../../processed_data/drug_data/DTI_Analysis"

file = sprintf("%s/Target_GO/MoA_Similarity.csv", dir)
GDSC_DFS_GO = read.csv(file)
GDSC_DFS_GO = GDSC_DFS_GO %>% into_pairwise

GDSC_DFS_GO = GDSC_DFS_GO %>% 
  mutate(Drug1=as.character(Drug1), 
         Drug2=as.character(Drug2))

file = sprintf("%s/Target_Seq/MoA_Similarity.csv", dir)
GDSC_DFS_Seq = read.csv(file)
GDSC_DFS_Seq = GDSC_DFS_Seq %>% into_pairwise

GDSC_DFS_Seq = GDSC_DFS_Seq %>% 
  mutate(Drug1=as.character(Drug1), 
         Drug2=as.character(Drug2))

by = c("Drug1"="Drug1", "Drug2"="Drug2")
merge_ = function(df1, df2) full_join(df1, df2, by=by)
GDSC_DFS_GO_Seq = Reduce(merge_, list(GDSC_DFS_GO, GDSC_DFS_Seq))

col = c("Drug1", "Drug2", "Dist_IC50", "PCC_IC50", "SCC_IC50", "Sim_IC50")
GDSC_DFS_GO_Seq = inner_join(GDSC_DFS_GO_Seq, Drug_IC50_Info[, col], by=by)


option_go = c("GO:BP", "GO:MF", "GO:CC")
xlab = sprintf("Target Similarity [%s Enrichment]", option_go)
dir = "../../processed_data/drug_data/DTI_Analysis/Target_GO"

GDSC_DFS_GOBP_Seq_ = GDSC_DFS_GO_Seq %>% 
  plot_dfs_multi(Sim_GOBP, xlab=xlab[1], dir=dir, save=T)
GDSC_DFS_GOMF_Seq_ = GDSC_DFS_GO_Seq %>% 
  plot_dfs_multi(Sim_GOMF, xlab=xlab[2], dir=dir, save=T)
GDSC_DFS_GOCC_Seq_ = GDSC_DFS_GO_Seq %>% 
  plot_dfs_multi(Sim_GOCC, xlab=xlab[3], dir=dir, save=T)

option_seq = c("Local Alignment", "Global Alignment")
xlab = sprintf("Target Similarity [%s]", option_seq)
dir = "../../processed_data/drug_data/DTI_Analysis/Target_Seq"

GDSC_DFS_Local_Seq_ = GDSC_DFS_GO_Seq %>% 
  plot_dfs_multi(Sim_Seq_Local, xlab=xlab[1], dir=dir, save=T)
GDSC_DFS_Global_Seq_ = GDSC_DFS_GO_Seq %>% 
  plot_dfs_multi(Sim_Seq_Global, xlab=xlab[2], dir=dir, save=T)

GDSC_DFS_GO_Seq = GDSC_DFS_GO_Seq %>% 
  subset(select=c(Drug1, Drug2, Sim_GOBP, Sim_GOMF, Sim_GOCC, Sim_Seq_Local, Sim_Seq_Global))



### 2-6. Target pathway annotated in GDSC

Anno_Drugs$Target_Pathway %>% is.na %>% sum   # 3
Anno_Drugs$Target_Pathway[is.na(Anno_Drugs$Target_Pathway)] = "Unclassified"

idx1 = match(Drug_IC50_Info$Drug1, Anno_Drugs$Drug_CID)
idx2 = match(Drug_IC50_Info$Drug2, Anno_Drugs$Drug_CID)

Drug_IC50_Info = Drug_IC50_Info %>% 
  mutate(Drug_Name1 = Anno_Drugs$Name[idx1], 
         Drug_Name2 = Anno_Drugs$Name[idx2], 
         Pathway1 = Anno_Drugs$Target_Pathway[idx1], 
         Pathway2 = Anno_Drugs$Target_Pathway[idx2], 
         Pathway_Equal = ifelse(Pathway1==Pathway2, Pathway1, "Not_Equal"))

mean(Drug_IC50_Info$Pathway_Equal=="Not_Equal", na.rm=T)   # 93.19%


dir = mkdir("../../processed_data/drug_data/DTI_Analysis/Target_Pathway")
option = c("Dist_IC50", "PCC_IC50", "Sim_IC50")
main = sprintf("%s/Target_Pathway & %s", dir, option)

ylab = c(bquote(ln(IC["50"])~Distance~"["*NRMSE*"]"), 
         bquote(ln(IC["50"])~Similarity~"["*PCC*"]"), 
         bquote(ln(IC["50"])~Similarity~"["*RBF*"]"))

# ylab = c(bquote(bold(atop(ln(IC["50"])~Distance~"["*NRMSE*"]"))), 
#          bquote(bold(atop(ln(IC["50"])~Similarity~"["*PCC*"]"))), 
#          bquote(bold(atop(ln(IC["50"])~Similarity~"["*RBF*"]"))))

Drug_IC50_Info %>% 
  boxplot_def(Pathway_Equal, Dist_IC50, main=main[1], 
              ylab=ylab[[1]], legend=F, vjust=1, hjust=1, 
              alpha=0.9, alpha_point=0.16, axis_tl=20, axis_tx=13.5, 
              width=32, height=14.4, dpi=1200, reorder=T, raster=T, save=T, save_svg=T)

Drug_IC50_Info %>% 
  boxplot_def(Pathway_Equal, PCC_IC50, main=main[2], 
              ylab=ylab[[2]], legend=F, vjust=1, hjust=1, 
              alpha=0.9, alpha_point=0.16, axis_tl=22.5, axis_tx=13.5, 
              width=36, height=14.4, dpi=1200, reorder=T, raster=T, save=T, save_svg=T)

Drug_IC50_Info %>% 
  boxplot_def(Pathway_Equal, Sim_IC50, main=main[3], 
              ylab=ylab[[3]], legend=F, vjust=1, hjust=1, 
              alpha=0.9, alpha_point=0.16, axis_tl=22.5, axis_tx=13.5, 
              width=36, height=14.4, dpi=1200, reorder=T, raster=T, save=T, save_svg=T)

col = c("Drug1", "Drug2", "Pathway1", "Pathway2", "Pathway_Equal", "Dist_IC50", "PCC_IC50")
GDSC_DFS_Target_Path_ = Drug_IC50_Info[, col] %>% 
  dplyr::rename(NRMSE_LN_IC50=Dist_IC50, PCC_LN_IC50=PCC_IC50)
  



##### 3. Save files...

merge_ = function(df1, df2) full_join(df1, df2, by=by)

by = c("Drug1"="Drug1", "Drug2"="Drug2")
col = c("Drug1", "Drug2", "Corr_DPI")
GDSC_DFS_RWR_Path = Reduce(merge_, list(GDSC_DFS_ssGSEA_0.3$DFS[, col], GDSC_DFS_ssGSEA_0.5$DFS[, col], 
                                        GDSC_DFS_ssGSEA_0.7$DFS[, col], GDSC_DFS_ssGSEA_0.9$DFS[, col], 
                                        GDSC_DFS_SING_0.3$DFS[, col], GDSC_DFS_SING_0.5$DFS[, col], 
                                        GDSC_DFS_SING_0.7$DFS[, col], GDSC_DFS_SING_0.9$DFS[, col]))

restart = c(0.3, 0.5, 0.7, 0.9)
col_sing = sprintf("PCC_RWR_SING_%s", restart)
col_ssgsea = sprintf("PCC_RWR_ssGSEA_%s", restart)
colnames(GDSC_DFS_RWR_Path)[3:10] = c(col_ssgsea, col_sing)


col = c("Drug1", "Drug2", "Corr_DPI")
GDSC_DFS_RWR_Hetero = Reduce(merge_, list(GDSC_DFS_RWRH_0.3$DFS[, col], GDSC_DFS_RWRH_0.5$DFS[, col], 
                                          GDSC_DFS_RWRH_0.7$DFS[, col], GDSC_DFS_RWRH_0.9$DFS[, col]))

restart = c(0.3, 0.5, 0.7, 0.9)
col = sprintf("PCC_RWR_Hetero_%s", restart)
colnames(GDSC_DFS_RWR_Hetero)[3:6] = col


GDSC_DFS_PathSim = merge_(GDSC_DFS_PathSim0, GDSC_DFS_PathSim1)
colnames(GDSC_DFS_PathSim)[3:4] = c("PathSim_0Hop", "PathSim_1Hop")

GDSC_MoA_Sim = Reduce(merge_, list(Drug_IC50_Info, GDSC_DFS_GO_Seq, 
                                   GDSC_DFS_RWR_Path, GDSC_DFS_RWR_Hetero, GDSC_DFS_PathSim))

GDSC_MoA_Sim = GDSC_MoA_Sim %>% 
  dplyr::rename(Drug_CID1=Drug1, Drug_CID2=Drug2) %>% 
  relocate(Drug_Name1, Drug_Name2, .after=Drug_CID2)   # 186192 [432*431]


dir = "../../processed_data/drug_data/DTI_Analysis"
file = sprintf("%s/GDSC_MoA_Similarity.csv", dir)
write.csv(GDSC_MoA_Sim, file=file, row.names=F)


supplementary = T
if (supplementary) {
  suppressMessages(library(openxlsx))
  # dir = "../../processed_data/drug_data/DTI_Analysis"
  
  ### [Supplementary Data] Supplementary Data 5
  GDSC_MoA_Sim_ = GDSC_MoA_Sim %>% 
    subset(Drug_CID1!=Drug_CID2) %>% 
    dplyr::rename(NRMSE_LN_IC50=Dist_IC50, PCC_LN_IC50=PCC_IC50, SCC_LN_IC50=SCC_IC50,
                 Tanimoto_Morgan_256bit=Tanimoto, Tanimoto_PubChem_881bit=Tanimoto_PC,
                 Target_Pathway1=Pathway1, Target_Pathway2=Pathway2, 
                 Target_Pathway_Equal=Pathway_Equal) %>% 
    subset(select=-Sim_IC50) %>% as.data.frame
  
  gdsc_cids = unique(GDSC_MoA_Sim_$Drug_CID1)
  Drug_Pair = combn(length(gdsc_cids), 2) %>% t
  Drug_Pair = data.frame(Drug_CID1 = gdsc_cids[Drug_Pair[, 1]], 
                         Drug_CID2 = gdsc_cids[Drug_Pair[, 2]])
  
  by = c("Drug_CID1", "Drug_CID2")
  GDSC_MoA_Sim_ = GDSC_MoA_Sim_ %>% right_join(Drug_Pair, by=by)   # 93096 [186192/2]
  GDSC_MoA_Sim_$Drug_CID1 %>% unique %>% length   # 431
  GDSC_MoA_Sim_$Drug_CID2 %>% unique %>% length   # 431
  
  sheets = "Supplementary Data 5"
  # file = sprintf("%s/%s.xlsx", dir, sheets)
  file = sprintf("%s.xlsx", sheets)
  write.xlsx(GDSC_MoA_Sim_, file=file, sheetName=sheets, rowNames=F)
}

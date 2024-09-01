#!/usr/bin/env Rscript

source("../functions.R")
suppressMessages(library(cogena))
suppressMessages(library(reshape2))
loadings()

args = commandArgs(trailingOnly=TRUE)
path_nth = args[1] %>% as.numeric   # Pathway to utilize
# num_k = args[2] %>% as.numeric      # Number of K in KNN
# cores = args[3] %>% as.numeric      # Number of CPU utilized

num_k = 5
cores = T

gsva_def = function(Omics, Path_List, method="gsva", 
                    filt_genes=T, do_qnorm=F, cores=T, ...) {
  
  # Method : gsva or ssgsea
  # Input : Omics [Cell x Genes]
  # Output : Omics [Cell x Pathways]
  # Quantile Normalization is not recommended in ssGSEA
  
  suppressMessages(library(GSVA))
  suppressMessages(library(stats))
  gene_list = Path_List %>% unlist %>% unique
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  
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

calc_sep_score = function(Network, geneset_A, geneset_B) {
  
  # Calculate separation score of a pair of pathways
  # This function requires the installation of dnet, igraph package
  # Users must use the function separation_score, rather than this function
  
  inf_omit = function(x) x[!is.infinite(x)]
  dist_AA = distances(Network, geneset_A, geneset_A, mode="out") %>% as.numeric %>% inf_omit %>% mean
  dist_BB = distances(Network, geneset_B, geneset_B, mode="out") %>% as.numeric %>% inf_omit %>% mean
  dist_AB = distances(Network, geneset_A, geneset_B, mode="out") %>% as.numeric %>% inf_omit %>% mean
  
  if (dist_AA==0 | is.na(dist_AA)) dist_AA = NA
  if (dist_BB==0 | is.na(dist_BB)) dist_BB = NA
  if (is.nan(dist_AB)) dist_AB = NA
  
  # Separation score is only calculated when all dists are calculated
  # If geneset_A are not connected, dist_AA = 0 > NA
  # If geneset_B are not connected, dist_BB = 0 > NA
  # If geneset_AB are not connected, dist_AB = NaN > NA
  
  ss_AB = dist_AB - (dist_AA + dist_BB) / 2
  ss_info = c(ss_AB, dist_AA, dist_BB, dist_AB)
  return(ss_info)
}

separation_score = function(Network, Path_List, min_genes=NULL, verbose=T, cores=T) {
  
  # Calculate separation score of all pairwise pathways
  # This function requires the installation of dnet, igraph package
  # The Network should contain two columns named Node1, Node2
  # Otherwise, the first two columns are recognized as Nodes
  
  # suppressMessages(library(dnet))
  suppressMessages(library(igraph))
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  cond1 = "Node1" %in% colnames(Network)
  cond2 = "Node2" %in% colnames(Network)
  
  if (!cond1 | !cond2) {
    Network = Network[, 1:2]
    colnames(Network) = c("Node1", "Node2")
  }
  
  genes_net = union(Network$Node1, Network$Node2)
  genes_path = Path_List %>% unlist %>% unique
  genes_int = intersect(genes_net, genes_path)
  
  Path_List_Ori = Path_List
  int_ratio = length(genes_int) / length(genes_path)
  Path_List = Path_List %>% sapply(function(x) x[x %in% genes_int])
  
  sprintf("# Pathway size range : [%s, %s]", 
          min(sapply(Path_List, length)), max(sapply(Path_List, length))) %>% print
  
  sprintf("# Total Genes utilized : %.2f%% [%s/%s]", 
          int_ratio*100, length(genes_int), length(genes_path)) %>% print
  
  if (is.numeric(min_genes)) {
    idx = sapply(Path_List, length)>=min_genes
    Path_List = Path_List[idx]
    sprintf("# Total Pathway utilized : %.2f%% [%s/%s]", 
            100 * sum(idx) / length(idx), sum(idx), length(idx)) %>% print
  }
  
  Network = Network %>% dplyr::select(from=Node1, to=Node2)
  Network = Network %>% graph_from_data_frame(directed=T)
  nC2 = combn(length(Path_List), 2) %>% t   
  # [n_path x (n_path-1)/2, 2]
  
  if (cores<2) {
    Network_SS = data.frame()
    for (i in 1:nrow(nC2)) {
      path_A = names(Path_List)[nC2[i, 1]]
      path_B = names(Path_List)[nC2[i, 2]]
      geneset_A = Path_List[[path_A]]
      geneset_B = Path_List[[path_B]]
      
      ss_info = calc_sep_score(Network, geneset_A, geneset_B)
      Network_SS_ = c(path_A, path_B, ss_info)
      Network_SS = Network_SS %>% rbind(Network_SS_)
    }
  } else {
    cluster = multicores(cores=cores)
    Network_SS = foreach (i=1:nrow(nC2), .combine=rbind) %dopar% {
      path_A = names(Path_List)[nC2[i, 1]]
      path_B = names(Path_List)[nC2[i, 2]]
      geneset_A = Path_List[[path_A]]
      geneset_B = Path_List[[path_B]]
      
      ss_info = calc_sep_score(Network, geneset_A, geneset_B)
      Network_SS_ = c(path_A, path_B, ss_info)
      return(Network_SS_)
    }
    Network_SS = Network_SS %>% as.data.frame
    stopCluster(cluster)
  }
  
  rownames(Network_SS) = NULL
  Network_SS[, 3:6] = Network_SS[, 3:6] %>% sapply(as.numeric)
  colnames(Network_SS) = c("Pathway1", "Pathway2", "Dist", "Dist_P1", "Dist_P2", "Dist_P1_P2")
  
  inlen = function(x, y) intersect(x, y) %>% length
  num_overlap = nC2 %>% apply(1, function(x) Reduce(inlen, Path_List_Ori[x]))
  Network_SS$Num_Overlap = num_overlap
  
  Path_Nums = sapply(Path_List_Ori, length)
  idx1 = match(Network_SS$Pathway1, names(Path_Nums))
  idx2 = match(Network_SS$Pathway2, names(Path_Nums))
  
  Network_SS = Network_SS %>% 
    mutate(Num1 = Path_Nums[idx1], Num2 = Path_Nums[idx2],
           Ratio_Overlap = Num_Overlap / (Num1+Num2-Num_Overlap)) %>% 
    relocate(Num_Overlap, .before=Ratio_Overlap)
  
  # Undirect the Pathway-Pathway Interaction Graph
  Network_SS_Rev = Network_SS
  colnames(Network_SS_Rev)[1:2] = c("Pathway2", "Pathway1")
  Network_SS = Network_SS %>% rbind(Network_SS_Rev)
  
  return(Network_SS)
}

knn_graph = function(Network, k=5, mode="min", col_stat="Dist", filt_na=F) {
  
  # Complete KNN graph from Network
  # Network should have the following columns
  # Pathway1, Pathway2, col_stat [Dist or Corr]
  
  Network_Rev = Network %>% dplyr::rename(Pathway1=Pathway2, Pathway2=Pathway1)
  Network = Network %>% rbind(Network_Rev) %>% distinct(Pathway1, Pathway2, .keep_all=T)
  
  filter_net = ifelse(mode=="max", slice_max, slice_min)
  if (filt_na) Network = Network %>% subset(!is.na(object(col_stat)))
  
  Network_KNN = Network %>% group_by(Pathway1) %>% 
    filter_net(object(col_stat), n=k) %>% as.data.frame   # n x (n-1) > n x k
  
  return(Network_KNN)
}

knn_graph_ovl = function(Network, k=5) {
  
  # Complete KNN graph from Network with separation score
  # Network should have the following columns
  # Pathway1, Pathway2, Dist, Ratio_Overlap
  
  # If tied with or do not have separation scores [Dist],
  # Filter the pairs by pathway overlap ratios [Ratio_Overlap]
  
  Network_Rev = Network %>% dplyr::rename(Pathway1=Pathway2, Pathway2=Pathway1)
  Network = Network %>% rbind(Network_Rev) %>% distinct(Pathway1, Pathway2, .keep_all=T)
  
  Network_KNN = Network %>% group_by(Pathway1) %>% 
    arrange(Dist, desc(Ratio_Overlap), .by_group=T) %>% 
    do(head(., n=k)) %>% as.data.frame
  
  sprintf("KNN Edges : %s", nrow(Network_KNN)) %>% print
  return(Network_KNN)
}


# Call Pathway Data
dir = "../../raw_data/MSigDB"

path_list = c("BIOCARTA", "KEGG", "C4_CM", "PID", "WikiPath")
gmt_list = c("c2.cp.biocarta.v2023.1.Hs.entrez.gmt", 
             "c2.cp.kegg.v2023.1.Hs.entrez.gmt", 
             "c4.cm.v2023.1.Hs.entrez.gmt", 
             "c2.cp.pid.v2023.1.Hs.entrez.gmt", 
             "c2.cp.wikipathways.v2023.1.Hs.entrez.gmt")

path_name = path_list[path_nth]
file = sprintf("%s/%s", dir, gmt_list[path_nth])
Path_List = gmt2list(file)

# Call Network Data
dir = "../../processed_data/net_data/RegNetwork"
file = sprintf("%s/RegNetwork_Ent.csv", dir)
Network_Reg = read.csv(file)   # 371910

dir = "../../processed_data/net_data/STRING"
file = sprintf("%s/STRING_Filt_Ent.csv", dir)
Network_STR7 = read.csv(file)   # 479908
Network_STR9 = Network_STR7 %>% subset(Combined_Score>=900)   # 236392


# Implement GSVA in SANGER RNA-Seq Data
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Ent.csv", dir)
SANGER_RNA = read.csv(file, row.names=1, check.names=F)
SANGER_RNA_GSVA = SANGER_RNA %>% gsva_def(Path_List, filt_genes=T, method="gsva", cores=cores)

dir = "../../processed_data/cell_data"
dir = mkdir(sprintf("%s/%s", dir, path_name))
file = sprintf("%s/SANGER_RNA_GSVA.csv", dir)
write.csv(SANGER_RNA_GSVA, file=file, row.names=T)


# GSVA Correlation KNN Graph
RNA_Corr = SANGER_RNA_GSVA %>% cor %>% reshape2::melt() %>% as.data.frame
colnames(RNA_Corr) = c("Pathway1", "Pathway2", "Corr")
RNA_Corr = RNA_Corr %>% subset(Pathway1!=Pathway2)

RNA_Corr_KNN3 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=3)
RNA_Corr_KNN5 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=5)
RNA_Corr_KNN7 = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=7)

dir = "../../processed_data/net_data"
dir = mkdir(sprintf("%s/%s", dir, path_name))
file = sprintf("%s/KNN%s_RNA_Corr.csv", dir, c("", 3, 5, 7))

write.csv(RNA_Corr, file=file[1], row.names=F)
write.csv(RNA_Corr_KNN3, file=file[2], row.names=F)
write.csv(RNA_Corr_KNN5, file=file[3], row.names=F)
write.csv(RNA_Corr_KNN7, file=file[4], row.names=F)


# Pathway-Network KNN Graph
Network_Reg$Node1 = Network_Reg$Node1 %>% as.character
Network_Reg$Node2 = Network_Reg$Node2 %>% as.character
Network_STR7$Node1 = Network_STR7$Node1 %>% as.character
Network_STR7$Node2 = Network_STR7$Node2 %>% as.character
Network_STR9$Node1 = Network_STR9$Node1 %>% as.character
Network_STR9$Node2 = Network_STR9$Node2 %>% as.character

Path_Score_Reg = separation_score(Network_Reg, Path_List, min_genes=NULL, cores=cores)
Path_Score_STR7 = separation_score(Network_STR7, Path_List, min_genes=NULL, cores=cores)
Path_Score_STR9 = separation_score(Network_STR9, Path_List, min_genes=NULL, cores=cores)

net_name = c("RegNetwork", "STRING_700", "STRING_900")
KNN_Graph_Reg = Path_Score_Reg %>% knn_graph_ovl(k=num_k)
KNN_Graph_STR7 = Path_Score_STR7 %>% knn_graph_ovl(k=num_k)
KNN_Graph_STR9 = Path_Score_STR9 %>% knn_graph_ovl(k=num_k)


dir = "../../processed_data/net_data"
dir = mkdir(sprintf("%s/%s", dir, path_name))

file = sprintf("%s/SS_%s.csv", dir, net_name)
write.csv(Path_Score_Reg, file=file[1], row.names=F)
write.csv(Path_Score_STR7, file=file[2], row.names=F)
write.csv(Path_Score_STR9, file=file[3], row.names=F)

file = sprintf("%s/KNN%s_%s.csv", dir, num_k, net_name)
write.csv(KNN_Graph_Reg, file=file[1], row.names=F)
write.csv(KNN_Graph_STR7, file=file[2], row.names=F)
write.csv(KNN_Graph_STR9, file=file[3], row.names=F)


# Graph for RGCN
col = c("Pathway1", "Pathway2")
KNN_Graph_Reg_ = KNN_Graph_Reg[, col] %>% mutate(Edge_Type="RegNetwork")
KNN_Graph_STR7_ = KNN_Graph_STR7[, col] %>% mutate(Edge_Type="STRING_700")
KNN_Graph_STR9_ = KNN_Graph_STR9[, col] %>% mutate(Edge_Type="STRING_900")

RNA_Corr_KNN_ = RNA_Corr %>% knn_graph(col_stat="Corr", mode="max", k=num_k)
RNA_Corr_KNN_ = RNA_Corr_KNN_[, col] %>% mutate(Edge_Type="RNA_Corr")

KNN_STR7_Reg_Corr = Reduce(rbind, list(KNN_Graph_Reg_, KNN_Graph_STR7_, RNA_Corr_KNN_))
KNN_STR9_Reg_Corr = Reduce(rbind, list(KNN_Graph_Reg_, KNN_Graph_STR9_, RNA_Corr_KNN_))

file1 = sprintf("%s/KNN%s_STR7_Reg_Corr.csv", dir, num_k)
file2 = sprintf("%s/KNN%s_STR9_Reg_Corr.csv", dir, num_k)
write.csv(KNN_STR7_Reg_Corr, file=file1, row.names=F)
write.csv(KNN_STR9_Reg_Corr, file=file2, row.names=F)

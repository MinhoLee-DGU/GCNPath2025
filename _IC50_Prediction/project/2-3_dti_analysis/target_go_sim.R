#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(coop))
suppressMessages(library(cogena))
suppressMessages(library(reshape2))

source("../functions.R")
loadings()

dir = "../../processed_data/drug_data/STiTCH"
file = sprintf("%s/GDSC_DTI.csv", dir)
GDSC_DTI = read.csv(file)

GDSC_DTI_List = split(GDSC_DTI$Target, GDSC_DTI$Drug_CID)
GDSC_DFS_GO = combn(length(GDSC_DTI_List), 2) %>% t %>% as.data.frame

colnames(GDSC_DFS_GO) = c("Drug1", "Drug2")
GDSC_DFS_GO$Drug1 = names(GDSC_DTI_List)[GDSC_DFS_GO$Drug1]
GDSC_DFS_GO$Drug2 = names(GDSC_DTI_List)[GDSC_DFS_GO$Drug2]



### 2-4. Target GO semantic similarity

sim_go = function(DFS_Info, DTI_List, GOBP, GOMF, GOCC, cores=12) {
  
  cluster = multicores(cores=cores, type="PSOCK")
  on.exit(stopCluster(cluster))
  
  Sim_GO = foreach(i=1:nrow(DFS_Info), .combine=rbind, 
                   .packages=c("GOSemSim", "org.Hs.eg.db"),
                   .export=c("DFS_Info", "DTI_List", "GOBP", "GOMF", "GOCC")) %dopar% {
     
    cid1 = DFS_Info$Drug1[i]
    cid2 = DFS_Info$Drug2[i]
    
    sim_gobp = clusterSim(DTI_List[[cid1]], DTI_List[[cid2]], 
                          semData=GOBP, measure="Rel", combine="BMA")
    sim_gomf = clusterSim(DTI_List[[cid1]], DTI_List[[cid2]], 
                          semData=GOMF, measure="Rel", combine="BMA")
    sim_gocc = clusterSim(DTI_List[[cid1]], DTI_List[[cid2]], 
                          semData=GOCC, measure="Rel", combine="BMA")
    
    c(sim_gobp, sim_gomf, sim_gocc)
  }
  
  rownames(Sim_GO) = NULL
  colnames(Sim_GO) = c("Sim_GOBP", "Sim_GOMF", "Sim_GOCC")
  
  DFS_Info = DFS_Info %>% cbind(Sim_GO)
  return(DFS_Info)
}


suppressMessages(library(GOSemSim))
suppressMessages(library(org.Hs.eg.db))

GOBP_Human = godata("org.Hs.eg.db", keytype="SYMBOL", ont="BP")
GOMF_Human = godata("org.Hs.eg.db", keytype="SYMBOL", ont="MF")
GOCC_Human = godata("org.Hs.eg.db", keytype="SYMBOL", ont="CC")

# GOBP_Human@keys %>% length   # 66091
# GDSC_DTI_List %>% sapply(function(x) any(x %in% GOBP_Human@keys)) %>% sum   # 300
# GDSC_DTI_List %>% sapply(function(x) any(x %in% GOMF_Human@keys)) %>% sum   # 300
# GDSC_DTI_List %>% sapply(function(x) any(x %in% GOCC_Human@keys)) %>% sum   # 300


cores = T
GDSC_DFS_GO = GDSC_DFS_GO %>% 
  sim_go(GDSC_DTI_List, GOBP_Human, GOMF_Human, GOCC_Human, cores=cores)

dir = mkdir("../../processed_data/drug_data/DTI_Analysis/Target_GO")
file = sprintf("%s/MoA_Similarity.csv", dir)
write.csv(GDSC_DFS_GO, file=file, row.names=F)

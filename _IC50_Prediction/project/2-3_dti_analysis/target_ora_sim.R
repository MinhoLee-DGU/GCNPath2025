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

dir = "../../processed_data/net_data/STRING"
file = sprintf("%s/STRING_Filt_Sym.csv", dir)
STRING = fread(file)
STRING = STRING %>% as.data.frame

GDSC_DTI_List = split(GDSC_DTI$Target, GDSC_DTI$Drug_CID)
GDSC_DFS_STiTCH = combn(length(GDSC_DTI_List[1:3]), 2) %>% t %>% as.data.frame

colnames(GDSC_DFS_STiTCH) = c("Drug1", "Drug2")
GDSC_DFS_STiTCH$Drug1 = names(GDSC_DTI_List)[GDSC_DFS_STiTCH$Drug1]
GDSC_DFS_STiTCH$Drug2 = names(GDSC_DTI_List)[GDSC_DFS_STiTCH$Drug2]



### 2-6. ORA of Target + 1hop
# [2020] DrugSimDB, A comprehensive integrated drug similarity resource for in-silico drug repositioning and beyond
# doi.org/10.1093/bib/bbaa126

suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))

get_ora = function(geneset, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont="BP", 
                   pval=0.05, qval=0.05, min_genes=10, max_genes=100) {
  
  ORA = enrichGO(gene=geneset, OrgDb=OrgDb, keyType=keyType, 
                 ont=ont, pAdjustMethod="BH", 
                 pvalueCutoff=pval, qvalueCutoff=qval, 
                 minGSSize=min_genes, maxGSSize=max_genes)
  
  ora_id = ORA@result %>% subset(pvalue<pval & qvalue<qval) %>% pull(ID)
  return(ora_id)
}

sim_ora = function(DFS_Info, DTI_List, STRING, GOBP, GOMF, GOCC, cores=T, ...) {
  
  ORA_1Hop_GOBP = list()
  ORA_1Hop_GOMF = list()
  ORA_1Hop_GOCC = list()
  
  for (i in 1:length(DTI_List)) {
    cid = names(DTI_List)[i]
    targets = DTI_List[[i]]
    STRING_Temp = STRING %>% subset(Node1 %in% targets | Node2 %in% targets)
    targets_1hop = union(STRING_Temp$Node1, STRING_Temp$Node2)
    
    ORA_1Hop_GOBP[[cid]] = get_ora(targets_1hop, OrgDb=org.Hs.eg.db, ont="BP", ...)
    ORA_1Hop_GOMF[[cid]] = get_ora(targets_1hop, OrgDb=org.Hs.eg.db, ont="MF", ...)
    ORA_1Hop_GOCC[[cid]] = get_ora(targets_1hop, OrgDb=org.Hs.eg.db, ont="CC", ...)
    
    sprintf("# Total %s targets+1hop in %s", length(targets_1hop), cid) %>% print
    sprintf("# Total %s GO-BP in %s", length(ORA_1Hop_GOBP[[cid]]), cid) %>% print
    sprintf("# Total %s GO-MF in %s", length(ORA_1Hop_GOMF[[cid]]), cid) %>% print
    sprintf("# Total %s GO-CC in %s", length(ORA_1Hop_GOCC[[cid]]), cid) %>% print
  }
  
  num_gobp = ORA_1Hop_GOBP %>% sapply(length)
  num_gomf = ORA_1Hop_GOMF %>% sapply(length)
  num_gocc = ORA_1Hop_GOCC %>% sapply(length)
  
  paste_ = function(x) paste0(x, collapse="|")
  gobp_list = ORA_1Hop_GOBP %>% sapply(paste_)
  gomf_list = ORA_1Hop_GOMF %>% sapply(paste_)
  gocc_list = ORA_1Hop_GOCC %>% sapply(paste_)
  
  ORA_Info = data.frame(Drug_CID=names(DTI_List),
                        Num_GOBP=num_gobp, Num_GOMF=num_gomf, Num_GOCC=num_gocc, 
                        GOBP=gobp_list, GOMF=gomf_list, GOCC=gocc_list)
  
  
  cluster = multicores(cores=cores, type="PSOCK")
  on.exit(stopCluster(cluster))
  
  Sim_ORA = foreach(i=1:nrow(DFS_Info), .combine=rbind, 
                    .packages=c("GOSemSim", "org.Hs.eg.db"),
                    .export=c("ORA_1Hop_GOBP", "ORA_1Hop_GOMF", "ORA_1Hop_GOCC", 
                              "DFS_Info", "GOBP", "GOMF", "GOCC")) %dopar% {
                                
    cid1 = DFS_Info$Drug1[i]
    cid2 = DFS_Info$Drug2[i]
    
    sim_gobp = mgoSim(ORA_1Hop_GOBP[[cid1]], ORA_1Hop_GOBP[[cid2]], 
                      semData=GOBP, measure="Rel", combine="BMA")
    sim_gomf = mgoSim(ORA_1Hop_GOMF[[cid1]], ORA_1Hop_GOMF[[cid2]], 
                      semData=GOMF, measure="Rel", combine="BMA")
    sim_gocc = mgoSim(ORA_1Hop_GOCC[[cid1]], ORA_1Hop_GOCC[[cid2]], 
                      semData=GOCC, measure="Rel", combine="BMA")
    
    c(sim_gobp, sim_gomf, sim_gocc)
  }
  
  rownames(Sim_ORA) = NULL
  colnames(Sim_ORA) = c("Sim_GOBP_ORA", "Sim_GOMF_ORA", "Sim_GOCC_ORA")
  DFS_Info = DFS_Info %>% cbind(Sim_ORA)
  
  DFS_Info = list(DFS_Info=DFS_Info, ORA_Info=ORA_Info)
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
GDSC_DFS_STiTCH = GDSC_DFS_STiTCH %>% 
  sim_ora(GDSC_DTI_List[1:3], STRING, GOBP_Human, GOMF_Human, GOCC_Human, cores=cores)

GDSC_ORA_STiTCH = GDSC_DFS_STiTCH$ORA_Info
GDSC_DFS_STiTCH = GDSC_DFS_STiTCH$DFS_Info

dir = mkdir("../../processed_data/drug_data/DTI_Analysis/Target_ORA")

file = sprintf("%s/ORA_Info.csv", dir)
write.csv(GDSC_ORA_STiTCH, file=file, row.names=F)

file = sprintf("%s/MoA_Similarity.csv", dir)
write.csv(GDSC_DFS_STiTCH, file=file, row.names=F)

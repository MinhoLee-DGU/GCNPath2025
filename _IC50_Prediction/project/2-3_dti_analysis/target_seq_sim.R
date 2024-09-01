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

file = sprintf("%s/STRING_Target.csv", dir)
STRING_Target = read.csv(file)


GDSC_DTI_List = split(GDSC_DTI$Target, GDSC_DTI$Drug_CID)
GDSC_DFS_STiTCH = combn(length(GDSC_DTI_List), 2) %>% t %>% as.data.frame

colnames(GDSC_DFS_STiTCH) = c("Drug1", "Drug2")
GDSC_DFS_STiTCH$Drug1 = names(GDSC_DTI_List)[GDSC_DFS_STiTCH$Drug1]
GDSC_DFS_STiTCH$Drug2 = names(GDSC_DTI_List)[GDSC_DFS_STiTCH$Drug2]



### 2-5. Target Seq Similarity
# Calculate normalized sequence similarities between targets
# Smith-Waterman [Local] & Needleman-Wunsch [Global]
# [2011] SITAR, Combining Drug and Gene Similarity Measures for Drug-Target Elucidation
# doi.org/10.1089/cmb.2010.0213

sw_norm = function(seq1, seq2, submat, type="local", gap_open=10, gap_ext=4, score_only=T) {
  
  seq1 = AAString(seq1)
  seq2 = AAString(seq2)
  
  sw12 = pairwiseAlignment(seq1, seq2, substitutionMatrix=submat, type=type, 
                           gapOpening=gap_open, gapExtension=gap_ext, scoreOnly=score_only)
  
  sw11 = pairwiseAlignment(seq1, seq1, substitutionMatrix=submat, type=type, 
                           gapOpening=gap_open, gapExtension=gap_ext, scoreOnly=score_only)
  
  sw22 = pairwiseAlignment(seq2, seq2, substitutionMatrix=submat, type=type, 
                           gapOpening=gap_open, gapExtension=gap_ext, scoreOnly=score_only)
  
  return(sw12/(sqrt(sw11)*sqrt(sw22)))
}

cluster_sw_norm = function(seq1_list, seq2_list, submat, ...) {
  
  sw = c()
  for (seq1 in seq1_list) {
    for (seq2 in seq2_list) {
      sw = sw %>% c(sw_norm(seq1, seq2, submat, ...))
    }
  }
  
  # Target Similarity by BMA [Best Match Average]
  SW = sw %>% matrix(length(seq1_list), length(seq2_list), byrow=T)
  sw = c(apply(SW, 1, max), apply(SW, 2, max)) %>% mean
  return(sw)
}


suppressMessages(library(ensembldb))
suppressMessages(library(Biostrings))
suppressMessages(library(EnsDb.Hsapiens.v86))

# Get protein sequences of GDSC targets
gdsc_target = GDSC_DTI$Target_Symbol %>% unique   # 2188
idx = match(gdsc_target, STRING_Target$SYMBOL)
gdsc_target_ens = STRING_Target$ENSEMBL_ID[idx]

GDSC_Target_Ens = data.frame(Target=gdsc_target, Ensembl=gdsc_target_ens)
GDSC_Target_Ens$Ensembl %>% is.na %>% sum   # 0

edb = EnsDb.Hsapiens.v86
Ensembl = proteins(edb, filter=ProteinIdFilter(GDSC_Target_Ens$Ensembl), return.type="DataFrame")

idx = match(GDSC_Target_Ens$Ensembl, Ensembl$protein_id)
GDSC_Target_Ens$Sequence = Ensembl$protein_sequence[idx]
GDSC_Target_Ens$Ensembl[is.na(GDSC_Target_Ens$Sequence)]
# ENSP00000378704  ENSP00000387888  ENSP00000368924  ENSP00000387760  ENSP00000483056


# Get protein sequences of GDSC targets [Manual]
dir = "../../processed_data/drug_data/DTI_Analysis"
file = sprintf("%s/Ensembl_NA.txt", dir)
Ensembl_NA = read.table(file, sep="\t", header=F)

colnames(Ensembl_NA) = c("Ensembl", "Sequence")
idx = match(GDSC_Target_Ens$Ensembl, Ensembl_NA$Ensembl)
GDSC_Target_Ens$Sequence = ifelse(is.na(idx), GDSC_Target_Ens$Sequence, Ensembl_NA$Sequence[na.omit(idx)])
GDSC_Target_Ens$Sequence %>% is.na %>% sum   # 0


# Calculate sequence similarity [Smith-Waterman]
sim_seq = function(DFS_Info, DTI_List, Target_Ens, SubMat=data("BLOSUM62"),
                   col_seq="Sequence", delete_ukn_aa=T, cores=12) {
  
  if (delete_ukn_aa) {
    aa_submat = SubMat %>% rownames
    aa_seq = Target_Ens[[col_seq]] %>% strsplit("") %>% unlist %>% unique
    aa_seq_ukn = aa_seq[!(aa_seq %in% aa_submat)]   # U [selenocysteine]
    
    if (aa_seq_ukn>0) {
      aa_ukn = paste0(aa_seq_ukn, collapse="|")
      aa_ukn_ = paste(aa_seq_ukn, collapse=",")
      Target_Ens[[col_seq]] = gsub(aa_ukn, "", Target_Ens[[col_seq]])
      sprintf("# Detected unknown amino acids not in SubMat... [%s]", aa_ukn_) %>% print
    } else print("# No unknown amino acids to detect...")
  }
  
  cluster = multicores(cores=cores)
  on.exit(stopCluster(cluster))
  
  Sim_Seq = foreach(i=1:nrow(DFS_Info), .combine=rbind) %dopar% {
    cid1 = DFS_Info$Drug1[i]
    cid2 = DFS_Info$Drug2[i]
    target1_list = DTI_List[[cid1]]
    target2_list = DTI_List[[cid2]]
    
    seq1_list = Target_Ens[[col_seq]][Target_Ens$Target %in% target1_list]
    seq2_list = Target_Ens[[col_seq]][Target_Ens$Target %in% target2_list]
    sim_seq_local = cluster_sw_norm(seq1_list, seq2_list, SubMat, type="local")
    sim_seq_global = cluster_sw_norm(seq1_list, seq2_list, SubMat, type="global")
    c(sim_seq_local, sim_seq_global)
  }
  
  rownames(Sim_Seq) = NULL
  colnames(Sim_Seq) = c("Sim_Seq_Local", "Sim_Seq_Global")
  DFS_Info = DFS_Info %>% cbind(Sim_Seq)
  return(DFS_Info)
}


cores = T
data("BLOSUM62")
GDSC_DFS_STiTCH = GDSC_DFS_STiTCH %>% 
  sim_seq(GDSC_DTI_List, GDSC_Target_Ens, BLOSUM62, cores=cores)

dir = mkdir("../../processed_data/drug_data/DTI_Analysis/Target_Seq")
file = sprintf("%s/MoA_Similarity.csv", dir)
write.csv(GDSC_DFS_STiTCH, file=file, row.names=F)


supplementary = T
if (supplementary) {
  suppressMessages(library(openxlsx))
  dir = "../../processed_data/drug_data/DTI_Analysis"
  
  ### [Supplementary Data] Supplementary Data 4
  idx = match(GDSC_DTI$Target_Symbol, GDSC_Target_Ens$Target)
  GDSC_DTI_ = GDSC_DTI %>% mutate(Target_Sequence=GDSC_Target_Ens$Sequence[idx])
  
  GDSC_DTI_$Drug_CID %>% unique %>% length        # 302
  GDSC_DTI_$Target_Symbol %>% unique %>% length   # 2188
  GDSC_DTI_$Target_Sequence %>% is.na %>% sum     # 0
  GDSC_DTI_[, c("Drug_CID", "Target_Symbol")] %>% distinct %>% dim   # 6777
  
  sheets = "Supplementary Data 4"
  file = sprintf("%s/%s.xlsx", dir, sheets)
  write.xlsx(GDSC_DTI_, file=file, sheetName=sheets, rowNames=F)
}

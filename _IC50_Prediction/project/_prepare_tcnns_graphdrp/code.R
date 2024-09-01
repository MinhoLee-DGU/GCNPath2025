#!/usr/bin/env Rscript

##### 1. Preparation

source("../functions.R")
suppressMessages(library(reshape2))
loadings()

# dir = "../../processed_data/cell_data/SANGER_Passports"
# file = sprintf("%s/Anno_Cells.csv", dir)
# Anno_Cells = read.csv(file)



##### 2. Process MUT & CNV [tCNNS, GraphDRP]

dir = "../../raw_data/GDSC"
file1 = sprintf("%s/PANCANCER_Genetic_features_Tue Jun 20 06_41_49 2023.csv", dir)
file2 = sprintf("%s/PANCANCER_Genetic_features_Tue Jun 20 06_44_02 2023.csv", dir)

MUT_CNV1 = fread(file1)   # 699665
MUT_CNV2 = fread(file2)   # 699665

colnames(MUT_CNV1) = gsub(" ", "_", colnames(MUT_CNV1))
colnames(MUT_CNV2) = gsub(" ", "_", colnames(MUT_CNV2))

# Genetic Features of GDSC1 and GDSC2 were the same!
identical(MUT_CNV1, MUT_CNV2)   # T
MUT_CNV = MUT_CNV1
MUT_CNV[, .N, by=.(IS_Mutated, Recurrent_Gain_Loss)]
MUT_CNV$IS_Mutated %>% is.na %>% sum   # 0

MUT_CNV = dcast(MUT_CNV, COSMIC_ID~Genetic_Feature, value.var="IS_Mutated", fill=0)
MUT_CNV = MUT_CNV %>% as.data.frame
rownames(MUT_CNV) = MUT_CNV$COSMIC_ID
MUT_CNV = MUT_CNV[, -1]   # 970 x 735


# Get the order of cell features in GraphDRP 
# First, get the "mut_dict" file with preprocess.py in GraphDRP github code
# And convert the "mut_dict" object into TXT format [mut_dict.py]

GraphDRP_Feat = read.csv("GraphDRP_Feat.txt", sep="\t")
identical(GraphDRP_Feat$Index, 0:734)   # T
all(colnames(MUT_CNV) %in% GraphDRP_Feat$Feature)   # F

# All feats are converted from . into -
setdiff(colnames(MUT_CNV), GraphDRP_Feat$Feature)
setdiff(GraphDRP_Feat$Feature, colnames(MUT_CNV))

col = gsub("\\.", "-", colnames(MUT_CNV))
all(col %in% GraphDRP_Feat$Feature)   # F
setdiff(col, GraphDRP_Feat$Feature)   # EWSR1-ERG_mut
setdiff(GraphDRP_Feat$Feature, col)   # EWRS1-ERG_mut

# EWSR1-ERG is right [ERSR1 = EWS RNA Binding Protein 1]
GraphDRP_Feat$Feature[GraphDRP_Feat$Feature=="EWRS1-ERG_mut"] = "EWSR1-ERG_mut"
all(col %in% GraphDRP_Feat$Feature)   # T
colnames(MUT_CNV) = col
MUT_CNV = MUT_CNV[, GraphDRP_Feat$Feature]

dir = mkdir("../../do_better/_materials/GraphDRP")
file = sprintf("%s/Genetic_Features.csv", dir)
write.csv(MUT_CNV, file=file, row.names=T)



### 3. Similarity with the original data

file = sprintf("PANCANCER_Genetic_feature.csv")   # Original data from GraphDRP/Data
MUT_CNV_Ori = fread(file)   # 714055

MUT_CNV_Ori = dcast(MUT_CNV_Ori, cosmic_sample_id~genetic_feature, value.var="is_mutated")
MUT_CNV_Ori = MUT_CNV_Ori %>% as.data.frame
rownames(MUT_CNV_Ori) = MUT_CNV_Ori$cosmic_sample_id

MUT_CNV_Ori = MUT_CNV_Ori[, -1]   # 990 x 735
GraphDRP_Feat_ = GraphDRP_Feat
GraphDRP_Feat_$Feature[GraphDRP_Feat_$Feature=="EWSR1-ERG_mut"] = "EWRS1-ERG_mut"

MUT_CNV_Ori = MUT_CNV_Ori[, GraphDRP_Feat_$Feature]
colnames(MUT_CNV_Ori) = GraphDRP_Feat$Feature
identical(colnames(MUT_CNV_Ori), colnames(MUT_CNV))   # T


# All identical in intesecting cell-lines!!!
#           Reference
# Prediction      0      1
#          0 670103      0
#          1      0  28092

cells = intersect(rownames(MUT_CNV_Ori), rownames(MUT_CNV))   # 968
data_ref = MUT_CNV_Ori[cells, ] %>% unlist %>% factor(levels=c(0, 1))
data_new = MUT_CNV[cells, ] %>% unlist %>% factor(levels=c(0, 1))
caret::confusionMatrix(data_ref, data_new)

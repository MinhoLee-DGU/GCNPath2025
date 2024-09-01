#!/usr/bin/env Rscript

##### 1. Preparation

source("../functions.R")
suppressMessages(library(cogena))
suppressMessages(library(graphite))
suppressMessages(library(reshape2))
loadings()

dir = "../../processed_data/cell_data/SANGER_Passport"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)



##### 2-1. Process cell network data [STRING]

dir = "../../raw_data/STRING"
file1 = sprintf("%s/9606.protein.links.v11.5.txt", dir)
file2 = sprintf("%s/9606.protein.info.v11.5.txt", dir)

STRING = read.csv(file1, header=T, sep=" ")
STRING_Anno = read.csv(file2, header=T, sep="\t", quote="")

colnames(STRING) = c("Node1", "Node2", "Combined_Score")
colnames(STRING_Anno)[1:2] = c("ENSEMBL_ID", "SYMBOL")
colnames(STRING_Anno) = colnames(STRING_Anno) %>% toupper
STRING_Filt = STRING %>% subset(Combined_Score>=700)   # 11938498 > 505968


# Ensembl ID > Symbol
idx1 = match(STRING_Filt$Node1, STRING_Anno$ENSEMBL_ID)
idx2 = match(STRING_Filt$Node2, STRING_Anno$ENSEMBL_ID)

STRING_Filt_Sym = STRING_Filt %>% mutate(Node1=STRING_Anno$SYMBOL[idx1],
                                         Node2=STRING_Anno$SYMBOL[idx2])
STRING_Filt_Sym = STRING_Filt_Sym[complete.cases(STRING_Filt_Sym), ]   # 505968 > 505968

# Ensembl ID > Entrez ID
idx1 = match(STRING_Filt_Sym$Node1, Anno_Genes$HGNC_SYMBOL)
idx2 = match(STRING_Filt_Sym$Node2, Anno_Genes$HGNC_SYMBOL)

STRING_Filt_Ent = STRING_Filt_Sym %>% mutate(Node1=Anno_Genes$ENTREZ_ID[idx1] %>% as.character,
                                             Node2=Anno_Genes$ENTREZ_ID[idx2] %>% as.character)

STRING_Filt_Ent = STRING_Filt_Ent[complete.cases(STRING_Filt_Ent), ]   # 505968 > 479204

union(STRING_Filt_Sym$Node1, STRING_Filt_Sym$Node2) %>% length   # 16812
union(STRING_Filt_Ent$Node1, STRING_Filt_Ent$Node2) %>% length   # 16327

# Save the files
dir = mkdir("../../processed_data/net_data/STRING")
file1 = sprintf("%s/STRING_Filt_Sym.csv", dir)
file2 = sprintf("%s/STRING_Filt_Ent.csv", dir)
file3 = sprintf("%s/Anno_Genes.csv", dir)

write.csv(STRING_Filt_Sym, row.names=F, file=file1)
write.csv(STRING_Filt_Ent, row.names=F, file=file2)
write.csv(STRING_Anno, row.names=F, file=file3)

file = sprintf("%s/STRING.RData", dir)
save(STRING_Filt_Sym, STRING_Filt_Ent, STRING_Anno, file=file)



##### 2-2. Process cell network data [RegNetwork]

dir = "../../raw_data/RegNetwork"
file = sprintf("%s/human.source", dir)
RegNetwork = read.csv(file, sep="\t", header=F)   # 372774

dir = "../../raw_data/RegNetwork"
file = sprintf("%s/new_kegg.human.reg.direction.txt", dir)
RegNetwork_New = read.csv(file, sep=" ")   # 3954

RegNetwork_Sym = RegNetwork[, c(1, 3)]
RegNetwork_Ent = RegNetwork[, c(2, 4)]
colnames(RegNetwork_Sym) = c("Node1", "Node2")
colnames(RegNetwork_Ent) = c("Node1", "Node2")

RegNetwork_New_Sym = RegNetwork_New[, c(1, 3)]   
RegNetwork_New_Ent = RegNetwork_New[, c(2, 4)]
colnames(RegNetwork_New_Sym) = c("Node1", "Node2")
colnames(RegNetwork_New_Ent) = c("Node1", "Node2")

# Intersection between RegNetwork & STRING
# How much information redundancy does RegNetwork have? [very little]
RegNetwork_Sym = RegNetwork_Sym %>% rbind(RegNetwork_New_Sym) %>% distinct   # 374180
RegNetwork_Ent = RegNetwork_Ent %>% rbind(RegNetwork_New_Ent) %>% distinct   # 371910

# Genes are abundant in both datasets
genes_regnet = union(RegNetwork_Ent$Node1, RegNetwork_Ent$Node2)     # 23101
genes_string = union(STRING_Filt_Ent$Node1, STRING_Filt_Ent$Node2)   # 16082
intersect(genes_regnet, genes_string) %>% length   # 15354

# RegNetwork is very different from STRING
dplyr::intersect(RegNetwork_Ent, STRING_Filt_Ent[, 1:2]) %>% nrow   # 10036  
# 10036 / 371910 [2.70%, RegNetwork]
# 10036 / 479204 [2.09%, STRING_700]

# Save files
dir = mkdir("../../processed_data/net_data/RegNetwork")
file1 = sprintf("%s/RegNetwork_Sym.csv", dir)
file2 = sprintf("%s/RegNetwork_Ent.csv", dir)

write.csv(RegNetwork_Sym, row.names=F, file=file1)
write.csv(RegNetwork_Ent, row.names=F, file=file2)

file = sprintf("%s/RegNetwork.RData", dir)
save(RegNetwork_Sym, RegNetwork_Ent, file=file)

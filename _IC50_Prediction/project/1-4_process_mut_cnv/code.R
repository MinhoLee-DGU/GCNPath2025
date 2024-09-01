#!/usr/bin/env Rscript

##### 1. Preparation

suppressMessages(library(cogena))
source("../functions.R")
loadings()

dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/mutations_all_20230202.csv", dir)
MUT = fread(file, na.strings="")   # 10050692 x 13

file = sprintf("%s/WES_pureCN_CNV_genes_20221213.csv", dir)
CNV = fread(file, na.strings=c("", "NA"))

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Genes.csv", dir)
Anno_Genes = read.csv(file)

# dir = "../../processed_data/path_data"
# file = sprintf("%s/KEGG_Total_Ent.gmt", dir)
# KEGG_List = gmt2list(file)
# geneset = KEGG_List %>% unlist %>% unique   # 5365



##### 2-1. Mutation


MUT_List$model_id %>% ulen   # 1357
MUT_List$model_id %>% table %>% mean %>% round   # 5309
# MUT_List$model_id %>% table %>% hist


breaks = seq(0, 1, 0.05)
dir = mkdir("~/GDSC_Last/Intermediate_1_Cell/MUT")

file = sprintf("%s/MUT.RData", dir)
file1 = sprintf("%s/MUT.csv", dir)
file2 = sprintf("%s/MUT_Ent.csv", dir)

write.csv(MUT, file1, quote=F)
write.csv(MUT_Ent, file2, quote=F)
save(MUT, MUT_Ent, file=file)



##### 2-2. CNV

into_num = function(df, fill_NA=NULL) {
  rows = rownames(df)
  df = df %>% apply(2, as.numeric) %>% as.data.frame
  if(!is.null(fill_NA)) df[is.na(df)] = fill_NA
  rownames(df) = rows
  return(df)
}

rownames(CNV) = CNV$V1
CNV = CNV[-4, -1]

rownames(CNV_C) = CNV_C$V1
CNV_C = CNV_C[-4, -1]

CNV = CNV %>% t %>% as.data.frame
CNV = CNV[!duplicated(CNV$model_id), ]
rownames(CNV) = CNV$model_id
# [Sanger ID] 1371 > 1265
CNV = CNV[, -c(1:3)]

CNV_C = CNV_C %>% t %>% as.data.frame
CNV_C = CNV_C[!duplicated(CNV_C$model_id), ]
rownames(CNV_C) = CNV_C$model_id
# [Sanger ID] 1371 > 1265
CNV_C = CNV_C[, -c(1:3)]

CNV_C[CNV_C=="Deletion"] = 0
CNV_C[CNV_C=="Loss"] = 1
CNV_C[CNV_C=="Neutral"] = 2
CNV_C[CNV_C=="Gain"] = 3
CNV_C[CNV_C=="Amplification"] = 4

CNV[] = CNV %>% sapply(as.numeric) %>% as.data.frame
CNV_C[] = CNV_C %>% sapply(as.numeric) %>% as.data.frame
# The CNV data is chosen considering each cell's ploidy information
# To prevent the outliers, CNV categorization is recommended 

# Val = 2log2(C/Ploidy) * 2
# if Val == 0: Category = 'Deletion'
# if Val == 1: Category = 'Loss'
# if Val == 2: Category = 'Neutral'
# if Val == 3: Category = 'Gain'
# if Val == 4: Category = 'Amplification'

# https://www.biostars.org/p/120741/
# https://www.nature.com/articles/s41467-019-13779-x
# https://depmap.sanger.ac.uk/documentation/datasets/copy-number/
# https://forum.depmap.org/t/what-is-relative-copy-number-copy-number-ratio/104/5

breaks = seq(0, 1300, 1)

CNV = CNV %>% into_num
CNV %>% is.na %>% sum  # 573363
CNV %>% apply(2, na) %>% hist(breaks=breaks)

CNV_C = CNV_C %>% into_num
CNV_C %>% is.na %>% sum  # 385703
CNV_C %>% apply(2, na) %>% hist(breaks=breaks)

genes_NA = CNV %>% apply(2, na)
genes_NA = genes_NA[genes_NA!=0] %>% names    # 777
CNV = CNV[, !(colnames(CNV) %in% genes_NA)]   # 18801 > 17724
CNV_Ent = CNV %>% into_entrez(Anno_Genes)     # 17224 > 16877

genes_NA = CNV_C %>% apply(2, na)
genes_NA = genes_NA[genes_NA!=0] %>% names          # 1124
CNV_C = CNV_C[, !(colnames(CNV_C) %in% genes_NA)]   # 18348 > 17724
CNV_C_Ent = CNV_C %>% into_entrez(Anno_Genes)       # 17224 > 16877


CNV %>% apply(1, mean) %>% range %>% round(3)     # 1.686 5.469
CNV %>% apply(2, mean) %>% range %>% round(3)     # 1.627 6.302
CNV_C %>% apply(1, mean) %>% range %>% round(3)   # 1.512 2.420
CNV_C %>% apply(2, mean) %>% range %>% round(3)   # 1.228 2.547

breaks1 = seq(1.5, 5.5, 0.25)
breaks2 = seq(1.5, 6.5, 0.25)
breaks3 = seq(1.2, 2.6, 0.1)
breaks4 = seq(1.2, 2.6, 0.1)

dir = mkdir("~/GDSC_Last/Intermediate_1_Cell/CNV")
file1 = sprintf("%s/CNV [Mean per Cell]", dir)
file2 = sprintf("%s/CNV [Mean per Gene]", dir)
file3 = sprintf("%s/CNV Category [Mean per Cell]", dir)
file4 = sprintf("%s/CNV Category [Mean per Gene]", dir)

CNV %>% apply(1, mean) %>% hist_def(file1, breaks=breaks1, show_info=T, show_title=F)
CNV %>% apply(2, mean) %>% hist_def(file2, breaks=breaks2, show_info=T, show_title=F)
CNV %>% apply(1, mean) %>% hist_def(file1, breaks=breaks1, show_info=T, show_title=T)
CNV %>% apply(2, mean) %>% hist_def(file2, breaks=breaks2, show_info=T, show_title=T)

CNV_C %>% apply(1, mean) %>% hist_def(file3, breaks=breaks3, show_info=T, show_title=F)
CNV_C %>% apply(2, mean) %>% hist_def(file4, breaks=breaks4, show_info=T, show_title=F)
CNV_C %>% apply(1, mean) %>% hist_def(file3, breaks=breaks3, show_info=T, show_title=T)
CNV_C %>% apply(2, mean) %>% hist_def(file4, breaks=breaks4, show_info=T, show_title=T)

file = "CNV.RData"
write.csv(CNV, sprintf("%s/CNV.csv", dir), quote=F)
write.csv(CNV_Ent, sprintf("%s/CNV_Ent.csv", dir), quote=F)
write.csv(CNV_C, sprintf("%s/CNV_Category.csv", dir), quote=F)
write.csv(CNV_C_Ent, sprintf("%s/CNV_Category_Ent.csv", dir), quote=F)
save(CNV, CNV_Ent, CNV_C, CNV_C_Ent, file=sprintf("%s/%s", dir, file))

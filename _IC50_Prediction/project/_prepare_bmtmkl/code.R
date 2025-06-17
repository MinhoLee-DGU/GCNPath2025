#!/usr/bin/env Rscript

##### 1. Packages and Data

suppressMessages(library(cogena))
suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))
suppressMessages(library(impute))

source("../functions.R")
loadings()

dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/Anno_Cells.csv", dir)
Anno_Cells = read.csv(file)

# Geneset [C2, CP]
dir = "../../raw_data/MSigDB"

file = sprintf("%s/c2.all.v2024.1.Hs.symbols.gmt", dir)
Path_C2 = gmt2list(file)

file = sprintf("%s/c2.cp.v2024.1.Hs.symbols.gmt", dir)
Path_CP = gmt2list(file)




##### 2. Pre-process Data

# 1. Proteome
dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/Protein_matrix_averaged_20250211.tsv", dir)
Proteome = fread(file, sep="\t")   # 899 x 144

cells = Proteome$V2[4:nrow(Proteome)]
genes = Proteome[2, 3:ncol(Proteome)] %>% as.character
Proteome = Proteome[4:nrow(Proteome), 3:ncol(Proteome)] %>% as.data.frame   # 948 x 8453

rownames(Proteome) = cells
colnames(Proteome) = genes

Proteome[Proteome==""] = 0
Proteome[is.na(Proteome)] = 0
Proteome[] = sapply(Proteome, as.numeric)


# 2. Expression
dir = "../../processed_data/cell_data/GDSC"
file = sprintf("%s/RNA_Array_Sym.csv", dir)
RNA_Array = fread_def(file)   # 1006 x 17419


# 3. Methylation
dir = "../../raw_data/GDSC"
file = sprintf("%s/F2_METH_CELL_DATA.txt", dir)
Methylation = fread_def(file, sep="\t")   # 14726 x 1029
Methylation = Methylation %>% t %>% as.data.frame   # 14726 x 1029

file = sprintf("%s/methSampleId_2_cosmicIds.xlsx", dir)
Anno_Meth = read.xlsx(file)   # 1029 x 20

ids_meth = paste(Anno_Meth$Sentrix_ID, Anno_Meth$Sentrix_Position, sep="_")
idx = match(rownames(Methylation), ids_meth)
cells = Anno_Meth$Sample_Name[idx]

cells %>% is.na %>% sum       # 0
cells %>% unique %>% length   # 1028 [1029]

idx = duplicated(cells)
Methylation = Methylation[!idx, ]
rownames(Methylation) = cells[!idx]

idx = match(rownames(Methylation), Anno_Cells$MODEL_NAME)
cells = Anno_Cells$SANGER_MODEL_ID[idx]
cells %>% is.na %>% sum   # 108

Methylation = Methylation[!is.na(cells), ]   # 920 x 14726
rownames(Methylation) = na.omit(cells)


# 4. RNA-seq
dir = "../../processed_data/cell_data/SANGER_Passports"
file = sprintf("%s/TPM_Sym.csv", dir)
RNA_Seq = fread_def(file)   # 1431 x 37600


# 5. CNV
dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/WES_pureCN_CNV_genes_20221213.csv", dir)
CNV = fread(file)   # 24488104 x 24

CNV$gene_mean %>% na.omit %>% range   # -31.219798   7.426159
CNV$gene_mean %>% na.omit %>% mean    # -0.0324042

mean_ = function(x) mean(x, na.rm=T)
CNV = acast(CNV, model_id~symbol, value.var="gene_mean", fun.aggregate=mean_, fill=0)   # 1253 x 18348
CNV[is.na(CNV)] = 0


# 6. Exome-seq
dir = "../../raw_data/SANGER_Passports"
file = sprintf("%s/mutations_all_20230202.csv", dir)
MUT = fread(file)   # 10050692 x 13

any_ = function(x) as.integer(any(!is.na(x)))
MUT = MUT %>% reshape2::acast(model_id~gene_symbol, value.var="gene_symbol",
                              fun.aggregate = function(x) as.integer(length(x)>0), fill=0)

MUT[is.na(MUT)] = 0   # 1435 x 23189




##### 3. Process Data

# 0. Intersecting cell-lines...
cells = Reduce(intersect, list(rownames(Proteome), rownames(RNA_Array), 
                               rownames(Methylation), rownames(RNA_Seq), 
                               rownames(CNV), rownames(MUT)))   # 808

dir = mkdir("../../benchmark_test/BMTMKL/_data")
file = sprintf("%s/cell_names.txt", dir)
write.table(cells, file=file, row.names=F, col.names=F, quote=F)


# 1. Real-valued genomic views
# Gaussian Kernel

calc_sim <- function(x, y, sigma = 1) {
  sigma <- length(x)
  dist <- sum((x - y)^2)
  sim <- exp(-dist / (2 * sigma))
  return(sim)
}

calc_sim_mat <- function(X) {
  # dist calculates row-wise Euclidean distances
  sigma <- ncol(X)
  dist <- as.matrix(dist(X))
  sim <- exp(-dist^2 / (2 * sigma))
  return(sim)
}

# # Toy Examples
# X = matrix(1:9, nrow=3)
# calc_sim_mat(X)
# calc_sim(as.numeric(X[1, ]), as.numeric(X[2, ]))   # Validation


# 2. Binary-valued genomic views
# Jaccard Kernel

calc_sim_mat_jaccard <- function(X) {
  # dist calculates row-wise Jaccard distances
  dist <- as.matrix(dist(X, method="binary"))
  sim <- 1 - dist
  return(sim)
}


### 0. Original View

Prot_View = Proteome[cells, ] %>% scale %>% calc_sim_mat   # 948 x 948
RNA_Array_View = RNA_Array[cells, ] %>% scale %>% calc_sim_mat
RNA_Seq_View = RNA_Seq[cells, ] %>% scale %>% calc_sim_mat
Meth_View = Methylation[cells, ] %>% scale %>% calc_sim_mat
CNV_View = CNV[cells, ] %>% scale %>% calc_sim_mat
MUT_View = MUT[cells, ] %>% calc_sim_mat_jaccard


### 1. Geneset View

aggregate_genesets <- function(rna_matrix, genesets, method = c("mean", "max")) {
  method <- match.arg(method)
  
  result <- matrix(NA, nrow = nrow(rna_matrix), ncol = length(genesets))
  rownames(result) <- rownames(rna_matrix)
  colnames(result) <- names(genesets)
  
  i_na = c()
  for (i in seq_along(genesets)) {
    genes <- genesets[[i]]
    genes_in_data <- intersect(genes, colnames(rna_matrix))
    
    if (length(genes_in_data) == 0) {
      # warning(paste("No matching genes found for gene set:", names(genesets)[i]))
      i_na = i_na %>% c(i)
      next
    }
    
    submatrix <- rna_matrix[, genes_in_data, drop = FALSE]
    
    if (method == "mean") {
      result[, i] <- rowMeans(submatrix, na.rm = TRUE)
    } else if (method == "max") {
      result[, i] <- apply(submatrix, 1, max, na.rm = TRUE)
    }
  }
  
  if (length(i_na)!=0) {
    sprintf("# The number of genesets without matching genes : %s", length(i_na)) %>% print
    result = result[, -i_na]
  }
  
  return(result)
}

# Methylation data were generated by loci, not by genes
# Methylation data could not be utilized in geneset views

EXP_C2_View = RNA_Array[cells, ] %>% aggregate_genesets(Path_C2, method="mean")
EXP_C2_View = EXP_C2_View %>% scale %>% calc_sim_mat

EXP_CP_View = RNA_Array[cells, ] %>% aggregate_genesets(Path_CP, method="mean")
EXP_CP_View = EXP_CP_View %>% scale %>% calc_sim_mat

CNV_C2_View = CNV[cells, ] %>% aggregate_genesets(Path_C2, method="max")
CNV_C2_View = CNV_C2_View %>% scale %>% calc_sim_mat

CNV_CP_View = CNV[cells, ] %>% aggregate_genesets(Path_CP, method="max")
CNV_CP_View = CNV_CP_View %>% scale %>% calc_sim_mat

MUT_C2_View = MUT[cells, ] %>% aggregate_genesets(Path_C2, method="max")
MUT_C2_View = MUT_C2_View %>% calc_sim_mat_jaccard

MUT_CP_View = MUT[cells, ] %>% aggregate_genesets(Path_CP, method="max")
MUT_CP_View = MUT_CP_View %>% calc_sim_mat_jaccard


### 2. Combination View [PARADIGM & gene-wise multiplication]
# The following sites did not work...
# https://sbenz.github.com/Paradigm
# https://sbenz.github.io/Paradigm
# https://github.com/sbenz/Paradigm
# http://paradigm.five3genomics.com

# Methylation data were generated by loci, not by genes
# Methylation data could not be utilized in gene-wise multiplication views
genes_expr_copy = intersect(colnames(RNA_Array), colnames(CNV))    # 15449
intersect(colnames(CNV), colnames(Methylation)) %>% length         # 0
intersect(colnames(RNA_Array), colnames(Methylation)) %>% length   # 0

Expr_CNV_View = RNA_Array[cells, genes_expr_copy] * CNV[cells, genes_expr_copy]   # 808 x 15449
Expr_CNV_View = Expr_CNV_View %>% scale %>% calc_sim_mat


### 3. Present-Absent View

# Mutation data are already processed in discrete values [0 or 1]
RNA_Seq_PA_View = RNA_Seq[cells, ]
RNA_Seq_PA_View[RNA_Seq_PA_View!=0] = 1
RNA_Seq_PA_View = RNA_Seq_PA_View %>% calc_sim_mat_jaccard




##### 4. Process Data [Total 14 View]

dir = mkdir("../../benchmark_test/BMTMKL/_data")
file = sprintf("%s/View.RData", dir)
save(Prot_View, RNA_Array_View, RNA_Seq_View, Meth_View, CNV_View, MUT_View, 
     EXP_C2_View, EXP_CP_View, CNV_C2_View, CNV_CP_View, MUT_C2_View, MUT_CP_View,
     Expr_CNV_View, RNA_Seq_PA_View, file=file)

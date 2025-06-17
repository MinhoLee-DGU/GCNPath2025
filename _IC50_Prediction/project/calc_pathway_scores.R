
calc_path_def = function(Omics, Path_List, cell_row=T, method="gsva", 
                         filt_genes=T, cores=T, fill_na=0, 
                         direction=T, center=T, stableGenes=NULL, ...) {
  
  # Method : UCell
  # Input : Omics [Cells x Genes if cell_row=T]
  # Input : Omics [Genes x Cells if cell_row=F]
  # Output : Omics [Cells x Pathways]
  
  gene_list = unique(unlist(Path_List))
  n_before_gset = length(gene_list)
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  
  if (cell_row) {
    n_before = ncol(Omics)
    n_after = sum(colnames(Omics) %in% gene_list)
    if (filt_genes) Omics = Omics[, colnames(Omics) %in% gene_list]
  } else {
    n_before = nrow(Omics)
    n_after = sum(rownames(Omics) %in% gene_list)
    if (filt_genes) Omics = Omics[rownames(Omics) %in% gene_list, ]
  }
  
  int_ratio = round(100 * n_after/n_before_gset, 2)
  print(sprintf("# Omics data contain %s%% pathway genes... [%s / %s]", int_ratio, n_after, n_before_gset))
  if (filt_genes) print(sprintf("# Genes are filtered from omics data... [%s > %s]", n_before, n_after))
  
  method = tolower(method)
  if (cell_row) Omics = as.matrix(t(Omics))
  if (method %in% c("gsva", "ssgsea")) {
    Omics_Path = calc_gsva(Omics, Path_List, method=method, cores=cores, ...)
  } else if (method=="singscore") {
    Omics_Path = calc_singscore(Omics, Path_List, direction=direction, 
                                center=center, stableGenes=stableGenes, ...)
  } else if (method=="aucell") {
    # TODO
    suppressMessages(library(AUCell))
  } else if (method=="ucell") {
    Omics_Path = calc_ucell(Omics, Path_List, cores=cores, ...)
  } else {
    stop("Choose methods : GSVA, ssGSEA, AUCell, UCell, UniPath")
  }
  
  methods_crow = c("ucell")
  if (cell_row & !(method %in% methods_crow)) {
    Omics_Path = as.data.frame(t(Omics_Path))
  }
    
  if (all(names(Path_List) %in% colnames(Omics_Path))) {
    print(sprintf("# All %s pathways are processed...", length(Path_List)))
  } else {
    col_na = setdiff(names(Path_List), colnames(Omics_Path))
    Omics_Path[, col_na] = fill_na
    Omics_Path = Omics_Path[, names(Path_List)]
    print(sprintf("# The following %s pathways are not processed...", length(col_na)))
    print(sprintf("# Those pathways in all samples were filled with %s...", fill_na))
  }
  return(Omics_Path)
}

calc_gsva = function(Omics, Path_List, method="gsva", cores=T, ...) {
  # Input : Omics [Genes x Cells]
  # Output : Omics [Pathways x Cells]
  
  suppressMessages(library(GSVA))
  if (method=="gsva") {
    param = gsvaParam(Omics, Path_List, ...)
  } else param = ssgseaParam(Omics, Path_List, ...)
  
  BPPARAM = BiocParallel::MulticoreParam(workers=cores)
  Omics_Path = gsva(param, BPPARAM=BPPARAM)
  return(Omics_Path)
}

calc_singscore = function(Omics, Path_List, direction=T, center=T, stableGenes=NULL, ...) {
  # Input : Omics [Genes x Cells]
  # Output : Omics [Pathways x Cells]
  
  suppressMessages(library(GSEABase))
  suppressMessages(library(singscore))
  
  path_names = names(Path_List)
  Path_List = Path_List %>% lapply(GeneSet)
  Omics = Omics %>% rankGenes(stableGenes=stableGenes)
  
  for (i in 1:length(path_names)) Path_List[[i]]@setName = path_names[[i]]
  Omics_Path = multiScore(Omics, upSetColc=Path_List, knownDirection=direction, centerScore=center, ...)
  Omics_Path = as.data.frame(Omics_Path[[1]])
  return(Omics_Path)
}

calc_ucell = function(Omics, Path_List, cores=T, ...) {
  # Input : Omics [Genes x Cells]
  # Output : Omics [Cells x Pathways]
  
  suppressMessages(library(UCell))
  ranks = StoreRankings_UCell(Omics, ncores=cores)
  Omics_Path = ScoreSignatures_UCell(features=Path_List, precalc.ranks=ranks, ncores=cores, ...)
  
  colnames(Omics_Path) = gsub("_UCell", "", colnames(Omics_Path))
  Omics_Path = as.data.frame(Omics_Path)
  return(Omics_Path)
}

unify_genes = function(Omics, genes, value=0, cell_row=T) {
  if (!is.data.frame(Omics)) Omics = as.data.frame(Omics)
  genes_omics = if (cell_row) colnames(Omics) else rownames(Omics)
  genes_na = genes[!(genes %in% genes_omics)]
  
  if (length(genes_na)==length(genes)) {
    print("No gene matching, check the categories of gene names")
  } else {
    print(sprintf("The number of genes missing : %s", len(genes_na)))
  }
  
  if (cell_row) {
    Omics[, genes_na] = value
    Omics = Omics[, genes]
  } else {
    Omics[genes_na, ] = value
    Omics = Omics[genes, ]
  }
  return(Omics)
}

combat_def = function(Omics_List, tcga_list, db_list, 
                      ref_batch=1, par_prior=T, cell_row=T, cores=T, ...) {
  
  # Omics_List [List of data.frame] : Omics_Source, Omics_Target1, Omics_Target2, ...
  # tcga_list [List of character] : tcga_source, tcga_target1, tcga_target2, ...
  # db_list [Character] : db_source, db_target1, db_target2, ... 
  
  suppressMessages(library(sva))
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  BPPARAM = BiocParallel::MulticoreParam(workers=cores)
  
  t_def = function(df) as.data.frame(t(df))
  if (cell_row) Omics_List = lapply(Omics_List, t_def)
  cell_list = lapply(Omics_List, colnames)
  
  n_data = sapply(Omics_List, ncol)
  Omics_List = Reduce(cbind, Omics_List)
  
  if (is.null(tcga_list)) {
    Mod = NULL
  } else {
    Pheno = data.frame(Tissue=unlist(tcga_list))
    rownames(Pheno) = colnames(Omics_List)
    Mod = model.matrix(~as.factor(Tissue), data=Pheno)
  }
  
  # Reference-batch ComBat with covariates
  batch = rep(1:length(db_list), n_data)
  Omics_List = ComBat(dat=Omics_List, batch=batch, mod=Mod, 
                      par.prior=par_prior, ref.batch=ref_batch, BPPARAM=BPPARAM, ...)
  
  Omics_List = t_def(Omics_List)
  Omics_List = lapply(cell_list, function(x) Omics_List[x, ])
  return(Omics_List)
}

pca_exp = function(EXP, db=NULL, tcga_code=NULL, db_lvl=NULL, rank=10) {
  
  sd_zero = (sapply(EXP, sd)==0)
  if (any(sd_zero)) sprintf("# Num of SD=0 : %s", sum(sd_zero)) %>% print
  
  EXP_PCA = EXP[, !sd_zero] %>% prcomp(center=T, scale=T, rank.=rank)
  EXP_PCA = EXP_PCA$x %>% as.data.frame
  
  if (!is.null(tcga_code)) {
    EXP_PCA = EXP_PCA %>% mutate(TCGA_Code=tcga_code)
  }
  
  if (!is.null(db)) {
    EXP_PCA = EXP_PCA %>% mutate(Database=db) %>% 
      relocate(Database, .before=everything())
  }
  
  if (!is.null(db_lvl)) {
    EXP_PCA$Database = EXP_PCA$Database %>% factor(levels=db_lvl)
  }
  
  return(EXP_PCA)
}

plot_pca_batch = function(Omics_List, db_list=NULL, tcga_list=NULL, tcga_color=NULL, 
                          main=NULL, add_db=NULL, add_tcga=NULL, legend="Tissue", 
                          size=0.5, width=16.5, width2=16.5, height=13.5, height2=13.5, 
                          save=T, save_svg=F, return_pca=T, ...) {
  
  # Omics_List [List of data.frame] : Omics_Source, Omics_Target1, Omics_Target2, ...
  # tcga_list [List of character] : tcga_source, tcga_target1, tcga_target2, ...
  # db_list [Character] : db_source, db_target1, db_target2, ... 
  
  # Annotation Label [TCGA Code]
  if (!is.null(tcga_list)) {
    tcga_code = tcga_list %>% unlist
    cond = identical(sapply(tcga_list, length), sapply(Omics_List, nrow))
    if (!cond) stop("Check the TCGA Codes...")
  } else {
    tcga_code = NULL
  }
  
  # Annotation Label [DB]
  n_samples = sapply(Omics_List, nrow)
  db = rep(db_list, n_samples)
  
  # Compress Omics Data by PCA
  col = Reduce(intersect, lapply(Omics_List, colnames))
  rbind_ = function(df1, df2) rbind(df1[, col], df2[, col])
  Omics_PCA = Reduce(rbind_, Omics_List) %>% 
    pca_exp(db=db, tcga_code=tcga_code, db_lvl=db_list)
  
  # Set X and Y Ranges
  xlim = Omics_PCA$PC1 %>% range
  ylim = Omics_PCA$PC2 %>% range
  xlim = c(xlim[1]-0.1, xlim[2]+0.1)
  ylim = c(ylim[1]-0.1, ylim[2]+0.1)
  
  # Plot PCA [Database]
  Omics_PCA %>% arrange(Database) %>% 
    plot_def(PC1, PC2, color=Database, shape=Database, main=main, add=add_db,
             legend="Database", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
             size=size, width=width, height=height, save=save, save_svg=save_svg, ...)
  
  if (!is.null(tcga_color)) {
    for (i in 1:length(db_list)) {
      main_ = sprintf("%s (%s colored by %s)", main, db_list[i], legend)
      Omics_PCA %>% subset(Database==db_list[i] & TCGA_Code %in% tcga_color) %>% 
        plot_def(PC1, PC2, color=TCGA_Code, shape=TCGA_Code, main=main_, add=add_tcga,
                 legend=legend, xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
                 size=size, width=width2, height=height2, save=save, save_svg=save_svg, ...)
    }
  }
  
  if (return_pca) return(Omics_PCA)
}

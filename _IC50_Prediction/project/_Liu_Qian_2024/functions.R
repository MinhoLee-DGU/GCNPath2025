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
                          main=NULL, add_db=NULL, add_tcga=NULL, legend="TCGA Code", 
                          size=0.5, width=16.5, width2=16.5, height=13.5, height2=13.5, 
                          save=T, save_svg=F, return_pca=T, ...) {
  
  # Omics_List [List of data.frame] : Omics_Source, Omics_Target1, Omics_Target2, ...
  # tcga_list [List of character] : tcga_source, tcga_target1, tcga_target2, ...
  # db_list [Character] : db_source, db_target1, db_target2, ... 
  
  # Annotation Label [TCGA Code]
  if (!is.null(tcga_list)) {
    tcga_code = tcga_list %>% unlist
    cond = identical(unname(sapply(tcga_list, length)), unname(sapply(Omics_List, nrow)))
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
      add_tcga = c(add_tcga, list(labs(shape=legend)))
      main_ = sprintf("%s (%s colored by %s)", main, db_list[i], legend)
      Omics_PCA %>% subset(Database==db_list[i] & TCGA_Code %in% tcga_color) %>% 
        plot_def(PC1, PC2, color=TCGA_Code, shape=TCGA_Code, main=main_, add=add_tcga,
                 legend=legend, xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2", 
                 size=size, width=width2, height=height2, save=save, save_svg=save_svg, ...)
    }
  }
  
  if (return_pca) return(Omics_PCA)
}

plot_pca_tcga = function(Omics_List, Omics_List_CB, db_list=NULL, 
                         tcga_list=NULL, main=NULL, idx_tcga=NULL, 
                         width=18.6, width2=18.6, height=13.5) {
  
  if (!is.null(idx_tcga)) {
    Omics_List[[2]] = Omics_List[[2]][idx_tcga, ]
    Omics_List_CB[[2]] = Omics_List_CB[[2]][idx_tcga, ]
    if (!is.null(tcga_list)) tcga_list[[2]] = tcga_list[[2]][idx_tcga]
  }
  
  if (!is.null(tcga_list)) {
    tissue_top7 = tcga_list[[2]] %>% 
      table %>% sort(decreasing=T) %>% head(7) %>% names
  } else tissue_top7 = NULL
  
  
  shape_db = c(21, 24)
  color_db = c("brown1", "royalblue1")
  names(color_db) = db_list
  
  shape_tcga = 0:6
  color_tcga = c("brown1", "coral", "gold", "seagreen3", "royalblue1", "mediumorchid", "black")
  if (!is.null(tissue_top7) && length(tissue_top7)<7) color_tcga = color_tcga[1:length(color_tcga)]
  
  add_db = list(scale_color_manual(values=color_db), 
                scale_shape_manual(values=shape_db), 
                theme(legend.key.size=unit(1, "cm")), 
                guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
  
  add_tcga = list(scale_color_manual(values=color_tcga), 
                  scale_shape_manual(values=shape_tcga), 
                  theme(legend.key.size=unit(1, "cm")), 
                  guides(shape=guide_legend(override.aes=list(size=3, alpha=1))))
  
  PCA = Omics_List %>% plot_pca_batch(
    db_list, tcga_list, tcga_color=tissue_top7,
    main=main[1], size=1.5, alpha=0.25, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, width2=width2, height=height,
    add_db=add_db, add_tcga=add_tcga, save=T, save_svg=T)
  
  PCA_CB = Omics_List_CB %>% plot_pca_batch(
    db_list, tcga_list, tcga_color=tissue_top7, 
    main=main[2], size=1.5, alpha=0.25, legend_tl=20, legend_tx=18,
    axis_tl=24, axis_tx=18, width=width, width2=width2, height=height, 
    add_db=add_db, add_tcga=add_tcga, save=T, save_svg=T)
  
  PCA_List = list(PCA, PCA_CB)
  names(PCA_List) = c("ComBat_X", "ComBat_O")
  return(PCA_List)
}
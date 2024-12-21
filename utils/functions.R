#!/usr/bin/env Rscript

len = length

not = function(x) "!"(x)

ifelse_def = function(cond, x, y) {
  ifelse(cond, return(x), return(y))
}

floor_def = function(x, digit=0) {
  floor(x * 10**digit) / 10**digit
}

ceil_def = function(x, digit=0) {
  ceiling(x * 10**digit) / 10**digit
}

append_def = function(A, ...) {
  add_list = list(...)
  A = c(A, add_list)
  return(A)
}

standard_cols = function(x, symbols=c("\\.", "-", " ", "\\(", "\\)", "#"), specific=NULL) {
  
  if (!is.null(symbols)) {
    symbols = paste0(symbols, collapse="|")
    x = gsub(symbols, "_", x)
  }
  
  toupper_first = function(a) sub("(.)", "\\U\\1", a, perl=T)
  x = gsub("[_]{2, }", "_", x) %>% strsplit("_") %>% sapply(toupper_first)
  
  if (!is.null(specific)) {
    replace_words = function(a, specific) {
      a_lower = a %>% tolower
      s_lower = specific %>% tolower
      idx = match(a_lower, s_lower, incomparables=NA)
      a[!is.na(idx)] = specific[na.omit(idx)]
      return(a)
    }
    x = x %>% sapply(function(a) replace_words(a, specific))
  }
  
  x = x %>% sapply(function(a) paste0(a[a!=""], collapse="_"))
  return(x)
}

ulen = function(x, len=T) {
  if (len) {
    length(unique(na.omit(x)))
  } else unique(na.omit(x))
}

inlen = function(..., len=T) {
  args = list(...)
  if (len) {
    length(Reduce(intersect, args))
  } else Reduce(intersect, args)
}

union_len = function(..., len=T) {
  args = list(...)
  if (len) {
    length(Reduce(union, args))
  } else Reduce(union, args)
}

na = function(x, not=F) {
  if (!not) {
    sum(is.na(x))
  } else sum(!is.na(x))
}

object = function(x, envir=NULL) {
  if (is.null(envir)) envir=parent.frame()
  eval(parse(text=x), envir=envir)
}

in_rank = function(x, n, top=T) {
  if (n<0) stop("Parameter prob not adequete...")
  if (n>=length(x)) n = length(x)
  if (n>=0 & n<1) n = floor(length(x)*n) 
  cond = ifelse_def(top, x>=sort(x, decreasing=T)[n], x<=sort(x)[n])
  if (n==0) cond = rep(F, length(x))
  return(cond)
}

in_top = function(x, n) {
  in_rank(x, n, top=T)
}

in_bottom = function(x, n) {
  in_rank(x, n, top=F)
}

dup = function(df, col) {
  col = deparse(substitute(col))
  df_sub = df %>% group_by_at(col) %>% 
    filter(n()>1) %>% as.data.frame
  return(df_sub)
}

go = function(dir) {
  if (is.character(dir)) {
    dir = strsplit(dir, "/") %>% unlist
    for (i in 1:length(dir)) {
      if (!(dir[i] %in% list.files()) & !(dir[i]=="..")) {
        dir.create(dir[i])
      }
      setwd(dir[i])
    }
  } else if (is.numeric(dir)) {
    for (i in 1:dir) setwd("..")
  } else print("Error!!!")
}

mkdir = function(dir) {
  dir_ori = getwd()
  suppressWarnings(go(dir))
  setwd(dir_ori)
  return(dir)
}

h = function(x, what=6) {
  if (is.vector(x) | is.factor(x) | is.character(x)) {
    # numeric, integer, character, list
    range_ = min(what, length(x))
    x[1:range_]
  } else if (is.table(x)) {
    # table
    range_ = min(what, length(x))
    x[1:range_]
  } else {
    # matrix, data.frame
    range_row = min(what, nrow(x))
    range_col = min(what, ncol(x))
    x[1:range_row, 1:range_col]
  }
}

sh = h

fread_def = function(file, check_names=F, col_numeric=F, ...) {
  # fread the file into data.frame
  df = fread(file, check.names=check_names, ...)
  df = df %>% data.frame(row.names=df[[1]], check.names=check_names)
  df = df[, -1]
  
  if (col_numeric) {
    colnames(df) = unlist(df[1, ])
    df = df[-1, ]
  }
  return(df)
}

loadings = function(packages=NULL) {
  packages_base = c("dplyr", "ggplot2", "caret", "tidyr", "data.table")
  packages = c(packages, packages_base)
  invisible(lapply(packages, library, quietly=T, warn.conflict=F, character.only=T))
}

grep_def = function(x, words, and=T, not=F) {
  
  if (not) {
    # Not include words
    x = words %>% sapply(function(word) x[-grep(word, x)])
  } else {
    # MUST include words
    x = words %>% sapply(function(word) grep(word, x, value=T))
  }
  
  if (and) {
    # Logic AND
    x = Reduce(intersect, x)
  } else {
    # Logic OR
    x = Reduce(union, x)
  }
  
  return(x)
}

text_repel_def = function(df, label, label_row=F, rectangle=F,
                          direction="both", size=4, force=1, overlap=10, ...) {
  
  suppressMessages(library(ggrepel))
  
  if (!label_row) {
    mapping = aes(label={{label}})
  } else {
    df$Labels = rownames(df)
    mapping = aes(label={{Labels}})
  }
  
  text_func = ifelse_def(rectangle, geom_label_repel, geom_text_repel)
  text = text_func(data=df, mapping=mapping, direction=direction, 
                   size=size, force=force, max.overlaps=overlap, ...)
  
  return(text)
}

set_mapping = function(geom_object, df, var=NULL, category=NULL) {
  
  if (!is.null(var) & !identical(var, "NULL")) {
    if (length(vars)!=1 | length(category)!=1) {
      stop("Category and its variable is sould be 1:1...")
    }
    # var_contained = any(grepl(paste0(colnames(df), collapse="|"), var))
    # geom_object$mapping[[category]] = quo(sym(var))
    
    if (var %in% colnames(df)) {
      geom_object$mapping[[category]] = sym(var)
    } else {
      geom_object$aes_params[[category]] = var
    }
  }
  return(geom_object)
}

set_margin_lg = function(margin_lg, pos_legend=NULL, ratio_title=2.5) {
  # Assume margin_lg is given as numeric
  if (is.null(pos_legend) || pos_legend=="right") {
    margin_lgl = margin(b=margin_lg, unit="cm")
    margin_lgx = margin(b=margin_lg, l=margin_lg, unit="cm")
  } else if (pos_legend=="bottom") {
    margin_lgl = margin(r=margin_lg*ratio_title, unit="cm")
    margin_lgx = margin(r=margin_lg, l=margin_lg, unit="cm")
  } else {
    stop("This function supports legend position [bottom, right]...")
  }
  return(list(margin_lgl, margin_lgx))
}

save_fig = function(pl, main, width, height, dpi=400, svg=F, ggplot=T, ...) {
  
  main = gsub("\n", " ", main)
  file_svg = sprintf("%s.svg", main)
  file_png = sprintf("%s.png", main)
  
  if (ggplot) {
    ggsave(file_png, pl, width=width, height=height, dpi=dpi,...)
  } else {
    png(file_png, width=width, height=height, units="cm", res=dpi)
    grid::grid.newpage()
    grid::grid.draw(pl$gtable)
    dev.off()
  }
  
  if (svg) {
    if (ggplot) {
      ggsave(file_svg, pl, width=width/2.54, height=height/2.54, 
             dpi=dpi, device=grDevices::svg, family="sans")
    } else {
      grDevices::svg(file_svg, width=width/2.54, height=height/2.54, family="sans")
      grid::grid.newpage()
      grid::grid.draw(pl$gtable)
      dev.off()
    }
  }
}

save_fig_ggpubr = function(pl, main, width=20, height=15, dpi=1200, svg=F, ...) {
  
  size = c(width, height)
  size_in = size / 2.54
  size_px = size_in * dpi
  
  main = gsub("\n", " ", main)
  file_svg = sprintf("%s.svg", main)
  file_png = sprintf("%s.png", main)
  ggexport(pl, filename=file_png, width=size_px[1], height=size_px[2], res=dpi, ...) %>% suppressMessages
  
  if (svg) {
    ggexport(pl, filename=file_svg, width=size_in[1], height=size_in[2]) %>% suppressMessages
  }
}

hist_def = function(df, x=NULL, main=NULL, breaks=NULL, scale_x=NULL, 
                    xlab=NULL, ylab=NULL, text=NULL, add=NULL, legend=F, pos_legend=NULL,
                    plot_tl=22.5, axis_tl=18, axis_tx=15, text_tx=15, legend_tl=10, legend_tx=10, 
                    margin=0.5, margin_pl=0.25, margin_lg=0.4, text_round=3, dist=1, 
                    alpha=1, alpha_dist=0.25, width=16, height=13.5, dpi=400, text_ratio=1, size_line=1,
                    plot_face="plain", axis_face="plain", legend_face="plain", 
                    color="black", fill="grey", color_dist="black", fill_dist="red", 
                    text_info=T, hist=T, vline=F, show_title=F, force_bold=F, save=F, save_svg=T) {
  
  title = NULL
  e_text = element_text
  stopifnot(is.numeric(dist) | hist)
  filt_dash = function(x) gsub("\\\"", "", x)
  
  x_string = deparse(substitute(x)) %>% filt_dash
  fill_string = deparse(substitute(fill)) %>% filt_dash
  color_string = deparse(substitute(color)) %>% filt_dash
  fill_string_d = deparse(substitute(fill_dist)) %>% filt_dash
  color_string_d = deparse(substitute(color_dist)) %>% filt_dash
  
  if (is.data.frame(df)) {
    df$Var_Temp = df[[x_string]]
  } else {
    df = df %>% na.omit %>% as.data.frame %>% setNames("Var_Temp")
  }
  
  mean_x = df$Var_Temp %>% mean(na.omit=T) %>% round(text_round)
  median_x = df$Var_Temp %>% median(na.omit=T) %>% round(text_round)
  if (text_info) text = sprintf("Mean & Median : %s & %s\n", mean_x, median_x)
  
  if (is.null(breaks)) {
    binwidth = (max(df$Var_Temp) - min(df$Var_Temp)) / 30
    breaks = seq(min(df$Var_Temp), max(df$Var_Temp), binwidth)
  } else {
    binwidth = breaks[2] - breaks[1]
  }
  
  if (show_title) { 
    title = sprintf("%s\n", main)
    title = title %>% strsplit("/") %>% unlist
    title = title[length(title)]
  }
  
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  margin_x = margin(t=margin, b=margin, unit="cm")
  margin_y = margin(l=margin, r=margin, unit="cm")
  
  if (is.numeric(margin_lg)) margin_lg = set_margin_lg(margin_lg, pos_legend, ratio_title=2.5)
  margin_lgl = margin_lg[[1]]
  margin_lgx = margin_lg[[-1]]
  
  if (force_bold) {
    plot_face = "bold"; axis_face = "bold"; legend_face = "bold";
  }
  
  plot_tl = plot_tl * text_ratio
  axis_tl = axis_tl * text_ratio
  axis_tx = axis_tx * text_ratio
  legend_tl = legend_tl * text_ratio
  legend_tx = legend_tx * text_ratio
  
  pl = ggplot(df, aes(Var_Temp)) + 
    theme_classic() + labs(x=xlab, y=ylab, title=title, caption=text) + 
    theme(plot.margin = margin_pl,
          plot.title = e_text(size=plot_tl, face=plot_face, hjust=0.5),
          axis.title = e_text(size=axis_tl, face=axis_face, hjust=0.5), 
          axis.text.x = e_text(size=axis_tx, face=axis_face, margin=margin_x), 
          axis.text.y = e_text(size=axis_tx, face=axis_face, margin=margin_y), 
          legend.title = e_text(size=legend_tl, face=legend_face, margin=margin_lgl),
          legend.text = e_text(size=legend_tx, face=legend_face, margin=margin_lgx),
          plot.caption = e_text(size=text_tx)) + xlim(range(breaks))
  
  if (hist) {
    pl_hist = geom_histogram(breaks=breaks, alpha=alpha)
    pl_hist = pl_hist %>% set_mapping(df, fill_string, "fill")
    pl_hist = pl_hist %>% set_mapping(df, color_string, "colour")
    pl = pl + pl_hist
  }
  
  if (is.numeric(dist)) {
    pl_dist = geom_density(aes(y=after_stat(density)*nrow(df)*binwidth*dist), alpha=alpha_dist)
    pl_dist = pl_dist %>% set_mapping(df, fill_string_d, "fill")
    pl_dist = pl_dist %>% set_mapping(df, color_string_d, "colour")
    pl = pl + pl_dist
  }
  
  if (!is.null(scale_x)) pl = pl + scale_x_continuous(breaks=scale_x)
  if (is.numeric(vline)) pl = pl + 
    sapply(vline, function(x) geom_vline(xintercept=x, color="red", linetype=2, linewidth=size_line))
  
  if (legend==F) pl = pl + theme(legend.position = "none")
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  if (!is.null(pos_legend)) pl = pl + theme(legend.position=pos_legend)
  
  if (save) {
    save_fig(pl, main, svg=save_svg, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

plot_def = function(df, x, y, main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, legend=NULL, 
                    text=NULL, add=NULL, color="black", shape=NULL, color_line="red", size=1, alpha=1, stroke=1,
                    plot_tl=22.5, axis_tl=18, axis_tx=15, legend_tl=16.5, legend_tx=15, 
                    margin=0.5, margin_pl=0.25, margin_lg=0.4, width=15, height=15, dpi=400, text_ratio=1,
                    plot_face="plain", axis_face="plain", legend_face="plain", raster_dev="ragg_png", pos="identity", 
                    pos_legend=NULL, show_title=F, xy_line=F, force_bold=F, raster=F, save=F, save_svg=T) {
  
  title = NULL
  e_text = element_text
  filt_dash = function(x) gsub("\\\"", "", x)
  if (raster) suppressMessages(library(ggrastr))
  
  if (show_title) { 
    title = sprintf("%s\n", main)
    title = title %>% strsplit("/") %>% unlist
    title = title[length(title)]
  }
  
  mapping = aes(x={{x}}, y={{y}})
  size = deparse(substitute(size)) %>% filt_dash
  color = deparse(substitute(color)) %>% filt_dash
  shape = deparse(substitute(shape)) %>% filt_dash
  size = tryCatch(as.numeric(size), warning=function(e) return(size))
  shape = tryCatch(as.numeric(shape), warning=function(e) return(shape))
  
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  margin_x = margin(t=margin, b=margin, unit="cm")
  margin_y = margin(l=margin, r=margin, unit="cm")
  
  if (is.numeric(margin_lg)) margin_lg = set_margin_lg(margin_lg, pos_legend, ratio_title=2.5)
  margin_lgl = margin_lg[[1]]
  margin_lgx = margin_lg[[-1]]
  
  if (is.null(legend) || legend==T) legend = color
  xlab = ifelse_def(!is.null(xlab), xlab, deparse(substitute(x)))
  ylab = ifelse_def(!is.null(ylab), ylab, deparse(substitute(y)))
  
  # if (is.null(legend)) legend = color
  # if (legend==T) legend = color
  
  if (force_bold) {
    plot_face = "bold"; axis_face = "bold"; legend_face = "bold" 
  }
  
  plot_tl = plot_tl * text_ratio
  axis_tl = axis_tl * text_ratio
  axis_tx = axis_tx * text_ratio
  legend_tl = legend_tl * text_ratio
  legend_tx = legend_tx * text_ratio
  
  pl_point = geom_point(alpha=alpha, stroke=stroke, position=pos)
  pl_point = pl_point %>% set_mapping(df, size, "size")
  pl_point = pl_point %>% set_mapping(df, color, "colour")
  pl_point = pl_point %>% set_mapping(df, shape, "shape")
  
  if (raster) pl_point = pl_point %>% rasterise(dpi=dpi, dev=raster_dev)
  
  pl = ggplot(df, mapping) + pl_point + 
    theme_classic() + labs(x=xlab, y=ylab, title=title, color=legend) + 
    theme(plot.margin = margin_pl,
          plot.title = e_text(size=plot_tl, face=plot_face, hjust=0.5),
          axis.title = e_text(size=axis_tl, face=axis_face, hjust=0.5), 
          axis.text.x = e_text(size=axis_tx, face=axis_face, margin=margin_x), 
          axis.text.y = e_text(size=axis_tx, face=axis_face, margin=margin_y),
          legend.title = e_text(size=legend_tl, face=legend_face, margin=margin_lgl),
          legend.text = e_text(size=legend_tx, face=legend_face, margin=margin_lgx))
  
  if (!is.null(xlim)) pl = pl + xlim(xlim)
  if (!is.null(ylim)) pl = pl + ylim(ylim)
  if (xy_line) pl = pl + geom_abline(slope=1, intercept=0, color=color_line, lty=2)
  
  if (legend==F) pl = pl + theme(legend.position="none")
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  if (!is.null(pos_legend)) pl = pl + theme(legend.position=pos_legend)
  
  if (save) {
    save_fig(pl, main, svg=save_svg, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

boxplot_def = function(df, x, y, fill=NULL, main=NULL, xlab=NULL, ylab=NULL, 
                       legend=NULL, add=NULL, pos=NULL, pos_point=NULL, pos_legend=NULL, 
                       plot_tl=22.5, axis_tl=15, axis_tx=12, legend_tl=12, legend_tx=12, 
                       angle=30, alpha=0.5, alpha_point=0.5, lwd=0.5, vjust=0.5, hjust=0.5, size_point=1,
                       margin=0.5, margin_pl=0.25, margin_lg=0.4, width=16, height=12, dpi=400, text_ratio=1,
                       plot_face="plain", axis_face="plain", legend_face="plain", raster_dev="ragg_png",
                       point=T, violin=F, trim=F, reorder=F, reorder_median=T, reorder_rev=F,
                       show_title=F, force_bold=F, raster=F, save=F, save_svg=T) {
  
  title = NULL
  e_text = element_text
  if (is.null(pos)) pos = position_dodge(width=0.8)
  if (is.null(pos_point)) pos_point = pos
  if (raster) suppressMessages(library(ggrastr))
  
  x_string = deparse(substitute(x))
  y_string = deparse(substitute(y))
  fill_string = deparse(substitute(fill))
  reorder_string = deparse(substitute(reorder))
  
  if (is.null(legend)) legend = fill_string
  if (legend==T) legend = fill_string
  axis_x = !is.numeric(df[, x_string])
  
  if (reorder) {
    if (axis_x) {
      df[, x_string] = factor(df[, x_string])
    } else df[, y_string] = factor(df[, y_string])
  }
  
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  margin_x = margin(t=margin, b=margin, unit="cm")
  margin_y = margin(l=margin, r=margin, unit="cm")
  
  if (is.numeric(margin_lg)) margin_lg = set_margin_lg(margin_lg, pos_legend, ratio_title=2.5)
  margin_lgl = margin_lg[[1]]
  margin_lgx = margin_lg[[-1]]
  
  reorder_def1 = function(x, ...) stats::reorder(x, ..., FUN=mean)
  reorder_def2 = function(x, ...) stats::reorder(x, ..., FUN=median)
  reorder_def = ifelse_def(!reorder_median, reorder_def1, reorder_def2)
  reorder_def = ifelse_def(!reorder_rev, reorder_def, function(x, ...) rev(reorder_def))
  
  if (fill_string=="NULL" & !reorder) {
    mapping = aes(x={{x}}, y={{y}}, fill={{x}})
    mapping_t = aes(x={{y}}, y={{x}}, fill={{x}})
  } else if (fill_string!="NULL" & !reorder) {
    mapping = aes(x={{x}}, y={{y}}, fill={{fill}})
    mapping_t = aes(x={{y}}, y={{x}}, fill={{fill}})
  } else if (fill_string=="NULL" & reorder) {
    mapping = aes(x=reorder_def({{x}}, {{y}}), y={{y}}, fill=reorder_def({{x}}, {{y}}))
    mapping_t = aes(x=reorder_def({{y}}, {{x}}), y={{x}}, fill=reorder_def({{y}}, {{x}}))
  } else {
    mapping = aes(x=reorder_def({{x}}, {{y}}), y={{y}}, fill=reorder_def({{fill}}, {{x}}))
    mapping_t = aes(x=reorder_def({{y}}, {{x}}), y={{x}}, fill=reorder_def({{fill}}, {{x}}))
  }
  
  if (show_title) {
    title = sprintf("%s\n", main)
    title = title %>% strsplit("/") %>% unlist
    title = title[length(title)]
  }
  
  if (!axis_x) mapping_t = mapping
  pl = ggplot(df, mapping)
  size_outlier = ifelse(!point, size_point, 0)
  
  if (violin) {
    pl_main = geom_violin(alpha=alpha, lwd=lwd, na.rm=T, position=pos, outlier.size=size_outlier, trim=trim)
  } else {
    pl_main = geom_boxplot(alpha=alpha, lwd=lwd, na.rm=T, position=pos, outlier.size=size_outlier)
  }
  
  if (raster) pl_main = pl_main %>% rasterise(dpi=dpi, dev=raster_dev)
  pl = pl + pl_main
  
  if (force_bold) {
    plot_face = "bold"; axis_face = "bold"; legend_face = "bold" 
  }
  
  plot_tl = plot_tl * text_ratio
  axis_tl = axis_tl * text_ratio
  axis_tx = axis_tx * text_ratio
  legend_tl = legend_tl * text_ratio
  legend_tx = legend_tx * text_ratio
  
  if (point) {
    pl_point = geom_point(size=size_point, alpha=alpha_point, position=pos_point, show.legend=F)
    if (raster) pl_point = pl_point %>% rasterise(dpi=dpi, dev=raster_dev)
    pl = pl + pl_point
  }
  
  pl = pl + theme_classic() + labs(x=xlab, y=ylab, title=title, fill=legend) +
    theme(plot.margin = margin_pl,
          plot.title = e_text(size=plot_tl, face=plot_face, hjust=0.5),
          axis.title = e_text(size=axis_tl, face=axis_face, hjust=0.5),
          axis.text.y = e_text(size=axis_tx, face=axis_face, margin=margin_y),
          axis.text.x = e_text(size=axis_tx, face=axis_face, margin=margin_x,
                               angle=angle, vjust=vjust, hjust=hjust), 
          legend.title = e_text(size=legend_tl, face=legend_face, margin=margin_lgl),
          legend.text = e_text(size=legend_tx, face=legend_face, margin=margin_lgx))
  
  if (legend==F) pl = pl + theme(legend.position="none")
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  if (!is.null(pos_legend)) pl = pl + theme(legend.position=pos_legend)
  
  if (save) {
    save_fig(pl, main, svg=save_svg, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

grid_def = function(df, x, y, fill, main=NULL, xlab=NULL, ylab=NULL, color=NULL, color_line="white", 
                    size_line=0.5, legend=NULL, pos_legend=NULL, add=NULL, label=NULL, limits=NULL,
                    round=3, size=6, angle=30, margin=0.5, margin_pl=0.25, margin_lg=0.4, 
                    plot_tl=22.5, axis_tl=18, axis_tx=15, legend_tl=16.5, legend_tx=16.5, 
                    width=15, height=15, dpi=400, text_ratio=1, 
                    plot_face="plain", axis_face="plain", legend_face="plain", text_face="plain", 
                    trans=F, mean_summary=T, show_title=F, force_bold=F, save=F, save_svg=T) {
  
  title = NULL
  e_text = element_text
  
  x_string = deparse(substitute(x))
  y_string = deparse(substitute(y))
  fill_string = deparse(substitute(fill))
  
  df[, x_string] = df[, x_string] %>% as.factor
  df[, y_string] = df[, y_string] %>% as.factor
  if (is.null(legend)) legend = fill_string
  
  if (!is.expression(legend)) {
    if (legend==T) legend = fill_string
  }
  
  if (is.null(label)) {
    if (is.null(round)) {
      label = aes(label=object(fill_string))
    } else label = aes(label=round(object(fill_string), round))
  }
  
  if (is.null(color)) {
    color = cm.colors(25)
    midpoint = (min(df[[fill_string]]) + max(df[[fill_string]])) / 2
    color = scale_fill_gradient2(low=color[1], high=color[length(color)], midpoint=midpoint, limits=limits)
  }
  
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  margin_x = margin(t=margin, b=margin, unit="cm")
  margin_y = margin(l=margin, r=margin, unit="cm")
  
  if (is.numeric(margin_lg)) margin_lg = set_margin_lg(margin_lg, pos_legend, ratio_title=2.5)
  margin_lgl = margin_lg[[1]]
  margin_lgx = margin_lg[[-1]]
  
  if (mean_summary) {
    columns = c(x_string, y_string)
    df = df %>% group_by(across(all_of(columns))) %>% 
      summarise(fill_mean=mean(object(fill))) %>% as.data.frame
    colnames(df) = c(x_string, y_string, fill_string)
  }
  
  if (show_title) {
    title = sprintf("%s\n", main)
    title = title %>% strsplit("/") %>% unlist
    title = title[length(title)]
  }
  
  if (trans) {
    mapping = aes(x={{y}}, y={{x}}, fill={{fill}})
  } else {
    mapping = aes(x={{x}}, y={{y}}, fill={{fill}})
  }
  
  if (force_bold) {
    plot_face = "bold"; axis_face = "bold"; 
    legend_face = "bold" ; text_face = "bold"; 
  }
  
  plot_tl = plot_tl * text_ratio
  axis_tl = axis_tl * text_ratio
  axis_tx = axis_tx * text_ratio
  legend_tl = legend_tl * text_ratio
  legend_tx = legend_tx * text_ratio
  
  pl = ggplot(df, mapping) + 
    geom_tile(color=color_line, linewidth=size_line) + theme_classic() +
    labs(x=xlab, y=ylab, title=title, fill=legend) + color + 
    geom_text(label, size=size, fontface=text_face, check_overlap=T) +
    theme(plot.margin = margin_pl,
          plot.title = e_text(size=plot_tl, face=plot_face, hjust=0.5),
          axis.title = e_text(size=axis_tl, face=axis_face),
          axis.text.y = e_text(size=axis_tx, margin=margin_y, face=axis_face),
          axis.text.x = e_text(size=axis_tx, margin=margin_x, face=axis_face, angle=angle, hjust=1, vjust=1), 
          legend.title = e_text(size=legend_tl, face=legend_face, margin=margin_lgl),
          legend.text = e_text(size=legend_tx, face=legend_face, margin=margin_lgx))
  
  if (legend==F) pl = pl + theme(legend.position = "none")
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  if (!is.null(pos_legend)) pl = pl + theme(legend.position=pos_legend)
  
  if (save) {
    save_fig(pl, main, svg=save_svg, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

# venn_def = function(..., title=NULL, labels=NULL, col=NULL, is_list=F, show_label=T, 
#                     text_size=12, set_name_size=12, width=27, height=27, save=F) {
#   
#   suppressMessages(library(ggvenn))
#   Args = list(...)
#   Args = Args %>% sapply(unique)
#   
#   if (is.null(col)) {
#     if (length(Args)==2) {
#       col = c("red", "blue")
#     } else if (length(Args)==3) {
#       col = c("red", "green", "blue")
#     } else if (length(Args)==4) {
#       col = c("red", "green", "yellow", "blue")
#     } else col=NULL
#   }
#   
#   names(Args) = labels
#   if (show_label) {
#     pl = ggvenn(Args, fill_color=col, fill_alpha=0.3, 
#                 stroke_size=0.8, text_size=text_size, set_name_size=set_name_size)
#   } else {
#     pl = ggvenn(Args, fill_color=col, fill_alpha=0.3, 
#                 stroke_size=0.8, text_size=text_size, set_name_size=0)
#   }
#   
#   if (save) {
#     file = sprintf("%s.png", title)
#     ggsave(file, pl, width=width, height=height, units="cm")
#   } else print(pl)
# }

venn_def = function(..., main=NULL, labels=NULL, col=NULL, add=NULL, 
                    alpha=0.75, label_alpha=0, text_size=6, set_size=6, 
                    margin_pl=0.25, margin_lg=0.4, legend_tl=12, legend_tx=12, 
                    text_ratio=1, width=20, height=16, dpi=400, 
                    legend=T, force_bold=F, save=F, save_svg=T) {
  
  suppressMessages(library(ggVennDiagram))
  e_text = element_text
  if (is.null(legend)) legend = F
  
  Args = list(...)
  if (length(Args)==1) Args = Args[[1]]
  Args = Args %>% lapply(unique)
  
  if (!is.null(labels)) {
    names(Args) = labels
  } else {
    if (is.null(names(Args))) {
      labels = sprintf("Set%s", 1:length(Args))
      names(Args) = labels
    }
  }
  
  if (is.null(col)) {
    col = scale_fill_gradientn(colors=alpha(c("white", "red"), alpha))
  }
  
  set_size = set_size * text_ratio
  text_size = text_size * text_ratio
  legend_tl = legend_tl * text_ratio
  legend_tx = legend_tx * text_ratio
  
  margin_lg = margin(b=margin_lg, unit="cm")
  margin_pl = unit(rep(margin_pl, 4), units="cm")
  
  pl = ggVennDiagram(Args, label_alpha=label_alpha, 
                     label_size=text_size, set_size=set_size)
  
  pl = pl + col + scale_color_manual(values=rep("black", length(Args)))
  if (isFALSE(legend)) pl = pl + theme(legend.position="none")
  if (is.character(legend)) pl = pl + labs(fill=legend)
  if (!is.null(add)) for (i in 1:length(add)) pl = pl + add[[i]]
  
  pl = pl + theme(plot.margin = margin_pl,
                  legend.title = e_text(size=legend_tl, margin=margin_lg),
                  legend.text = e_text(size=legend_tx, margin=margin_lg))
  
  if (save) {
    save_fig(pl, main, svg=save_svg, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

heatmap_def = function(df, main=NULL, Anno_Row=NULL, Anno_Col=NULL, breaks=NULL, 
                       color=NULL, color_type=NULL, color_bin=NULL, 
                       color_center="white", color_num="black", color_border="grey60", 
                       color_scale=1, text_ratio=1, text_row=15, text_col=15, text_num=15, 
                       tree_row=50, tree_col=50, width=20, height=20, dpi=400,  
                       clust="complete", clust_row=T, clust_col=T, show_row=F, show_col=T, 
                       show_num=F, center_zero=F, from_zero=F, 
                       scale_row=F, scale_col=F, save=F, save_svg=T, ...) {
  
  suppressMessages(library(pheatmap))
  if (scale_row & scale_col) stop("Both scale_row=T and scale_col=T not compatible...")
  if (center_zero & from_zero) stop("Both center_zero=T and from_zero=T not compatible...")
  
  if (scale_col) df = df %>% scale %>% as.data.frame
  if (scale_row) df = df %>% apply(1, scale) %>% t %>% as.data.frame
  dist_df = df %>% unlist %>% na.omit %>% as.numeric
  
  norm_zero = function(x) (0-min(x)) / (max(x)-min(x))
  norm_median = function(x) (median(x)-min(x)) / (max(x)-min(x))
  ratio_median = ifelse(!center_zero, norm_median(dist_df), norm_zero(dist_df))
  
  above_zero = from_zero & min(dist_df)>=0
  below_zero = from_zero & max(dist_df)<=0
  
  if (is.null(breaks)) {
    min_df = ifelse(above_zero, 0, min(dist_df))
    max_df = ifelse(below_zero, 0, max(dist_df))
    breaks = seq(min_df, max_df, length=color_scale*100)
  }
  
  if (is.null(color_type)) {
    color_type = c("royalblue2", "white", "firebrick2")
    color_type = ifelse_def(!above_zero, color_type, c("white", "firebrick2"))
    color_type = ifelse_def(!below_zero, color_type, c("royalblue2", "white"))
  }
  
  if (!from_zero) {
    range_low = 100*ratio_median
    range_high = (100-range_low)
    range_low = round(range_low*color_scale)
    range_high = round(range_high*color_scale)
    
    if (is.null(color_bin)) {
      color_bin = c(range_low, range_high)
      sprintf("# Color bin : [%s, %s]", range_low, range_high) %>% print
    } else {
      cond_no_1 = sum(color_bin)!=1
      cond_no_2 = length(color_bin)!=length(color_type)-1
      if (cond_no_1 | cond_no_2) {
        stop("Check the parameter color_bin...")
      } 
      
      if (center_zero) {
        idx_center = which(color_center==color_type)
        color_low = color_bin[1:(idx_center-1)]
        color_high = color_bin[idx_center:length(color_bin)]
        
        norm_ratio = function(x) x/sum(x)
        color_low = norm_ratio(color_low) * range_low
        color_high = norm_ratio(color_high) * range_high
        color_bin = c(color_low, color_high)
      } else {
        color_bin = round(color_bin * 100 * color_scale)
      }
    }
  } else {
    if (is.null(color_bin)) color_bin = 1
    color_bin = round(color_bin * 100 * color_scale)
  }
  
  if (is.null(color)) {
    color = c()
    for (i in 1:(length(color_type)-1)) {
      color_ith = colorRampPalette(c(color_type[i], color_type[i+1]))(color_bin[i])
      color = color %>% c(color_ith)
    }
  }
  
  text_row = text_row * text_ratio
  text_col = text_col * text_ratio
  text_num = text_num * text_ratio
  
  pl = df %>% pheatmap(color=color, breaks=breaks, clustering_method=clust, 
                       cluster_rows=clust_row, cluster_cols=clust_col, 
                       show_rownames=show_row, show_colnames=show_col,
                       annotation_row=Anno_Row, annotation_col=Anno_Col, 
                       fontsize_row=text_row, fontsize_col=text_col, 
                       display_numbers=show_num, fontsize_number=text_num, 
                       number_color=color_num, border_color=color_border, 
                       treeheight_row=tree_row, treeheight_col=tree_col, ...)
  
  if (save) {
    save_fig(pl, main, svg=save_svg, ggplot=F, width=width, height=height, dpi=dpi, units="cm")
  } else print(pl)
}

RMSE_Norm = function(pred, obs, na.rm=T) {
  # Normalized RMSE by interquartile range [NRMSE-IQR]
  rmse = RMSE(pred, obs, na.rm=na.rm)
  iqr_obs = quantile(obs, 3/4) - quantile(obs, 1/4)
  rmse_norm = rmse / iqr_obs
  return(rmse_norm)
}

corr_pair = function(df1, df2, by_row=T, into_wide=F, 
                     method="pearson", use="pairwise.complete.obs") {
  
  if (by_row) df1 = df1 %>% t
  if (by_row) df2 = df2 %>% t
  df_corr = cor(df1, df2, method=method, use=use) %>% as.data.frame
  
  if (!into_wide) {
    df_corr = df_corr %>% as.matrix %>% 
      reshape2::melt() %>% as.data.frame
    colnames(df_corr)[3] = c("Stat")
  }
  
  return(df_corr)
}

stat_pair = function(df1, df2, method=RMSE, by_row=T, into_wide=F, cores=NULL) {
  
  if (is.numeric(cores) | isTRUE(cores)) {
    # if NULL or numeric, foreach mode is activated
    df_stat = stat_pair_foreach(df1, df2, method=method, by_row=by_row, 
                                into_wide=into_wide, cores=cores)
  } else {
    # if F, foreach mode is not activated
    by_row = ifelse(by_row, 1, 2)
    df_stat = apply(df1, by_row, function(x) {
      apply(df2, by_row, function(y) method(x, y))
    })
    
    if (!into_wide) {
      df_stat = df_stat %>% t %>% reshape2::melt() %>% as.data.frame
      colnames(df_stat)[3] = c("Stat")
    }
  }
  
  return(df_stat)
}

stat_pair_foreach = function(df1, df2, method=RMSE, by_row=T, into_wide=F, cores=T) {
  
  ndim_df = ifelse_def(by_row, nrow, ncol)
  names_df = ifelse_def(by_row, rownames, colnames)
  if (isTRUE(cores)) cores = parallel::detectCores()/2
  
  tryCatch({
    cluster = multicores(cores=cores)
    on.exit(stopCluster(cluster))
    
    df_stat = foreach(i=1:ndim_df(df1), .combine=rbind) %dopar% {
      stat_info_ = foreach(j=1:ndim_df(df2), .combine=rbind) %do% {
        v1 = ifelse_def(by_row, as.numeric(df1[i, ]), as.numeric(df1[, i]))
        v2 = ifelse_def(by_row, as.numeric(df2[j, ]), as.numeric(df2[, j]))
        stat_info = c(names_df(df1)[i], names_df(df2)[j], method(v1, v2))
        return(stat_info)
      }
      return(stat_info_)
    }
  }, error = function(e) print("Unknown Error..."))
  
  df_stat = df_stat %>% as.data.frame
  
  rownames(df_stat) = NULL
  colnames(df_stat) = c("Var1", "Var2", "Stat")
  df_stat$Stat = df_stat$Stat %>% as.numeric
  
  if (into_wide) {
    df_stat = reshape2::acast(Var1~Var2, value.var="Stat")
    df_stat = df_stat %>% as.data.frame
  }
  
  return(df_stat)
}

stat_self = function(df, method=RMSE, by_row=T, into_wide=F, is_corr=F, cores=F,
                     method_corr="pearson", use_corr="pairwise.complete.obs") {
  
  if (!is_corr) {
    df_stat = stat_pair(df, df, by_row=by_row, into_wide=into_wide, method=method, cores=cores)
  } else {
    df_stat = corr_pair(df, df, by_row=by_row, into_wide=into_wide, method=method_corr, use=use_corr)
  }
  return(df_stat)
}

stat_pair_apply = function(df1, df2, method=RMSE, by_row=T) {
  
  names_ = ifelse_def(by_row, rownames, colnames)
  names_int = intersect(names_(df1), names_(df2))
  df1 = ifelse_def(by_row, df1[names_int, ], df1[, names_int])
  df2 = ifelse_def(by_row, df2[names_int, ], df2[, names_int])
  
  stat_row = function(x) method(as.numeric(df1[x, ]), as.numeric(df2[x, ]))
  stat_col = function(x) method(as.numeric(df1[, x]), as.numeric(df2[, x]))
  stat = sapply(names_int, ifelse_def(by_row, stat_row, stat_col))
  
  names(stat) = names_int
  return(stat)
}

write_xlsx_def = function(..., is_list=F, file=NULL, sheets=NULL, rowNames=F, colNames=T) {
  Args = list(...)
  if (is_list) Args = Args %>% unlist(recursive=F)
  args_class = Args %>% sapply(class)
  idx_0 = ifelse(args_class=="data.frame", Args %>% sapply(nrow), Args %>% sapply(length))
  Args[idx_0==0] = Args[idx_0==0] %>% lapply(function(df) x[1]="(No Results)")
  
  if (is.null(sheets)) {
    openxlsx::write.xlsx(Args, file=file, rowNames=rowNames, colNames=colNames)
  } else {
    openxlsx::write.xlsx(Args, file=file, rowNames=rowNames, colNames=colNames, sheetName=sheets)
  }
}

multicores = function(cores=T, type="FORK", outfile="") {
  suppressMessages(library(parallel))
  suppressMessages(library(doParallel))
  if (isTRUE(cores)) cores = detectCores()/2
  cluster = makeCluster(cores, type=type, outfile=outfile)
  registerDoParallel(cluster)
  return(cluster)
}


# clust_future = function(cores=T, type="multisession", outfile="") {
#   
#   suppressMessages(library(future))
#   suppressMessages(library(doFuture))
#   suppressMessages(library(parallel))
#   # sequential, multisession, multicore, cluster
#   
#   if (isTRUE(cores)) cores = availableCores()/2
#   cl = makeCluster(cores, outfile=outfile)
#   plan(type, workers=cl)
#   return(cl)
# }



##### Difference between missing & is.null
# 
# speak1 = function(x) {
#   if (!missing(x)) {
#     print(sprintf("x : %s", x))
#   } else print("x : Empty")
# }
# 
# speak2 = function(x) {
#   if (!is.null(x)) {
#     print(sprintf("x : %s", x))
#   } else print("x : Empty")
# }
# 
# speak3 = function(x=NULL) {
#   if (!is.null(x)) {
#     print(sprintf("x : %s", x))
#   } else print("x : Empty")
# }
# 
# speak1()      # x : Empty
# speak1(NULL)  # character(0)
# speak1(1)     # x : 1
# 
# speak2()      # [Error Messages]
# speak2(NULL)  # x : Empty
# speak2(1)     # x : 1
# 
# speak3()      # x : Empty
# speak3(NULL)  # x : Empty
# speak3(1)     # x : 1

# plot_ex = function(df, x, y, color="black", size=1) {
#   
#   mapping = aes(x={{x}}, y={{y}}, color={{color}})
#   in_col = deparse(substitute(color)) %in% colnames(df)
#   if (!in_col) mapping$colour = NULL
#   
#   size = deparse(substitute(size))
#   draw_draft = sprintf("geom_point(size=%s)", size)
#   pl = df %>% ggplot(mapping) + object(draw_draft, envir=df)
#   print(pl)
# }
# 
# score = rep(1:3, each=4)
# class = rep(c("A", "B"), each=6)
# df = data.frame(PC1=rnorm(12), PC2=rnorm(12), Class=class, Score=score)
# plot_ex(df, PC1, PC2, color=Class, size=1)
# plot_ex(df, PC1, PC2, color=Class, size=Score)
# 
# print_self = function(x="Hello") {
#   print(deparse(substitute(x)))
#   print(deparse(quote(x)))
# }
# 
# print_self()
# # substitute  "\"Hello\""
# # quote       "x"
# 
# print_self(V1)
# # substitute  "V1"
# # quote       "x"
# 
# 
# args = list(size=5, color="red")
# df %>% ggplot(aes(V1, V2)) + do.call(geom_point, args)
# 
# draw_draft = "geom_point(size=score)"
# df %>% ggplot(aes(PC1, PC2)) + object(draw_draft)

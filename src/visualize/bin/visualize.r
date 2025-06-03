#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

theme_ih2025 <- function(xlim = NULL, ylim = NULL) {
  f <- theme_minimal() %+replace%
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.background = element_rect(color = "black", fill = "transparent"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.byrow = TRUE
  )
  return(f)
}

scale_color_nk2023 <- function() {
  f <- scale_color_manual(
    values = c("Hth" = "#C76F6B",
               "Between Hth & Hth/Opa" = "#EFA900",
               "Hth/Opa" = "#EFC8B9",
               "Opa/Erm" = "#FEE699",
               "Erm/Ey" = "#79A68C",
               "Ey/Hbn" = "#82C9C5",
               "Hbn/Opa/Slp" = "#B4C6E6",
               "Slp/D" = "#A6A1CD",
               "D/BH-1" = "#A46690")
  )
  return(f)
}

scale_color_subsystem <- function() {
  f <- scale_color_manual(
    values = c(
      Color = "#7FC97F",
      Object = "#BEAED4",
      Motion = "#FDC086",
      Luminance = "#FFFF99",
      Unannotated = "#386CB0",
      Polarization = "#F0027F",
      Form = "#BF5B17"
    )
  )
  return(f)
}

scale_color_type <- function() {
  types <- unique(opc_anno$cell_type)
  hues <- seq(20, 380, length.out = length(types) + 1)[seq_along(types)]
  palette <- hsv(h = (hues %% 360) / 360, s = 0.5, v = 0.5)
  set.seed(1)
  palette <- sample(palette)
  names(palette) <- types
  return(
    scale_color_manual(values = palette)
  )
}

filter_temporal <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[temporal_label != "unknown" & Confident_annotation == "Y"]
  ann$temporal_label <- factor(
    ann$temporal_label,
    levels = unique(ann$temporal_label),
    labels = unique(ann$temporal_label)
  )
  coord <- merge(
    coord,
    ann[, .(cell_type, temporal_label, Notch)],
    by.x = by_x,
    by.y = "cell_type"
  )
  return(coord)
}

filter_subsystem_ann <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[func != "unknown" & Confident_annotation == "Y"]
  ann <- ann[func != "Unannotated"]
  ann$func <- factor(ann$func)
  coord <- merge(
    coord,
    ann[, .(cell_type, func, Notch)],
    by.x = by_x,
    by.y = "cell_type"
  )
  return(coord)
}

filter_new_type <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[new == "Y"]
  coord <- merge(
    coord,
    ann[, .(cell_type, func, Notch)],
    by.x = by_x,
    by.y = "cell_type"
  )
  coord$cell_type <- coord[[by_x]]
  return(coord)
}

rsubsample <- function(x, n = 1e4, seed = 1) {
  if (!is.null(nrow) && nrow(x) > n) {
    set.seed(seed)
    if (is.data.table(x)) {
      out <- x[sample(.N, n)]
    } else {
      out <- x[sample(seq_len(nrow(x)), n), , with = FALSE]
    }
    return(out)
  }
  return(x)
}

filsplit <- function(x, f, slimit = 100L) {
  slist <- split(x, f)
  to_keep <- vapply(slist, function(x) nrow(x) > slimit, FUN.VALUE = logical(1))
  return(slist[to_keep])
}

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

### TODO
if (interactive()) {
  argvs$np <- "ME_L"
  argvs$syn_type <- "post"
  argvs$use_preset <- "temporal"
  argvs$density <- "per_type"
  argvs$ann <- "data/visual_neurons_20250602.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$subsample <- 10000L
  argvs$sparse_limit <- 100L
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$subsample <- as.integer(argvs$subsample)
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  syn_path <- argvs$synf
}

## Generate


## Load plot metadata and presets
plot_meta <- fread(argvs$meta)[neuropil == argvs$np]
preset <- fread(argvs$preset)[preset == argvs$use_preset]

## Load annotations and coordinates
opc_anno <- fread(argvs$ann)
np_coord <- fread(syn_path)

## Dynamic layout generation
filter_func <- get(preset$filter_func)
color_func <- get(preset$palette)
scale_axis_1 <- get(plot_meta$axis_1_func)
scale_axis_2 <- get(plot_meta$axis_2_func)
x_axis <- paste(argvs$syn_type, plot_meta$x_axis, sep = "_")
y_axis <- paste(argvs$syn_type, plot_meta$y_axis, sep = "_")

## Filter by annotation and subsample
np_coord <- filter_func(np_coord, opc_anno, syn_type = argvs$syn_type)

if (preset$notch_split) {
  np_coord <- filsplit(np_coord, np_coord$Notch, slimit = argvs$sparse_limit)
} else {
  np_coord <- list(all = np_coord)
}

for (i in names(np_coord)) {
  np_plot <- rsubsample(np_coord[[i]], n = argvs$subsample)
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$use_preset, argvs$density, i,
    sep = "_"
  )
  
  ## Generate the dot plot
  dotp <- np_plot |>
    ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
    geom_point(aes(color = .data[[preset$color_by]])) +
    labs(color = preset$color_guide) +
    theme_ih2025() +
    scale_axis_1(limits = c(plot_meta$xmin, plot_meta$xmax)) +
    scale_axis_2(limits = c(plot_meta$ymin, plot_meta$ymax)) +
    color_func()
  
  ## Extract legends to prevent layout fluctuation
  legendsp <- get_plot_component(dotp, "guide-box", return_all = TRUE) 
  dotp <- dotp + theme(legend.position="none")
  
  if (is.list(legendsp)) {
    ## Drop empty elements
    to_keep <- sapply(legendsp, function(x) "gtable" %in% class(x))
    legendsp <- legendsp[to_keep]
    if (length(legendsp) == 1) {
      legendsp <- legendsp[[1]]
    } else {
      stop("Legend extraction error: There is more than 1 item.")
    }
  }
  
  ## Manual calculation and interpolation of density
  np_den <- np_coord[[i]][ , .(
    group = get(preset$color_by),
    type = get(paste(argvs$syn_type, "type", sep = "_")),
    depth = get(paste(argvs$syn_type, "rz", sep = "_"))
  )]
  
  if (plot_meta$zinvert) {
    np_den$depth <- np_den$depth * -1
  }
  
  np_den <- filsplit(np_den, np_den$group, slimit = argvs$sparse_limit)
  den_grid <- seq(plot_meta$minz, plot_meta$maxz, length.out = 1024)
  if (argvs$density == "asis") {
    np_den <- lapply(np_den, function(x) return(density(x$depth)))
    interpolated <- lapply(names(np_den), function(type) {
      d <- np_den[[type]]
      y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
      return(data.frame(x = den_grid, y = y_interp, type = type))
    })
    inter <- do.call(rbind, interpolated)
  } else {
    type_avg <- lapply(names(np_den), function(tn) {
      x <- np_den[[tn]]
      id_type <- filsplit(x, x$type, slimit = argvs$sparse_limit)
      den_type <- lapply(id_type, function(idt) return(density(idt$depth)))
      interpolated <- lapply(names(den_type), function(type) {
        d <- den_type[[type]]
        y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
        return(data.frame(x = den_grid, y_interp = y_interp, type = type))
      })
      inter <- do.call(rbind, interpolated)
      out <- data.frame(
        x = den_grid,
        y = tapply(inter$y_interp, inter$x, mean),
        type = tn
      )
      return(out)
    })
    inter <- do.call(rbind, type_avg)
  }
  
  denp <- inter |>
    ggplot(aes(x = x, y = y, color = type)) +
    geom_line() +
    theme_ih2025() +
    color_func() +
    guides(color = "none")
  
  if (plot_meta$zid == "y") denp <- denp + coord_flip()
  
  if (plot_meta$outlayout == "landscape") {
    outp <- (dotp | denp) +
      plot_layout(widths = c(plot_meta$outr1, plot_meta$outr2))
  } else {
    outp <- (dotp / denp) +
      plot_layout(heights = c(plot_meta$outr1, plot_meta$outr2))
  }
  ggsave(
    plot = outp, filename = paste0(out_prefix, ".pdf"),
    width = plot_meta$outd1, height = plot_meta$outd2
  )
  ggsave(
    plot = legendsp, paste0(out_prefix, "_legend.pdf")
  )
}



# TODO: Get test results
# kwt <- readRDS(args$test)
# syn_type <- args$syn
# color_by <- args$color_by
# min <- as.integer(args$min)






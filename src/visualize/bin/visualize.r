#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))












argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

### TODO
if (interactive()) {
  argvs$np <- "LOP_R"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "temporal_known"
  argvs$density <- "asis"
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

if (preset$do_highlight && preset$hl_type == "label") {
  np_coord[, highlight := get(preset$hl_col) == preset$hl_val]
}
if (preset$do_highlight && preset$hl_type == "path") {
  if (!file.exists(preset$hl_val)) {
    stop("Cannot find highlight label file.")
  } else {
    hl_labels <- readLines(preset$hl_val)
  }
  np_coord[, highlight := get(preset$hl_col) %in% hl_labels]
}

if (preset$notch_split) {
  np_coord <- filsplit(np_coord, np_coord$Notch, slimit = argvs$sparse_limit)
} else {
  np_coord <- list(all = np_coord)
}

for (i in names(np_coord)) {
  if (preset$do_highlight) {
    if (sum(np_coord[[i]]$highlight) <= argvs$sparse_limit) {
      next
    }
  } else {
    if (nrow(np_coord[[i]]) <= argvs$sparse_limit) {
      next
    }
  }

  
  np_plot <- rsubsample(np_coord[[i]], n = argvs$subsample)
  if (preset$do_highlight && !any(np_plot$highlight)) {
    # For rare synapses that are not subsampled, we add 5 synapses as
    # representative synapses for visualization.
    pos_dt <- np_coord[[i]][highlight == TRUE]
    if (nrow(pos_dt) > 5) {
      rep_syn <- rsubsample(pos_dt, n = 5)
    } else {
      rep_syn <- pos_dt
    }
    np_plot <- rbind(np_plot, rep_syn)
  }
  
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$use_preset, argvs$density, i,
    sep = "_"
  )
  
  ## Generate the dot plot
  if (preset$do_highlight) {
    dotp <- np_plot |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(
        aes(color = .data[[preset$color_by]], alpha = highlight), dpi = 450
      )) +
      guides(alpha = "none") +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.05))
  } else {
    dotp <- np_plot |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(aes(color = .data[[preset$color_by]])), dpi = 450)
  }
  dotp <- dotp +
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
  if (preset$do_highlight) {
    np_den <- np_coord[[i]][highlight == TRUE, .(
      group = get(preset$color_by),
      type = get(paste(argvs$syn_type, "type", sep = "_")),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    np_den <- np_coord[[i]][ , .(
      group = get(preset$color_by),
      type = get(paste(argvs$syn_type, "type", sep = "_")),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  
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






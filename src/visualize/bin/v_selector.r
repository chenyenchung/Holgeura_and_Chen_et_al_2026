#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

options("ggrastr.default.dpi" = 450)

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (interactive()) {
  argvs$np <- "LOP_R"
  argvs$syn_type <- "post"
  argvs$ann <- "data/visual_neurons_20250602.csv"
  argvs$ts <- "data/media-2.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$tslut <- "data/to_selector.csv"
  argvs$density <- "pertype"
  argvs$subsample <- 10000L
  argvs$sparse_limit <- 100L
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$subsample <- as.integer(argvs$subsample)
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  syn_path <- argvs$synf
}


## Load plot metadata and presets
plot_meta <- fread(argvs$meta)[neuropil == argvs$np]

# Generated 
x_axis <- paste(argvs$syn_type, plot_meta$x_axis, sep = "_")
y_axis <- paste(argvs$syn_type, plot_meta$y_axis, sep = "_")
scale_axis_1 <- get(plot_meta$axis_1_func)
scale_axis_2 <- get(plot_meta$axis_2_func)
color_func <- function() {
  f <- scale_color_manual(
    values = c(
      "Notch On" = "#008800", "Notch Off" = "#990099", "nil" = "grey75"
    ),
    breaks = c("Notch On", "Notch Off")
  )
  return(f)
}

## Load annotations, coordinates, and terminal selector expression
## (Ozel et al., 2022)
opc_anno <- fread(argvs$ann)
np_coord <- fread(syn_path)
ts <- fread(argvs$ts)
tslut <- fread(argvs$tslut)

## Validate that all data points fall within specified axis limits
x_values <- np_coord[[x_axis]]
y_values <- np_coord[[y_axis]]

# Check X axis bounds
x_out_of_bounds <- sum(x_values < plot_meta$xmin | x_values > plot_meta$xmax, na.rm = TRUE)
if (x_out_of_bounds > 0) {
  stop(sprintf("ERROR: %d data points fall outside X-axis limits [%.2f, %.2f]. Check plot_meta boundaries for %s.", 
               x_out_of_bounds, plot_meta$xmin, plot_meta$xmax, argvs$np))
}

# Check Y axis bounds  
y_out_of_bounds <- sum(y_values < plot_meta$ymin | y_values > plot_meta$ymax, na.rm = TRUE)
if (y_out_of_bounds > 0) {
  stop(sprintf("ERROR: %d data points fall outside Y-axis limits [%.2f, %.2f]. Check plot_meta boundaries for %s.", 
               y_out_of_bounds, plot_meta$ymin, plot_meta$ymax, argvs$np))
}

# Annotate and filter to keep only annotated synapses
npp <- filter_default(np_coord, opc_anno, ts, tslut, syn_type = argvs$syn_type)
np_plot <- rsubsample(npp, n = argvs$subsample)

ts_symbols <- colnames(np_plot)[
  # Assuming the only logical columns are all selector statuses
  vapply(np_plot, is.logical, logical(1))
]

for (i in ts_symbols) {
  if (!any(npp[[i]])) {
    next
  }
  
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$density, i,
    sep = "_"
  )
  
  np_plot$color_label <- np_plot$Notch
  np_plot$color_label[!np_plot[[i]]] <- "nil"
  
  if (all(!np_plot[[i]])) {
    # For rare synapses that are not subsampled, we add 5 synapses as
    # representative synapses for visualization.
    pos_dt <- npp[get(i) == TRUE]
    if (nrow(pos_dt) > 5) {
      rep_syn <- rsubsample(pos_dt, n = 5)
    } else {
      rep_syn <- pos_dt
    }
    rep_syn$color_label <- rep_syn$Notch
    rep_syn$color_label[!rep_syn[[i]]] <- "nil"
    np_plot <- rbind(np_plot, rep_syn)
  }
  
  ## Generate the dot plot
  dotp <- np_plot |>
    ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
    rasterize(geom_point(aes(color = color_label))) +
    labs(color = "Notch Status") +
    theme_ih2025() +
    color_func() +
    scale_axis_1(limits = c(plot_meta$xmin, plot_meta$xmax)) +
    scale_axis_2(limits = c(plot_meta$ymin, plot_meta$ymax))
  
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
  np_den <- npp[ , .(
    group = Notch,
    type = get(paste(argvs$syn_type, "type", sep = "_")),
    depth = get(paste(argvs$syn_type, "rz", sep = "_"))
  )]
  np_den$group[!npp[[i]]] <- "nil"
  
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
      if (length(den_type) > 0) {
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
      } else {
        out <- data.frame(
          x = den_grid,
          y = 0,
          type = tn
        )
      }
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
  outp <- outp +
    plot_annotation(
      title = i,
      subtitle = paste0(argvs$syn_type, "synapses")
    ) & 
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12, face = "italic")
    )
  
  ggsave(
    plot = outp, filename = paste0(out_prefix, ".pdf"),
    width = plot_meta$outd1, height = plot_meta$outd2
  )
  ggsave(
    plot = legendsp, paste0(out_prefix, "_legend.pdf")
  )
}

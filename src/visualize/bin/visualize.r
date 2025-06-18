#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

options("ggrastr.default.dpi" = 450)

# Validation functions
validate_config_file <- function(file_path, required_cols) {
  if (!file.exists(file_path)) stop(sprintf("Config file not found: %s", file_path))
  data <- fread(file_path)
  if (nrow(data) == 0) stop(sprintf("Config file is empty: %s", file_path))
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) stop(sprintf("Missing columns in %s: %s", file_path, paste(missing_cols, collapse=", ")))
  return(data)
}

validate_functions <- function(preset, plot_meta) {
  # Check filter function
  if (!exists(preset$filter_func, mode="function")) 
    stop(sprintf("Filter function not found: %s", preset$filter_func))
  
  # Check color function  
  if (!exists(preset$palette, mode="function"))
    stop(sprintf("Color function not found: %s", preset$palette))
    
  # Check axis functions
  if (!exists(plot_meta$axis_1_func, mode="function"))
    stop(sprintf("Axis function not found: %s", plot_meta$axis_1_func))
  if (!exists(plot_meta$axis_2_func, mode="function"))
    stop(sprintf("Axis function not found: %s", plot_meta$axis_2_func))
}

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

### TODO
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$np <- "LO_L"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "temporal_new"
  argvs$density <- "pertype"
  argvs$ann <- "data/visual_neurons_anno.csv"
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


## Load plot metadata and presets with validation
plot_meta_all <- validate_config_file(argvs$meta, c("neuropil", "axis_1_func", "axis_2_func", "x_axis", "y_axis", "xmin", "xmax", "ymin", "ymax", "minz", "maxz", "zinvert", "zid", "outlayout", "outr1", "outr2", "outd1", "outd2"))
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) stop(sprintf("No metadata found for neuropil: %s", argvs$np))

preset_all <- validate_config_file(argvs$preset, c("preset", "filter_func", "palette", "color_by", "color_guide", "do_highlight", "hl_type", "hl_col", "hl_val", "notch_split"))
preset <- preset_all[preset == argvs$use_preset]
if (nrow(preset) == 0) stop(sprintf("No preset found: %s", argvs$use_preset))

# Validate that required functions exist
validate_functions(preset, plot_meta)

## Load annotations and coordinates
opc_anno <- validate_config_file(argvs$ann, c("cell_type"))
if (!file.exists(syn_path)) stop(sprintf("Synapse coordinate file not found: %s", syn_path))
np_coord <- fread(syn_path)

## Validate required columns in coordinate data
required_coord_cols <- c(
  paste0(argvs$syn_type, "_type"),
  paste0(argvs$syn_type, "_rz"),
  paste(argvs$syn_type, plot_meta$x_axis, sep = "_"),
  paste(argvs$syn_type, plot_meta$y_axis, sep = "_")
)
missing_coord_cols <- setdiff(required_coord_cols, colnames(np_coord))
if (length(missing_coord_cols) > 0) {
  stop(sprintf("Missing required columns in coordinate data: %s", paste(missing_coord_cols, collapse=", ")))
}

## Dynamic layout generation
filter_func <- get(preset$filter_func)
color_func <- get(preset$palette)
scale_axis_1 <- get(plot_meta$axis_1_func)
scale_axis_2 <- get(plot_meta$axis_2_func)
x_axis <- paste(argvs$syn_type, plot_meta$x_axis, sep = "_")
y_axis <- paste(argvs$syn_type, plot_meta$y_axis, sep = "_")

## Filter by annotation and subsample
np_coord <- filter_func(np_coord, opc_anno, syn_type = argvs$syn_type)

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

# Filter out types with too few synapses for consistent visualization
type_col <- paste0(argvs$syn_type, "_type")
type_counts <- np_coord[, .N, by = c(type_col)]
keep_types <- type_counts[N >= argvs$sparse_limit, get(type_col)]
dropped_types <- type_counts[N < argvs$sparse_limit]
if (nrow(dropped_types) > 0) {
  message(sprintf("Dropped %d %s types with < %d synapses", 
                  nrow(dropped_types), argvs$syn_type, argvs$sparse_limit))
}
np_coord <- np_coord[get(type_col) %in% keep_types]

# Filter out types with insufficient variation for density calculation
depth_col <- paste0(argvs$syn_type, "_rz")
if (depth_col %in% colnames(np_coord)) {
  depth_variation <- np_coord[, .(unique_depths = length(unique(get(depth_col)))), by = c(type_col)]
  valid_types <- depth_variation[unique_depths >= 2, get(type_col)]
  if (length(valid_types) < nrow(depth_variation)) {
    message(sprintf("Dropped %d types with insufficient depth variation for density calculation", 
                    nrow(depth_variation) - length(valid_types)))
  }
  np_coord <- np_coord[get(type_col) %in% valid_types]
}

if (preset$do_highlight && preset$hl_type == "label") {
  if (!preset$hl_col %in% colnames(np_coord)) {
    stop(sprintf("Highlight column not found in data: %s", preset$hl_col))
  }
  np_coord[, highlight := get(preset$hl_col) == preset$hl_val]
}
if (preset$do_highlight && preset$hl_type == "path") {
  if (!file.exists(preset$hl_val)) {
    stop("Cannot find highlight label file.")
  } else {
    hl_labels <- readLines(preset$hl_val)
  }
  if (!preset$hl_col %in% colnames(np_coord)) {
    stop(sprintf("Highlight column not found in data: %s", preset$hl_col))
  }
  np_coord[, highlight := get(preset$hl_col) %in% hl_labels]
}


np_raw <- np_coord
if (grepl("^ME", argvs$np)) {
  np_coord <- rsubsample(np_coord, n = argvs$subsample * 4)
} else {
  np_coord <- np_coord[!grepl("^(Mi|Dm|Pm)", get(paste0(argvs$syn_type, "_type")))]
  np_coord <- rsubsample(np_coord, n = argvs$subsample * 2)
}


if (preset$notch_split) {
  required_notch_cols <- c("Notch", "ntype")
  missing_notch_cols <- setdiff(required_notch_cols, colnames(np_coord))
  if (length(missing_notch_cols) > 0) {
    stop(sprintf("Missing columns for notch_split: %s", paste(missing_notch_cols, collapse=", ")))
  }
  np_coord$notch_ntype <- paste(np_coord$Notch, np_coord$ntype, sep = "_")
  np_coord <- split(np_coord, np_coord$notch_ntype)
  np_raw$notch_ntype <- paste(np_raw$Notch, np_raw$ntype, sep = "_")
  np_raw <- split(np_raw, np_raw$notch_ntype)
} else {
  np_coord <- list(all = np_coord)
  np_raw <- list(all = np_raw)
}

for (i in names(np_coord)) {
  # Skip if no data after filtering
  if (nrow(np_coord[[i]]) == 0) {
    next
  }
  

  if (preset$do_highlight && !any(np_coord[[i]]$highlight)) {
    # For rare synapses that are not subsampled, we add 5 synapses as
    # representative synapses for visualization.
    pos_dt <- np_raw[[i]][highlight == TRUE]
    if (nrow(pos_dt) > 5) {
      rep_syn <- rsubsample(pos_dt, n = 5)
    } else {
      rep_syn <- pos_dt
    }
    np_coord[[i]] <- rbind(np_coord[[i]], rep_syn)
  }
  
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$use_preset, argvs$density, i,
    sep = "_"
  )
  
  ## Generate the dot plot
  if (preset$do_highlight) {
    dotp <- np_coord[[i]] |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(
        aes(color = .data[[preset$color_by]], alpha = highlight)
      )) +
      guides(alpha = "none") +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.05))
  } else {
    dotp <- np_coord[[i]] |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(aes(color = .data[[preset$color_by]])))
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
  
  np_den <- split(np_den, np_den$group)
  
  
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
    # No need to filter here since we already filtered at the type level
    if (length(np_den) < 1) {
      next
    }
    type_avg <- lapply(names(np_den), function(tn) {
      x <- np_den[[tn]]
      id_type <- split(x, x$type)
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
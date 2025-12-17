#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

# Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$connections <- 'data/connections_princeton_no_threshold.csv.gz'
  argvs$celltypes <- 'data/consolidated_cell_types.csv.gz'
  argvs$rotated_dir <- 'int/idv_mat'
  argvs$viz_meta <- 'data/viz_meta.csv'
  argvs$output_dir <- 'plots/lc_lplc_per_neuron'
  argvs$syn_threshold <- 5
} else {
  # Ensure utils.r is sourced
  if (!file.exists("./utils.r") && file.exists("src/utils.r")) {
    source("src/utils.r", chdir = FALSE)
  }
}

# Set default synapse threshold if not provided
if (is.null(argvs$syn_threshold)) {
  argvs$syn_threshold <- 5
} else {
  argvs$syn_threshold <- as.numeric(argvs$syn_threshold)
}

# Validate input files
if (is.null(argvs$connections) || !file.exists(argvs$connections)) {
  stop("Connection file not found: ", argvs$connections)
}

if (is.null(argvs$celltypes) || !file.exists(argvs$celltypes)) {
  stop("Cell type annotation file not found: ", argvs$celltypes)
}

if (is.null(argvs$viz_meta) || !file.exists(argvs$viz_meta)) {
  stop("Visualization metadata file not found: ", argvs$viz_meta)
}

if (is.null(argvs$rotated_dir) || !dir.exists(argvs$rotated_dir)) {
  stop("Rotated matrix directory not found: ", argvs$rotated_dir)
}

if (is.null(argvs$output_dir)) {
  argvs$output_dir <- 'plots/lc_lplc_per_neuron'
}

# Create output directory if it doesn't exist
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", argvs$output_dir))
}

cat("========================================\n")
cat("LC/LPLC Per-Neuron Partner Visualization\n")
cat("========================================\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Identify partner neurons for a given LC/LPLC type
#'
#' @param connections Connection data.table with columns: pre_root_id, post_root_id,
#'                    neuropil, syn_count, pre_type, post_type
#' @param target_lc_lplc Target LC/LPLC type (e.g., "LC10a")
#' @param syn_threshold Minimum synapse count to define partnership (default: 5)
#' @return Vector of partner neuron root_ids
identify_partner_neurons <- function(connections, target_lc_lplc, syn_threshold = 5) {
  # Filter for optic lobe neuropils
  ol_npils <- c("ME_L", "ME_R", "LO_L", "LO_R", "LOP_L", "LOP_R")

  # Get connections where post_type == target LC/LPLC and syn_count > threshold
  partners <- connections[
    post_type == target_lc_lplc &
    syn_count > syn_threshold &
    neuropil %in% ol_npils,
    unique(pre_root_id)
  ]

  return(partners)
}

#' Annotate synapses with partner status and perform stratified subsampling
#'
#' @param coord Synapse coordinate data.table
#' @param partner_neurons Vector of partner neuron root_ids
#' @param sample_size Total number of synapses to subsample
#' @return Subsampled data.table with is_partner column
annotate_and_subsample <- function(coord, partner_neurons, sample_size) {
  # Add binary is_partner flag
  coord[, is_partner := pre_root_id %in% partner_neurons]

  # Calculate original ratio
  n_partner <- sum(coord$is_partner)
  n_total <- nrow(coord)

  if (n_total == 0) {
    return(coord)
  }

  partner_ratio <- n_partner / n_total

  # If total is already less than sample size, return all
  if (n_total <= sample_size) {
    return(coord)
  }

  # Calculate stratified sample sizes
  n_partner_sample <- round(sample_size * partner_ratio)
  n_nonpartner_sample <- sample_size - n_partner_sample

  # Sample separately
  set.seed(1)
  partner_data <- coord[is_partner == TRUE]
  nonpartner_data <- coord[is_partner == FALSE]

  # Handle edge cases
  if (nrow(partner_data) == 0) {
    # No partners - sample all from non-partners
    result <- nonpartner_data[sample(.N, min(.N, sample_size))]
  } else if (nrow(nonpartner_data) == 0) {
    # All partners - sample all from partners
    result <- partner_data[sample(.N, min(.N, sample_size))]
  } else {
    # Normal case - stratified sampling
    partner_sample <- partner_data[sample(.N, min(.N, n_partner_sample))]
    nonpartner_sample <- nonpartner_data[sample(.N, min(.N, n_nonpartner_sample))]
    result <- rbind(partner_sample, nonpartner_sample)
  }

  return(result)
}

#' Create per-neuron partner visualization plot
#'
#' @param coord Annotated coordinate data.table with is_partner column
#' @param meta Visualization metadata for the neuropil
#' @param lc_lplc_type Target LC/LPLC type
#' @param neuropil Neuropil name
#' @return ggplot object
create_per_neuron_plot <- function(coord, meta, lc_lplc_type, neuropil) {
  # Extract axes from metadata
  x_axis <- paste0("pre_", meta$x_axis)
  y_axis <- paste0("pre_", meta$y_axis)

  # Validate axes exist
  if (!x_axis %in% colnames(coord) || !y_axis %in% colnames(coord)) {
    warning(sprintf("Missing axes for %s: %s or %s", neuropil, x_axis, y_axis))
    return(NULL)
  }

  # Count statistics
  n_partner <- sum(coord$is_partner)
  n_total <- nrow(coord)
  n_nonpartner <- n_total - n_partner

  # Sort so partner synapses (blue) are plotted on top
  setorder(coord, is_partner)  # FALSE first, TRUE last

  # Create scatter plot
  p <- ggplot(coord, aes_string(x = x_axis, y = y_axis, color = "is_partner")) +
    rasterise(
      geom_point(size = 0.5, alpha = 0.6),
      dpi = 450
    ) +
    scale_color_manual(
      values = c("TRUE" = "blue", "FALSE" = "grey92"),
      name = "Presynaptic Partner",
      labels = c("FALSE" = "Non-partner", "TRUE" = "Partner")
    ) +
    labs(
      title = sprintf("Per-Neuron Presynaptic Partners of %s in %s",
                      lc_lplc_type, neuropil),
      subtitle = sprintf("n = %d synapses | Partners: %d (%.1f%%) | Non-partners: %d (%.1f%%)",
                        n_total, n_partner, 100*n_partner/n_total,
                        n_nonpartner, 100*n_nonpartner/n_total)
    ) +
    coord_fixed(ratio = 1) +
    theme_ih2025()

  # Apply axis scaling functions from metadata (handles reversals)
  # Note: axis_1/2 refer to metadata ordering, not plot x/y axes
  if (grepl("reverse", meta$axis_1_func)) {
    min1 <- meta$max1
    max1 <- meta$min1
    meta$min1 <- min1
    meta$max1 <- max1
  }
  if (grepl("reverse", meta$axis_2_func)) {
    min2 <- meta$max2
    max2 <- meta$min2
    meta$min2 <- min2
    meta$max2 <- max2
  }
  scale_axis_1 <- get(meta$axis_1_func)
  scale_axis_2 <- get(meta$axis_2_func)

  p <- p +
    scale_axis_1(limits = c(meta$min1, meta$max1)) +
    scale_axis_2(limits = c(meta$min2, meta$max2))

  return(p)
}

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Load data
cat("Loading data...\n")
connections <- fread(argvs$connections)
cat(sprintf("  Loaded connections: %d rows\n", nrow(connections)))

type_ann <- fread(argvs$celltypes)
cat(sprintf("  Loaded cell type annotations: %d neurons\n", nrow(type_ann)))

viz_meta <- fread(argvs$viz_meta)
cat(sprintf("  Loaded visualization metadata: %d neuropils\n", nrow(viz_meta)))

# Normalize root IDs by stripping prefix
cat("\nNormalizing root IDs...\n")
connections[, pre_root_id := sub("^720575940", "", pre_root_id)]
connections[, post_root_id := sub("^720575940", "", post_root_id)]
type_ann[, root_id := sub("^720575940", "", root_id)]

# Annotate connections with cell types
cat("Annotating connections with cell types...\n")

# Merge presynaptic types
connections <- merge(connections, type_ann[, .(root_id, primary_type)],
                     by.x = "pre_root_id", by.y = "root_id", all.x = TRUE)
setnames(connections, "primary_type", "pre_type")

# Merge postsynaptic types
connections <- merge(connections, type_ann[, .(root_id, primary_type)],
                     by.x = "post_root_id", by.y = "root_id", all.x = TRUE)
setnames(connections, "primary_type", "post_type")

# Filter for LC/LPLC postsynaptic types
connections <- connections[grepl("^LC[0-9]+|^LPLC[0-9]+", post_type)]

cat(sprintf("  Filtered to %d LC/LPLC connections\n", nrow(connections)))

# Get unique LC/LPLC types
lc_lplc_types <- sort(unique(connections$post_type))
cat(sprintf("\nFound %d LC/LPLC types to visualize:\n", length(lc_lplc_types)))
cat(sprintf("  %s\n", paste(lc_lplc_types, collapse=", ")))

# Get neuropils from metadata
neuropils <- unique(viz_meta$neuropil)
cat(sprintf("\nProcessing %d neuropils: %s\n\n",
            length(neuropils), paste(neuropils, collapse=", ")))

# Main visualization loop
plots_created <- 0
plots_skipped <- 0

cat("Generating plots...\n")
cat(sprintf("Synapse threshold: >%d synapses\n", argvs$syn_threshold))
cat("Progress: [LC/LPLC type] × [neuropil]\n\n")

for (lc_lplc_type in lc_lplc_types) {
  cat(sprintf("Processing %s...\n", lc_lplc_type))

  # Identify partner neurons for this LC/LPLC type
  partner_neurons <- identify_partner_neurons(
    connections,
    lc_lplc_type,
    syn_threshold = argvs$syn_threshold
  )

  cat(sprintf("  Found %d partner neurons\n", length(partner_neurons)))

  for (np in neuropils) {
    # Load rotated coordinate matrix
    coord_file <- file.path(argvs$rotated_dir,
                            paste0(np, "_rotated.csv.gz"))

    if (!file.exists(coord_file)) {
      warning(sprintf("  %s: coordinate file not found (skipped)", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    coord <- fread(coord_file)

    # Determine sample size
    sample_size <- ifelse(grepl("^ME_", np), 40000, 20000)

    # Annotate and stratified subsample
    n_before <- nrow(coord)
    coord_sub <- annotate_and_subsample(coord, partner_neurons, sample_size)
    n_after <- nrow(coord_sub)

    # Skip if no data
    if (n_after == 0) {
      cat(sprintf("  %s: no data (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Get neuropil metadata
    meta <- viz_meta[viz_meta$neuropil == np, ]

    if (nrow(meta) == 0) {
      warning(sprintf("  %s: no metadata (skipped)", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Extract first row (should only be one)
    meta <- meta[1, ]

    # Create plot
    plot <- create_per_neuron_plot(
      coord_sub, meta, lc_lplc_type, np
    )

    if (is.null(plot)) {
      plots_skipped <- plots_skipped + 1
      next
    }

    # Save plot
    output_file <- file.path(
      argvs$output_dir,
      sprintf("%s_%s_per_neuron_partners.pdf", np, lc_lplc_type)
    )

    # Use dimensions from metadata (landscape: 10x6, portrait: 5x10)
    ggsave(output_file, plot, width = meta$outd1, height = meta$outd2)

    # Report progress
    n_partner <- sum(coord_sub$is_partner)
    if (n_before > n_after) {
      cat(sprintf("  %s: saved (%d→%d synapses, %.1f%% partners, subsampled)\n",
                  np, n_before, n_after, 100*n_partner/n_after))
    } else {
      cat(sprintf("  %s: saved (%d synapses, %.1f%% partners)\n",
                  np, n_after, 100*n_partner/n_after))
    }

    plots_created <- plots_created + 1
  }

  cat("\n")
}

# Summary
cat("========================================\n")
cat("Visualization complete!\n")
cat("========================================\n")
cat(sprintf("Plots created: %d\n", plots_created))
cat(sprintf("Plots skipped (no data): %d\n", plots_skipped))
cat(sprintf("Output directory: %s\n", argvs$output_dir))
cat(sprintf("\nTotal files expected: %d (%d LC/LPLC × %d neuropils)\n",
            length(lc_lplc_types) * length(neuropils),
            length(lc_lplc_types), length(neuropils)))

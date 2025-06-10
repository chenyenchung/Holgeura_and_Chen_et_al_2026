#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

#' Perform pairwise Kolmogorov-Smirnov tests between groups
#' @param data Data frame with depth and group columns
#' @param sparse_limit Minimum number of observations per group
#' @return Data frame with columns: group1, group2, pvalue, p_adj
perform_ks_tests <- function(data, sparse_limit = 100) {
  # Filter groups with sufficient data
  group_counts <- table(data$group)
  valid_groups <- names(group_counts)[group_counts > sparse_limit]
  
  if (length(valid_groups) < 2) {
    warning("Less than 2 groups have sufficient data for testing")
    return(data.frame(group1 = character(0), group2 = character(0), 
                     pvalue = numeric(0), p_adj = numeric(0)))
  }
  
  # Generate all pairwise combinations
  combinations <- combn(valid_groups, 2, simplify = FALSE)
  
  # Perform K-S tests
  results <- lapply(combinations, function(pair) {
    group1_data <- data[data$group == pair[1], "depth"]
    group2_data <- data[data$group == pair[2], "depth"]
    
    ks_result <- ks.test(group1_data, group2_data)
    
    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      pvalue = ks_result$p.value
    )
  })
  
  # Combine results
  ks_results <- do.call(rbind, results)
  
  # Apply multiple testing correction (Benjamini-Hochberg)
  ks_results$p_adj <- p.adjust(ks_results$pvalue, method = "BH")
  
  return(ks_results)
}

# Load helper functions from visualize.r
source("src/visualize/bin/visualize.r", chdir = FALSE)

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Interactive testing setup
if (interactive()) {
  argvs$np <- "LOP_R"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "temporal_known"
  argvs$density <- "asis"
  argvs$ann <- "data/visual_neurons_20250602.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$sparse_limit <- 100L
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  syn_path <- argvs$synf
}

# Load plot metadata and presets
plot_meta <- fread(argvs$meta)[neuropil == argvs$np]
preset <- fread(argvs$preset)[preset == argvs$use_preset]

# Load annotations and coordinates
opc_anno <- fread(argvs$ann)
np_coord <- fread(syn_path)

# Apply same filtering as visualization
filter_func <- get(preset$filter_func)
np_coord <- filter_func(np_coord, opc_anno, syn_type = argvs$syn_type)

# Handle highlighting if specified
if (preset$do_highlight && preset$hl_type == "label") {
  np_coord[, highlight := get(preset$hl_col) == preset$hl_val]
}
if (preset$do_highlight && preset$hl_type == "path") {
  if (file.exists(preset$hl_val)) {
    hl_labels <- readLines(preset$hl_val)
    np_coord[, highlight := get(preset$hl_col) %in% hl_labels]
  }
}

# Split by Notch if required
if (preset$notch_split) {
  np_coord_list <- filsplit(np_coord, np_coord$Notch, slimit = argvs$sparse_limit)
} else {
  np_coord_list <- list(all = np_coord)
}

# Perform K-S tests for each split
all_ks_results <- list()

for (i in names(np_coord_list)) {
  current_data <- np_coord_list[[i]]
  
  # Skip if insufficient data
  if (nrow(current_data) <= argvs$sparse_limit) {
    next
  }
  
  # Prepare data for K-S testing
  if (preset$do_highlight) {
    test_data <- current_data[highlight == TRUE, .(
      group = get(preset$color_by),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    test_data <- current_data[, .(
      group = get(preset$color_by),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  # Apply z-axis inversion if needed
  if (plot_meta$zinvert) {
    test_data$depth <- test_data$depth * -1
  }
  
  # Perform K-S tests
  ks_results <- perform_ks_tests(test_data, sparse_limit = argvs$sparse_limit)
  
  if (nrow(ks_results) > 0) {
    ks_results$split <- i
    ks_results$neuropil <- argvs$np
    ks_results$syn_type <- argvs$syn_type
    ks_results$preset <- argvs$use_preset
    all_ks_results[[i]] <- ks_results
  }
}

# Combine all results
if (length(all_ks_results) > 0) {
  final_results <- do.call(rbind, all_ks_results)
  
  # Generate output filename
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$use_preset, argvs$density,
    sep = "_"
  )
  output_file <- paste0(out_prefix, "_ks_tests.csv")
  
  # Write results
  fwrite(final_results, output_file)
  cat("K-S test results written to:", output_file, "\n")
  cat("Number of comparisons:", nrow(final_results), "\n")
  
  # Print summary
  sig_count <- sum(final_results$p_adj < 0.05, na.rm = TRUE)
  cat("Significant comparisons (p_adj < 0.05):", sig_count, "\n")
} else {
  cat("No valid comparisons could be performed\n")
}
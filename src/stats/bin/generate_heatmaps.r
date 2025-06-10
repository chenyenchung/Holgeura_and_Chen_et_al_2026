#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))

#' Generate combined heatmap with both p-values and effect sizes
#' @param test_results Data frame with test results
#' @param title_prefix Title prefix
#' @return List of ggplot objects
generate_combined_heatmaps <- function(test_results,
                                       title_prefix = "Quantile L2 Wasserstein Analysis") {
  
  conf_level <- if ("conf_level" %in% colnames(test_results)) {
    unique(test_results$conf_level)[1]
  } else {
    95
  }
  
  p1 <- generate_pvalue_heatmap(
    test_results,
    title = paste(title_prefix, "- Adjusted P-values"),
    show_values = TRUE
  )
  
  p2 <- generate_effect_size_heatmap(
    test_results,
    title = paste0(title_prefix, " - Effect Sizes (", conf_level, "% CI)")
  )
  
  return(list(pvalue_heatmap = p1, effect_heatmap = p2))
}


#' Generate enhanced p-value heatmap for Quantile L2 Wasserstein test results
#' @param test_results Data frame with columns: group1, group2, adj_pval, test_statistic_median
#' @param title Title for the heatmap
#' @param show_values Show p-values as text on heatmap
#' @return ggplot object
generate_pvalue_heatmap <- function(test_results, 
                                    title = "Quantile L2 Wasserstein Test P-values",
                                    show_values = TRUE) {
  
  # Get all unique groups
  all_groups <- unique(c(test_results$group1, test_results$group2))
  
  # Create a full matrix with all pairwise combinations
  pvalue_matrix <- matrix(NA, nrow = length(all_groups), ncol = length(all_groups))
  rownames(pvalue_matrix) <- all_groups
  colnames(pvalue_matrix) <- all_groups
  
  # Fill diagonal with 1 (groups compared to themselves)
  diag(pvalue_matrix) <- 1
  
  # Fill the matrix with adjusted p-values
  for (i in 1:nrow(test_results)) {
    g1 <- test_results$group1[i]
    g2 <- test_results$group2[i]
    p_adj <- test_results$adj_pval[i]
    
    # Fill both upper and lower triangles for symmetry
    pvalue_matrix[g1, g2] <- p_adj
    pvalue_matrix[g2, g1] <- p_adj
  }
  
  # Use original group order
  group_order <- all_groups
  
  # Convert to long format for ggplot and blank out upper triangle
  pvalue_df <- melt(pvalue_matrix, varnames = c("Group1", "Group2"), value.name = "adj_pval")
  pvalue_df <- pvalue_df[!is.na(pvalue_df$adj_pval), ]
  
  # Blank out upper triangle to remove redundancy
  pvalue_df$adj_pval[as.numeric(pvalue_df$Group1) > as.numeric(pvalue_df$Group2)] <- NA
  
  # Reorder factors based on clustering
  pvalue_df$Group1 <- factor(pvalue_df$Group1, levels = group_order)
  pvalue_df$Group2 <- factor(pvalue_df$Group2, levels = group_order)
  
  # Create significance categories
  pvalue_df$significance <- cut(pvalue_df$adj_pval, 
                                breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                labels = c("***", "**", "*", ".", "ns"),
                                include.lowest = TRUE)
  
  # Create heatmap
  p <- ggplot(pvalue_df, aes(x = Group1, y = Group2)) +
    geom_tile(aes(fill = adj_pval), color = "white", linewidth = 0.5) +
    geom_tile(data = pvalue_df[is.na(pvalue_df$adj_pval), ], fill = "white", color = "white") +
    scale_fill_gradientn(
      colors = c("red", "red", "white", "white"),
      values = c(0, 0.001, 0.05, 1),
      name = "Adj.\nP-value",
      trans = "log1p",
      breaks = c(0.001, 0.01, 0.05),
      labels = c("0.001", "0.01", "0.05"),
      limits = c(0, 0.1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    labs(title = title) +
    coord_fixed()
  
  # Add text annotations if requested
  if (show_values) {
    p <- p + geom_text(aes(label = significance), color = "black", size = 3)
  }
  
  return(p)
}

#' Generate effect size heatmap
#' @param test_results Data frame with test results
#' @param title Title for the heatmap
#' @param show_ci Show confidence intervals
#' @return ggplot object
generate_effect_size_heatmap <- function(test_results,
                                        title = "Effect Size (Median EMD)",
                                        show_ci = FALSE) {
  
  # Get all unique groups
  all_groups <- unique(c(test_results$group1, test_results$group2))
  
  # Create effect size matrix
  effect_matrix <- matrix(NA, nrow = length(all_groups), ncol = length(all_groups))
  rownames(effect_matrix) <- all_groups
  colnames(effect_matrix) <- all_groups
  
  # Fill diagonal with 0 (no difference within groups)
  diag(effect_matrix) <- 0
  
  # Fill the matrix with median effect sizes
  for (i in 1:nrow(test_results)) {
    g1 <- test_results$group1[i]
    g2 <- test_results$group2[i]
    effect <- test_results$test_statistic_median[i]
    
    # Fill both upper and lower triangles for symmetry
    effect_matrix[g1, g2] <- effect
    effect_matrix[g2, g1] <- effect
  }
  
  # Use original group order
  group_order <- all_groups
  
  # Convert to long format and blank out upper triangle
  effect_df <- melt(effect_matrix, varnames = c("Group1", "Group2"), value.name = "effect_size")
  effect_df <- effect_df[!is.na(effect_df$effect_size), ]
  
  # Blank out upper triangle to remove redundancy
  effect_df$effect_size[as.numeric(effect_df$Group1) > as.numeric(effect_df$Group2)] <- NA
  
  # Reorder factors
  effect_df$Group1 <- factor(effect_df$Group1, levels = group_order)
  effect_df$Group2 <- factor(effect_df$Group2, levels = group_order)
  
  # Create heatmap
  p <- ggplot(effect_df, aes(x = Group1, y = Group2, fill = effect_size)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_tile(data = effect_df[is.na(effect_df$effect_size), ], fill = "white", color = "white") +
    scale_fill_gradient(name = "Median\nEMD", low = "white", high = "red") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    labs(title = title) +
    coord_fixed()
  
  return(p)
}

#' Generate confidence interval plot
#' @param test_results Data frame with test results
#' @param title Plot title
#' @return ggplot object
generate_confidence_interval_plot <- function(test_results, 
                                             title = "Effect Sizes with Confidence Intervals") {
  
  # Create comparison labels
  test_results$comparison <- paste(test_results$group1, "vs", test_results$group2)
  
  # Order by effect size
  test_results <- test_results[order(test_results$test_statistic_median, decreasing = TRUE), ]
  test_results$comparison <- factor(test_results$comparison, levels = test_results$comparison)
  
  # Create significance indicator
  test_results$significant <- ifelse(test_results$adj_pval < 0.05, "Significant", "Not Significant")
  
  conf_level <- if ("conf_level" %in% colnames(test_results)) {
    unique(test_results$conf_level)[1]
  } else {
    95
  }
  
  p <- ggplot(test_results, aes(x = comparison, y = test_statistic_median, color = significant)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = test_statistic_lower, ymax = test_statistic_upper), 
                  width = 0.2, linewidth = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey60")) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0(title, " (", conf_level, "% CI)"),
      x = "Comparison",
      y = "Effect Size (Median EMD)",
      color = "Significance"
    )
  
  return(p)
}

#' Generate comprehensive statistical report plots
#' @param test_results Data frame with test results
#' @param output_prefix Output file prefix
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return List of generated file paths
generate_statistical_report <- function(test_results, output_prefix, 
                                       width = 10, height = 8) {
  
  generated_files <- c()
  
  # Generate combined heatmaps
  heatmaps <- generate_combined_heatmaps(test_results)
  
  # Save p-value heatmap
  pval_file <- paste0(output_prefix, "_pvalue_heatmap.pdf")
  ggsave(plot = heatmaps$pvalue_heatmap, filename = pval_file,
         width = width, height = height, dpi = 300)
  generated_files <- c(generated_files, pval_file)
  
  # Save effect size heatmap
  effect_file <- paste0(output_prefix, "_effect_size_heatmap.pdf")
  ggsave(plot = heatmaps$effect_heatmap, filename = effect_file,
         width = width, height = height, dpi = 300)
  generated_files <- c(generated_files, effect_file)
  
  # Save confidence interval plot
  ci_plot <- generate_confidence_interval_plot(test_results)
  ci_file <- paste0(output_prefix, "_confidence_intervals.pdf")
  ggsave(plot = ci_plot, filename = ci_file,
         width = width + 2, height = height, dpi = 300)
  generated_files <- c(generated_files, ci_file)
  
  return(generated_files)
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  argvs$input_file <- "LOP_L_pre_temporal_known_asis_quantile_wasserstein_enhanced.csv"
}

# Check if running as standalone script
if (!is.null(argvs$input_file)) {
  
  # Load test results
  if (!file.exists(argvs$input_file)) {
    stop("Input file not found: ", argvs$input_file)
  }
  
  test_results <- fread(argvs$input_file)
  
  # Set defaults
  output_prefix <- if (is.null(argvs$output_prefix)) {
    gsub("\\.csv$", "", basename(argvs$input_file))
  } else {
    argvs$output_prefix
  }
  
  width <- if (is.null(argvs$width)) 10 else as.numeric(argvs$width)
  height <- if (is.null(argvs$height)) 8 else as.numeric(argvs$height)
  
  # Generate plots
  cat("Generating statistical visualization plots...\n")
  generated_files <- generate_statistical_report(test_results, output_prefix, width, height)
  
  cat("Generated files:\n")
  for (file in generated_files) {
    cat("  ", file, "\n")
  }
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  cat("Total comparisons:", nrow(test_results), "\n")
  cat("Significant (adj_pval < 0.05):", sum(test_results$adj_pval < 0.05, na.rm = TRUE), "\n")
  
  if ("conf_level" %in% colnames(test_results)) {
    cat("Confidence level:", unique(test_results$conf_level), "%\n")
  }
  
  effect_sizes <- test_results$test_statistic_median
  cat("Effect size range:", round(min(effect_sizes, na.rm = TRUE), 4), "to", 
      round(max(effect_sizes, na.rm = TRUE), 4), "\n")
}

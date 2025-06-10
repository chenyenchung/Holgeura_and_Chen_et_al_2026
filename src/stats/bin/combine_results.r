#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

#' Combine multiple depth statistics result files
#' @param pattern Pattern to match result files (default "^results_.*\\.csv$")
#' @param output_file Output filename for combined results
#' @param summary_file Output filename for statistical summary
combine_depth_stats_results <- function(pattern = "^results_.*\\.csv$",
                                       output_file = "combined_depth_stats_results.csv",
                                       summary_file = "statistical_summary.txt") {
  
  # Find all result files
  result_files <- list.files(pattern = pattern)
  
  if (length(result_files) > 0) {
    # Combine all results
    all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)
    
    # Write combined results
    fwrite(all_results, output_file)
    
    # Generate summary statistics
    sink(summary_file)
    cat("Depth Statistics (Wasserstein + KS) Summary\n")
    cat("============================================\n\n")
    
    cat("Total comparisons:", nrow(all_results), "\n")
    cat("Unique neuropils:", length(unique(all_results$neuropil)), "\n")
    cat("Unique presets:", length(unique(all_results$preset)), "\n")
    cat("Confidence level:", unique(all_results$conf_level)[1], "%\n")
    cat("\n")
    
    # Wasserstein effect size summary
    cat("Wasserstein effect size summary:\n")
    cat("  Min:", round(min(all_results$test_statistic_median, na.rm = TRUE), 4), "\n")
    cat("  Q1: ", round(quantile(all_results$test_statistic_median, 0.25, na.rm = TRUE), 4), "\n")
    cat("  Median:", round(median(all_results$test_statistic_median, na.rm = TRUE), 4), "\n")
    cat("  Q3: ", round(quantile(all_results$test_statistic_median, 0.75, na.rm = TRUE), 4), "\n")
    cat("  Max:", round(max(all_results$test_statistic_median, na.rm = TRUE), 4), "\n\n")
    
    # KS p-value summary if available
    if ("ks_pval_raw_median" %in% colnames(all_results)) {
      ks_pvals <- all_results$ks_pval_raw_median
      valid_ks <- ks_pvals[!is.na(ks_pvals)]
      if (length(valid_ks) > 0) {
        cat("KS raw p-value summary:\n")
        cat("  Min:", round(min(valid_ks), 4), "\n")
        cat("  Q1: ", round(quantile(valid_ks, 0.25), 4), "\n")
        cat("  Median:", round(median(valid_ks), 4), "\n")
        cat("  Q3: ", round(quantile(valid_ks, 0.75), 4), "\n")
        cat("  Max:", round(max(valid_ks), 4), "\n")
        cat("KS subsample size:", unique(all_results$ks_subsample_size)[1], "\n")
        cat("KS iterations:", unique(all_results$ks_n_iterations)[1], "\n\n")
      }
    }
    
    # Top effect sizes
    cat("Top 10 largest Wasserstein effect sizes:\n")
    top_effects <- all_results[order(-test_statistic_median)][1:min(10, nrow(all_results))]
    for (i in 1:nrow(top_effects)) {
      cat(sprintf("  %s vs %s (Wasserstein = %.3f, CI: %.3f-%.3f)\n", 
                 top_effects$group1[i], top_effects$group2[i], 
                 top_effects$test_statistic_median[i], 
                 top_effects$test_statistic_lower[i], 
                 top_effects$test_statistic_upper[i]))
    }
    
    sink()
    
    cat("Combined", nrow(all_results), "results into", output_file, "\n")
    cat("Summary written to", summary_file, "\n")
  } else {
    # Create empty files if no results
    write.csv(data.frame(), output_file, row.names = FALSE)
    writeLines("No results to combine", summary_file)
    cat("No result files found matching pattern:", pattern, "\n")
  }
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Set defaults
pattern <- if (is.null(argvs$pattern)) "^results_.*\\.csv$" else argvs$pattern
output_file <- if (is.null(argvs$output)) "combined_depth_stats_results.csv" else argvs$output
summary_file <- if (is.null(argvs$summary)) "statistical_summary.txt" else argvs$summary

# Run the combination
combine_depth_stats_results(pattern, output_file, summary_file)
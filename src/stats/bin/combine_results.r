#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))

# Named vector for converting preset names to readable sheet names
preset_names <- c(
  "subsystem_known" = "Subsystem Known",
  "temporal_known" = "Temporal Known"
)

#' Combine multiple depth statistics result files into Excel workbook
#' @param pattern Pattern to match result files (default "^.*\\.csv$")
#' @param summary_file Output filename for statistical summary
#' @param excel_file Output filename for Excel file with multiple sheets
combine_depth_stats_results <- function(pattern = "^.*\\.csv$",
                                       summary_file = "statistical_summary.txt",
                                       excel_file = "combined_depth_stats_results.xlsx") {
  
  # Find all result files
  result_files <- list.files(pattern = pattern)
  
  if (length(result_files) > 0) {
    # Combine all results
    all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)
    
    # Export to Excel with multiple sheets split by preset
    if (!is.null(excel_file) && excel_file != "") {
      # Split data by preset
      preset_splits <- split(all_results, all_results$preset)
      
      # Create workbook
      wb <- createWorkbook()
      
      # Add each preset as a separate sheet
      for (preset in names(preset_splits)) {
        # Use readable name if available, otherwise use original
        sheet_name <- ifelse(preset %in% names(preset_names), 
                            preset_names[preset], 
                            preset)
        
        # Ensure sheet name is valid (Excel has limitations)
        sheet_name <- substr(gsub("[^[:alnum:][:space:]]", "_", sheet_name), 1, 31)
        
        # Add worksheet
        addWorksheet(wb, sheet_name)
        
        # Write data
        writeData(wb, sheet_name, preset_splits[[preset]])
        
        # Auto-size columns for better readability
        setColWidths(wb, sheet_name, cols = 1:ncol(preset_splits[[preset]]), widths = "auto")
      }
      
      # Save workbook
      saveWorkbook(wb, excel_file, overwrite = TRUE)
    }
    
    # Generate summary statistics
    sink(summary_file)
    
    # Determine analysis type based on columns
    is_broad_depth <- "p_value_exceeds_threshold" %in% colnames(all_results)
    
    if (is_broad_depth) {
      cat("Broad Depth Analysis Summary\n")
      cat("===========================\n\n")
      
      cat("Total analyses:", nrow(all_results), "\n")
      cat("Unique neuropils:", length(unique(all_results$neuropil)), "\n")
      cat("Unique presets:", length(unique(all_results$preset)), "\n")
      cat("Coefficient:", unique(all_results$coefficient)[1], "\n")
      cat("Bootstrap iterations:", unique(all_results$n_bootstrap)[1], "\n")
      cat("Confidence level:", unique(all_results$conf_level)[1], "%\n")
      cat("\n")
      
      # Direction summary
      direction_counts <- table(all_results$direction)
      cat("Direction summary:\n")
      for (dir in names(direction_counts)) {
        cat("  ", dir, ":", direction_counts[dir], "\n")
      }
      cat("\n")
      
      # P-value summary
      if ("p_value_exceeds_threshold" %in% colnames(all_results)) {
        p_vals <- all_results$p_value_exceeds_threshold
        valid_p <- p_vals[!is.na(p_vals)]
        if (length(valid_p) > 0) {
          cat("Raw P-value summary:\n")
          cat("  Min:", round(min(valid_p), 4), "\n")
          cat("  Median:", round(median(valid_p), 4), "\n")
          cat("  Max:", round(max(valid_p), 4), "\n")
          cat("  Significant (raw p < 0.05):", sum(valid_p < 0.05), "/", length(valid_p), "\n\n")
        }
      }
      
      # FDR corrected p-value summary
      if ("p_value_fdr" %in% colnames(all_results)) {
        fdr_vals <- all_results$p_value_fdr
        valid_fdr <- fdr_vals[!is.na(fdr_vals)]
        if (length(valid_fdr) > 0) {
          cat("FDR corrected P-value summary:\n")
          cat("  Min:", round(min(valid_fdr), 4), "\n")
          cat("  Median:", round(median(valid_fdr), 4), "\n")
          cat("  Max:", round(max(valid_fdr), 4), "\n")
          cat("  Significant (FDR p < 0.05):", sum(valid_fdr < 0.05), "/", length(valid_fdr), "\n\n")
        }
      }
      
      # Top significant results (prioritize FDR corrected if available)
      if ("p_value_fdr" %in% colnames(all_results)) {
        sig_results <- all_results[all_results$p_value_fdr < 0.05 & !is.na(all_results$p_value_fdr), ]
        p_col <- "p_value_fdr"
        p_label <- "FDR p"
      } else {
        sig_results <- all_results[all_results$p_value_exceeds_threshold < 0.05, ]
        p_col <- "p_value_exceeds_threshold"
        p_label <- "raw p"
      }
      
      if (nrow(sig_results) > 0) {
        cat(sprintf("Significant results (%s < 0.05):\n", p_label))
        sig_results <- sig_results[order(sig_results[[p_col]]), ]  # Sort by p-value
        for (i in 1:min(10, nrow(sig_results))) {
          cat(sprintf("  %s (%s = %.4f, direction = %s, distance_diff = %.3f)\n", 
                     sig_results$types_of_interest[i], 
                     p_label,
                     sig_results[[p_col]][i],
                     sig_results$direction[i],
                     sig_results$observed_distance_diff[i]))
        }
        cat("\n")
      }
      
    } else {
      cat("Depth Statistics (Wasserstein + KS) Summary\n")
      cat("============================================\n\n")
      
      cat("Total comparisons:", nrow(all_results), "\n")
      cat("Unique neuropils:", length(unique(all_results$neuropil)), "\n")
      cat("Unique presets:", length(unique(all_results$preset)), "\n")
      cat("Confidence level:", unique(all_results$conf_level)[1], "%\n")
      cat("\n")
      
      # Wasserstein effect size summary
      # Check for both naming conventions (wasserstein_* and test_statistic_*)
      median_col <- if ("wasserstein_median" %in% colnames(all_results)) {
        "wasserstein_median"
      } else if ("test_statistic_median" %in% colnames(all_results)) {
        "test_statistic_median"
      } else {
        NULL
      }
      
      if (!is.null(median_col)) {
        cat("Wasserstein effect size summary:\n")
        cat("  Min:", round(min(all_results[[median_col]], na.rm = TRUE), 4), "\n")
        cat("  Q1: ", round(quantile(all_results[[median_col]], 0.25, na.rm = TRUE), 4), "\n")
        cat("  Median:", round(median(all_results[[median_col]], na.rm = TRUE), 4), "\n")
        cat("  Q3: ", round(quantile(all_results[[median_col]], 0.75, na.rm = TRUE), 4), "\n")
        cat("  Max:", round(max(all_results[[median_col]], na.rm = TRUE), 4), "\n\n")
      }
      
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
      # Determine column names based on what's available
      median_col <- if ("wasserstein_median" %in% colnames(all_results)) {
        "wasserstein_median"
      } else if ("test_statistic_median" %in% colnames(all_results)) {
        "test_statistic_median"
      } else {
        NULL
      }
      
      lower_col <- if ("wasserstein_lower" %in% colnames(all_results)) {
        "wasserstein_lower"
      } else if ("test_statistic_lower" %in% colnames(all_results)) {
        "test_statistic_lower"
      } else {
        NULL
      }
      
      upper_col <- if ("wasserstein_upper" %in% colnames(all_results)) {
        "wasserstein_upper"
      } else if ("test_statistic_upper" %in% colnames(all_results)) {
        "test_statistic_upper"
      } else {
        NULL
      }
      
      if (!is.null(median_col) && !is.null(lower_col) && !is.null(upper_col)) {
        cat("Top 10 largest Wasserstein effect sizes:\n")
        top_effects <- all_results[order(-all_results[[median_col]])][1:min(10, nrow(all_results))]
        for (i in 1:nrow(top_effects)) {
          cat(sprintf("  %s vs %s (Wasserstein = %.3f, CI: %.3f-%.3f)\n", 
                     top_effects$group1[i], top_effects$group2[i], 
                     top_effects[[median_col]][i], 
                     top_effects[[lower_col]][i], 
                     top_effects[[upper_col]][i]))
        }
      }
    }
    
    sink()
    
    cat("Combined", nrow(all_results), "results\n")
    cat("Excel file written to", excel_file, "with", length(preset_splits), "sheets\n")
    cat("Summary written to", summary_file, "\n")
  } else {
    # Create empty files if no results
    writeLines("No results to combine", summary_file)
    cat("No result files found matching pattern:", pattern, "\n")
  }
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Set defaults
pattern <- if (is.null(argvs$pattern)) "_depth_stats\\.csv$" else argvs$pattern
summary_file <- if (is.null(argvs$summary)) "statistical_summary.txt" else argvs$summary
excel_file <- if (is.null(argvs$excel)) "combined_depth_stats_results.xlsx" else argvs$excel

# Run the combination
combine_depth_stats_results(pattern, summary_file, excel_file)

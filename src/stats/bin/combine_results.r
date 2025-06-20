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
#' @param pattern Pattern to match result files (default "^results_.*\\.csv$")
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

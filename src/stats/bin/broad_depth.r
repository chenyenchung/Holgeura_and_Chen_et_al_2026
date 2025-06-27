#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rcpp))

perform_broad_depth_analysis <- function(data, ref_superficial, ref_deep, 
                                       types_of_interest, coefficient = 0.25,
                                       n_bootstrap = 1000, conf_int = 95,
                                       threshold_type = "variable", seed = NULL) {
  
  ref_sup_data <- data[
    data$pre_type %in% ref_superficial | data$post_type %in% ref_superficial
  ]
  ref_deep_data <- data[
    data$pre_type %in% ref_deep | data$post_type %in% ref_deep
  ]
  interest_data <- data[data$group %in% types_of_interest]
  
  if (nrow(ref_sup_data) == 0 || nrow(ref_deep_data) == 0 || nrow(interest_data) == 0) {
    warning("Insufficient data for one or more groups")
    return(NULL)
  }
  
  # Extract depth values
  ref_sup_depths <- ref_sup_data$depth
  ref_deep_depths <- ref_deep_data$depth
  interest_depths <- interest_data$depth
  
  # Call C++ bootstrap function
  bootstrap_result <- perform_broad_depth_bootstrap(
    ref_sup_depths = ref_sup_depths,
    ref_deep_depths = ref_deep_depths,
    interest_depths = interest_depths,
    coefficient = coefficient,
    n_bootstrap = n_bootstrap,
    conf_int = conf_int,
    threshold_type = threshold_type,
    seed = seed
  )
  
  # Determine observed direction based on bootstrap results
  if (abs(bootstrap_result$observed_distance_diff) > bootstrap_result$observed_delta_thres) {
    observed_direction <- if (bootstrap_result$observed_distance_diff > 0) "superficial" else "deep"
  } else {
    observed_direction <- "neither"
  }
  
  result <- list(
    types_of_interest = paste(types_of_interest, collapse = ";"),
    ref_superficial = paste(ref_superficial, collapse = ";"),
    ref_deep = paste(ref_deep, collapse = ";"),
    n_types_interest = length(types_of_interest),
    n_ref_superficial = length(ref_superficial),
    n_ref_deep = length(ref_deep),
    n_samples_interest = nrow(interest_data),
    n_samples_superficial = nrow(ref_sup_data),
    n_samples_deep = nrow(ref_deep_data),
    observed_sup_d = bootstrap_result$observed_sup_d,
    observed_deep_d = bootstrap_result$observed_deep_d,
    observed_distance_diff = bootstrap_result$observed_distance_diff,
    observed_delta_thres_base = bootstrap_result$observed_delta_thres_base,
    observed_delta_thres = bootstrap_result$observed_delta_thres,
    bootstrap_distance_diff_median = bootstrap_result$bootstrap_distance_diff_median,
    bootstrap_distance_diff_lower = bootstrap_result$bootstrap_distance_diff_lower,
    bootstrap_distance_diff_upper = bootstrap_result$bootstrap_distance_diff_upper,
    bootstrap_delta_thres_base_median = bootstrap_result$bootstrap_delta_thres_base_median,
    bootstrap_delta_thres_base_lower = bootstrap_result$bootstrap_delta_thres_base_lower,
    bootstrap_delta_thres_base_upper = bootstrap_result$bootstrap_delta_thres_base_upper,
    coefficient = coefficient,
    p_closer_to_superficial = bootstrap_result$p_closer_to_superficial,
    p_closer_to_deep = bootstrap_result$p_closer_to_deep,
    direction = observed_direction,
    threshold_type = threshold_type,
    conf_level = conf_int,
    n_bootstrap = n_bootstrap
  )
  
  return(result)
}

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (interactive()) {
  source("src/utils.r")
  argvs$np <- "ME_L"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "temporal_known"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$sparse_limit <- 100L
  argvs$coefficient <- 0.25
  argvs$n_bootstrap <- 1000L
  argvs$conf_int <- 95.0
  argvs$cppsrc <- "src/stats/bin/broad_depth.cpp"
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  argvs$coefficient <- if (is.null(argvs$coefficient)) 0.25 else as.numeric(argvs$coefficient)
  argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
  argvs$conf_int <- if (is.null(argvs$conf_int)) 95.0 else as.numeric(argvs$conf_int)
  syn_path <- argvs$synf
  if (is.null(argvs$cppsrc)) argvs$cppsrc <- "./broad_depth.cpp"
  if (is.null(argvs$preset)) argvs$preset <- "./viz_preset.csv"
  if (is.null(argvs$use_preset)) stop("use_preset parameter is required")
}

# Compile C++ source
if (file.exists(argvs$cppsrc)) {
  Rcpp::sourceCpp(argvs$cppsrc)
} else {
  stop(sprintf("C++ source file not found: %s", argvs$cppsrc))
}

plot_meta_all <- validate_config_file(
  argvs$meta,
  c(
    "neuropil", "x_axis", "y_axis", "axis_1_func", "axis_2_func", "min1",
    "max1", "min2", "max2", "zid", "zinvert", "minz", "maxz", "outlayout",
    "outr1", "outr2", "outd1", "outd2"
  )
)
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) stop(sprintf("No metadata found for neuropil: %s", argvs$np))

preset_all <- validate_config_file(
  argvs$preset,
  c(
    "preset", "palette", "color_guide", "filter_func", "color_by",
    "notch_split", "do_highlight", "hl_type", "hl_col", "hl_val"
  )
)
preset <- preset_all[preset == argvs$use_preset]
if (nrow(preset) == 0) stop(sprintf("No preset found: %s", argvs$use_preset))

# Validate that required functions exist
validate_functions(preset, plot_meta)

opc_anno <- validate_config_file(argvs$ann, c("cell_type"))
np_coord <- fread(syn_path)

# Apply same filtering as visualization
filter_func <- get(preset$filter_func)
np_coord <- filter_type(
  np_coord, syn_type = argvs$syn_type, sparse_limit = argvs$sparse_limit
)
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
  np_coord$notch_ntype <- paste(np_coord$Notch, np_coord$ntype, sep = "_")
  np_coord_list <- split(np_coord, np_coord$notch_ntype, drop = TRUE)
} else {
  np_coord_list <- list(all = np_coord)
}

# Perform broad depth analysis for each split
all_test_results <- list()

for (i in names(np_coord_list)) {
  current_data <- np_coord_list[[i]]
  
  # Skip if insufficient data
  if (nrow(current_data) <= argvs$sparse_limit) {
    next
  }
  
  # Prepare test data using preset grouping
  if (preset$do_highlight) {
    test_data <- current_data[highlight == TRUE, .(
      group = get(preset$color_by),
      pre_type,
      post_type,
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    test_data <- current_data[, .(
      group = get(preset$color_by),
      pre_type,
      post_type,
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  # Apply z-axis inversion if needed
  if (plot_meta$zinvert) {
    test_data$depth <- test_data$depth * -1
  }
  
  # Get all available groups from filtered data
  available_groups <- unique(
    c(test_data$pre_type, test_data$post_type)
  )
  
  # Define reference groups based on neuropil (keeping existing logic for reference groups)
  if (grepl("^ME_", argvs$np)) {
    ref_superficial <- grep("^Dm", available_groups, value = TRUE)
    ref_deep <- grep("^Pm", available_groups, value = TRUE)
  } else if (grepl("^LOP_", argvs$np)) {
    ref_superficial <- grep("T(4|5)(a|b)$", available_groups, value = TRUE)
    ref_deep <- grep("T(4|5)(c|d)$", available_groups, value = TRUE)
  } else if (grepl("^LO_", argvs$np)) {
    ref_superficial <- c("Tm1", "Tm2", "Tm4")
    ref_deep <- c("Tm20", "Tm5a", "Tm5b", "Tm5c", "Tm5d", "Tm5e", "Tm5f")
    # Filter to only include available groups
    ref_superficial <- intersect(ref_superficial, available_groups)
    ref_deep <- intersect(ref_deep, available_groups)
  } else {
    # For other neuropils, use a generic approach or skip analysis
    ref_superficial <- character(0)
    ref_deep <- character(0)
  }
  
  # Get groups to test
  groups_to_test <- unique(test_data$group)
  
  # Perform bootstrap analysis for each group
  if (length(ref_superficial) > 0 && length(ref_deep) > 0 && length(groups_to_test) > 0) {
    split_results <- lapply(groups_to_test, function(type) {
      out <- perform_broad_depth_analysis(
        test_data,
        ref_superficial = ref_superficial,
        ref_deep = ref_deep,
        types_of_interest = type,
        coefficient = argvs$coefficient,
        n_bootstrap = argvs$n_bootstrap,
        conf_int = argvs$conf_int
      ) 
      return(out)
    })
    
    # Combine results for this split
    if (length(split_results) > 0) {
      split_result_df <- do.call(rbind.data.frame, split_results)
      split_result_df$split <- i
      all_test_results[[i]] <- split_result_df
    }
  } else {
    warning(paste("Insufficient reference groups or groups to test for bootstrap analysis in split:", i))
  }
}
# Combine all results from all splits
if (length(all_test_results) > 0) {
  result <- do.call(rbind, all_test_results)
  
  # Add common metadata if not already present
  if (!"neuropil" %in% colnames(result)) {
    result$neuropil <- argvs$np
  }
  if (!"syn_type" %in% colnames(result)) {
    result$syn_type <- argvs$syn_type
  }
  if (!"preset" %in% colnames(result)) {
    result$preset <- argvs$use_preset
  }
} else {
  result <- NULL
}

out_prefix <- paste(
  argvs$np, argvs$syn_type, argvs$use_preset,
  sep = "_"
)
output_file <- paste0(out_prefix, "_broad_depth.csv")

if (!is.null(result) && nrow(result) > 0) {
  
  cat("Broad depth analysis results written to:", output_file, "\n")
  cat("Number of analyses:", nrow(result), "\n")
  cat("Coefficient used:", unique(result$coefficient), "\n")
  cat("Bootstrap iterations:", unique(result$n_bootstrap), "\n")
  cat("Confidence level:", unique(result$conf_level), "%\n")
  
  sig_superficial <- result[result$p_closer_to_superficial < 0.05, ]
  sig_deep <- result[result$p_closer_to_deep < 0.05, ]
  cat("Significant results closer to superficial (p < 0.05):", nrow(sig_superficial), "\n")
  cat("Significant results closer to deep (p < 0.05):", nrow(sig_deep), "\n")
  
} else {
  result <- data.frame(
    types_of_interest = character(0),
    ref_superficial = character(0),
    ref_deep = character(0),
    n_types_interest = integer(0),
    n_ref_superficial = integer(0),
    n_ref_deep = integer(0),
    n_samples_interest = integer(0),
    n_samples_superficial = integer(0),
    n_samples_deep = integer(0),
    observed_sup_d = numeric(0),
    observed_deep_d = numeric(0),
    observed_distance_diff = numeric(0),
    observed_delta_thres_base = numeric(0),
    observed_delta_thres = numeric(0),
    bootstrap_distance_diff_median = numeric(0),
    bootstrap_distance_diff_lower = numeric(0),
    bootstrap_distance_diff_upper = numeric(0),
    bootstrap_delta_thres_base_median = numeric(0),
    bootstrap_delta_thres_base_lower = numeric(0),
    bootstrap_delta_thres_base_upper = numeric(0),
    coefficient = numeric(0),
    p_closer_to_superficial = numeric(0),
    p_closer_to_deep = numeric(0),
    direction = character(0),
    threshold_type = character(0),
    conf_level = numeric(0),
    n_bootstrap = integer(0),
    split = character(0),
    neuropil = character(0),
    syn_type = character(0),
    preset = character(0)
  )
  cat("No valid analyses could be performed\n")
}
fwrite(result, output_file)
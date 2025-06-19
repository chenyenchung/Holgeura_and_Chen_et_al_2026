#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

### TODO: Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$partner_data <- "partner_summary.csv"
  argvs$threshold <- 0.03
  argvs$output_dir <- "neuropil_partners"
} else {
  argvs$threshold <- as.numeric(argvs$threshold)
}

# Hard-coded types of interest
types_of_interest <- c("Tm5f", "Tm8a", "Li06", "TmY10", "TmY11", "Sm03", "Tm5d",
                      "Tm36", "Tm37", "Sm01", "Sm02", "Sm14", "Tm5e", "Tm33")

# Load partner data
cat("Loading partner data...\n")
partner_data <- read.csv(argvs$partner_data, colClasses = "character")
partner_data$number_of_synapses <- as.numeric(partner_data$number_of_synapses)
partner_data$number_of_neurons <- as.numeric(partner_data$number_of_neurons)

# Function to analyze partners for a neuropil using synapse count
analyze_neuropil_partners_synapses <- function(data, neuropil_name, types_of_interest, threshold) {
  # Filter for types of interest as presynaptic partners in this neuropil
  connections <- data %>%
    filter(neuropil == neuropil_name, pre %in% types_of_interest) %>%
    group_by(pre, post) %>%
    summarise(connection_sum = sum(number_of_synapses), .groups = "drop")
  
  # Calculate ratios and filter by threshold
  partner_ratios <- connections %>%
    group_by(pre) %>%
    mutate(connection_ratio = connection_sum / sum(connection_sum)) %>%
    filter(connection_ratio >= threshold)
  
  # Calculate "Others" category
  others_count <- connections %>%
    group_by(pre) %>%
    mutate(connection_ratio = connection_sum / sum(connection_sum)) %>%
    filter(connection_ratio < threshold) %>%
    summarize(n = n(), .groups = "drop")
  
  others_ratios <- partner_ratios %>%
    group_by(pre) %>%
    summarize(connection_ratio = 1 - sum(connection_ratio), .groups = "drop") %>%
    mutate(post = "Others")
  
  # Add count to Others label
  others_ratios <- left_join(others_ratios, others_count, by = "pre") %>%
    mutate(
      post = ifelse(is.na(n), "Others", paste0("Others (", n, ")")),
      n = NULL
    )
  
  # Combine main partners and others
  plot_data <- rbind(
    partner_ratios[, c("pre", "post", "connection_ratio")],
    others_ratios
  )
  
  return(plot_data)
}

# Function to create donut plot
create_donut_plot <- function(plot_data, neuropil_name, types_of_interest, threshold) {
  # Pool all data across types of interest for the neuropil
  pooled_data <- plot_data %>%
    group_by(post) %>%
    summarise(connection_ratio = sum(connection_ratio), .groups = "drop") %>%
    arrange(desc(connection_ratio)) %>%
    mutate(
      ymax = cumsum(connection_ratio),
      ymin = c(0, head(ymax, n = -1)),
      y = (ymin + ymax) / 2
    )
  
  # Create plot
  p <- ggplot(pooled_data, aes(
    ymin = ymin, ymax = ymax, group = post, 
    fill = connection_ratio, xmin = 3, xmax = 12
  )) +
    geom_rect(color = "black") +
    geom_label_repel(
      aes(x = 9, y = y,
          label = paste0(post, ": ", round(connection_ratio * 100, 1), "%")),
      nudge_x = 5, fill = "white", size = 2
    ) +
    scale_x_continuous(limits = c(0, 18)) +
    scale_fill_gradient(low = "lightblue", high = "blue",
                        limits = c(0, NA), labels = scales::percent) +
    labs(
      title = paste0("Main postsynaptic partners in ", neuropil_name),
      subtitle = paste0("For selected cell types: ", paste(types_of_interest, collapse = ", "), 
                       "\n(Partners accounting for < ", threshold * 100, "% are considered others)"),
      fill = "% of connections"
    ) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(p)
}

# Create output directory
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
}

# Process each neuropil
neuropil_names <- c("ME", "LO", "LOP")

for (neuropil_name in neuropil_names) {
  cat(sprintf("Processing %s neuropil...\n", neuropil_name))
  
  # Analyze partners using synapse count
  plot_data <- analyze_neuropil_partners_synapses(partner_data, neuropil_name, types_of_interest, argvs$threshold)
  
  # Create and save plot
  p <- create_donut_plot(plot_data, neuropil_name, types_of_interest, argvs$threshold)
  
  output_file <- file.path(argvs$output_dir, paste0(neuropil_name, "_partners.pdf"))
  ggsave(output_file, plot = p, width = 8, height = 6)
  
  cat(sprintf("Saved %s\n", output_file))
}

cat("Analysis complete!\n")
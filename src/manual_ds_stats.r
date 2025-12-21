library(ggplot2)
library(dplyr)
library(patchwork)
library(openxlsx)
source("src/utils.r")

boot_sheets_paths <- "int/stats/deep_superficial/combined_broad_depth_results.xlsx"
# openxlsx::getSheetNames(boot_sheets_paths)
## [1] "broad_known"     "broad_new"       "Subsystem Known" "subsystem_new"   
##     "temporal_all"    "Temporal Known"  "temporal_new"    "type_putative"
boot <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "Temporal Known")
)

boot$types_of_interest <- factor(
  boot$types_of_interest,
  levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
)

boot$syn_type <- factor(
  boot$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

# General plotting
dsplot <- function(
    stats, neuropil, split, title, bg_alpha = 0.5,
    star_size = 4, star_nudge = 2500,
    ylab = "Deep / Superficial Bias\n(+: Distal / -: Proximal)"
  ) {
  glossary <- c(
    "Notch Off_intrinsic" = "'Notch'^'Off'~'Interneurons'",
    "Notch Off_projection" = "'Notch'^'Off'~'Projection Neurons'",
    "Notch On_intrinsic" = "'Notch'^'On'~'Interneurons'",
    "Notch On_projection" = "'Notch'^'On'~'Projection Neurons'"
  )
  
  np_lut <- c(
    "ME_L" = "'Medulla'~",
    "LOP_L" = "'Lobula Plate'~",
    "LO_L" = "'Lobula'~"
  )
  deep_col <- ifelse(grepl("^ME", neuropil), "#E6EDE8", "#D9D0E3")
  sup_col <- ifelse(grepl("^ME", neuropil), "#D9D0E3", "#E6EDE8")
  
  
  p <- stats |>
    filter(neuropil == {{ neuropil }}, split %in% {{ split }}) %>%
    ggplot(aes(x = types_of_interest, color = types_of_interest)) +
    annotate(
      geom = "rect",
      fill = deep_col,
      ymin = -Inf,
      ymax = 0,
      xmin = 0,
      xmax = 9,
      alpha = bg_alpha
    ) +
    annotate(
      geom = "rect",
      fill = sup_col,
      ymin = 0,
      ymax = Inf,
      xmin = 0,
      xmax = 9,
      alpha = bg_alpha
    ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(aes(
      y = bootstrap_distance_diff_median,
      ymin = bootstrap_distance_diff_lower,
      ymax = bootstrap_distance_diff_upper
    ), size = 0.5) +
    geom_text(
      aes(y = bootstrap_distance_diff_median + star_nudge,
          label = ifelse(.data$significant_fdr, "*", "")),
      color = "black",
      size = star_size
    ) +
    theme(
      panel.background = element_blank(),
      strip.background = element_rect(fill = "transparent", color = "black"),
      strip.text = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold")
    ) +
    guides(color = "none") +
    labs(
      y = ylab,
      title = parse(text = paste0(
        np_lut[[neuropil]],
        ifelse(length(split) == 1, glossary[[split]], "'Projection'~'Neurons'")
        )
      )
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(-3.5e4, 3.5e4)) +
    scale_color_nk2023()
  
  if (length(split) == 1) {
    return(p + facet_grid(~syn_type))
  } else {
    return(p + facet_grid(split ~ syn_type))
  }
}

## Y1A — Medulla NotchOff interneurons: temporal shift
Y1A <- dsplot(
  boot, "ME_L", "Notch Off_intrinsic",
)
Y1A  

# Y1B — Medulla NotchOn projection neurons: presyn temporal effect; postsyn distal
Y1B <- dsplot(
  boot, "ME_L", "Notch On_projection",
)
Y1B

# Y1C — Medulla NotchOff projection neurons
Y1C <- dsplot(
  boot, "ME_L", "Notch Off_projection",
)
Y1C

# Y1D — Lobula NotchOn projection neurons: early superficial → late deep
Y1D <- dsplot(
  boot, "LO_L", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1D

# Y1E — Lobula NotchOff projection neurons: no monotonic targeting
Y1E <- dsplot(
  boot, "LO_L", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1E

# Y1F - Lobula plate TmY neurons: broadly distributed; exception TmY3 postsyn
Y1F <- dsplot(
  boot, "LOP_L", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1F

## New ones
boot_new <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "temporal_new")
)
boot_new$types_of_interest <- factor(
  boot_new$types_of_interest,
  levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
)

boot_new$syn_type <- factor(
  boot_new$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynase", "Postsynapse")
)

# Y1G - New Sm interneurons: weak/ambiguous bias
Y1G <- dsplot(
  boot_new, "ME_L", "Notch Off_intrinsic",
)
Y1G  

# Y1H - New lobula NotchOn Hbn/Opa/Slp types: deep-biased
Y1H <- dsplot(
  boot_new, "LO_L", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1H

# Y1I - Same cohort in medulla: subtly distal-biased
Y1I <- dsplot(
  boot_new, "ME_L", "Notch On_projection"
)
Y1I

# Y1J - New Erm/Ey NotchOff projection neurons in lobula: deep-biased
Y1J <- dsplot(
  boot_new, "LO_L", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1J

# Y1K - Same types in medulla: around serpentine; no prox/dist bias
Y1K <- dsplot(
  boot_new, "ME_L", "Notch Off_projection"
)
Y1K

Y1 <- list(Y1A, Y1B, Y1C, Y1D, Y1E, Y1F, Y1G, Y1H, Y1I, Y1J, Y1K)

Y1_p <- wrap_plots(Y1, ncol = 3) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "int/Supp_fig_Y1.pdf", plot = Y1_p, width = 8.5, height = 11)

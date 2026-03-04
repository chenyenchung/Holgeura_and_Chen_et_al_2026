library(ggplot2)
library(dplyr)
library(patchwork)
library(openxlsx)
library(RColorBrewer)
source("src/utils.r")

opc_anno <- read.csv("data/visual_neurons_anno.csv")
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
    star_size = 4, star_location = 1.15,
    ylab = "Deep / Superficial Bias\n(+: Distal / -: Proximal)"
  ) {
  glossary <- c(
    "Notch Off_intrinsic" = "'Notch'^'Off'~'Interneurons'",
    "Notch Off_projection" = "'Notch'^'Off'~'Projection Neurons'",
    "Notch On_intrinsic" = "'Notch'^'On'~'Interneurons'",
    "Notch On_projection" = "'Notch'^'On'~'Projection Neurons'",
    "all" = "'Putative'~'OPC'~'Neurons'"
  )
  
  np_lut <- c(
    "ME_R" = "'Medulla'~",
    "LOP_R" = "'Lobula Plate'~",
    "LO_R" = "'Lobula'~"
  )
  deep_col <- ifelse(grepl("^ME", neuropil), "#E6EDE8", "#D9D0E3")
  sup_col <- ifelse(grepl("^ME", neuropil), "#D9D0E3", "#E6EDE8")
  
  pad_width <- max(
    length(unique(stats$types_of_interest[stats$neuropil == neuropil & stats$split %in% split])) + 1,
    nlevels(stats$types_of_interest) + 1
  )
  
  
  p <- stats |>
    filter(neuropil == {{ neuropil }}, split %in% {{ split }}) %>%
    ggplot(aes(x = types_of_interest, color = types_of_interest)) +
    annotate(
      geom = "rect",
      fill = deep_col,
      ymin = -Inf,
      ymax = -0.5,
      xmin = 0,
      xmax = pad_width,
      alpha = bg_alpha
    ) +
    annotate(
      geom = "rect",
      fill = sup_col,
      ymin = 0.5,
      ymax = Inf,
      xmin = 0,
      xmax = pad_width,
      alpha = bg_alpha
    ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(aes(
      y = bootstrap_bias_ratio_median,
      ymin = bootstrap_bias_ratio_lower,
      ymax = bootstrap_bias_ratio_upper
    ), size = 0.5) +
    geom_text(
      aes(label = ifelse(.data$significant_fdr, "*", "")),
      y = star_location,
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
    scale_y_continuous(limits = c(-1, star_location + 0.1))
  
  uval <- unique(stats$types_of_interest[stats$neuropil == neuropil])
  tws <- c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
  is_temporal <- any(uval %in% tws)
  
  sbs <- c("Color", "Form", "Luminance", "Motion", "Object", "Polarization", "Unannotated")
  is_subsystem <- any(uval %in% sbs)
  
  if (length(split) == 1) {
    if (split != "all") {
      if (is_temporal) {
        p <- p + scale_color_nk2023()
      } else {
        p <- p + scale_color_subsystem()
      }
    } else {
      if (is_subsystem) {
        p <- p + scale_color_subsystem()
      } else {
        p <- p + scale_color_type()
      }
      
    }
    return(p + facet_grid(~syn_type))
  } else {
    if (is_temporal) {
      p <- p + scale_color_nk2023()
    } else {
      p <- p + scale_color_subsystem()
    }
    return(p + facet_grid(split ~ syn_type))
  }
}
## Y1A — Medulla NotchOn interneurons: temporal shift
Y1A <- dsplot(
  boot, "ME_R", "Notch On_intrinsic",
)
Y1A  

## Y1B — Medulla NotchOff interneurons: temporal shift
Y1B <- dsplot(
  boot, "ME_R", "Notch Off_intrinsic",
)
Y1B  

# Y1C — Medulla NotchOn projection neurons: presyn temporal effect; postsyn distal
Y1C <- dsplot(
  boot, "ME_R", "Notch On_projection",
)
Y1C

# Y1D — Medulla NotchOff projection neurons
Y1D <- dsplot(
  boot, "ME_R", "Notch Off_projection",
)
Y1D

# Y1E — Lobula NotchOn projection neurons: early superficial → late deep
Y1E <- dsplot(
  boot, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1E

# Y1F — Lobula NotchOff projection neurons: no monotonic targeting
Y1F <- dsplot(
  boot, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1F

# Y1G - Lobula plate TmY neurons: broadly distributed; exception TmY3 postsyn
Y1G <- dsplot(
  boot, "LOP_R", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1G

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

# Y1H - New Sm interneurons: weak/ambiguous bias
Y1H <- dsplot(
  boot_new, "ME_R", "Notch Off_intrinsic",
)
Y1H  

# Y1I - New lobula NotchOn Hbn/Opa/Slp types: deep-biased
Y1I <- dsplot(
  boot_new, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1I

# Y1J - Same cohort in medulla: subtly distal-biased
Y1J <- dsplot(
  boot_new, "ME_R", "Notch On_projection"
)
Y1J

# Y1K - New Erm/Ey NotchOff projection neurons in lobula: deep-biased
Y1K <- dsplot(
  boot_new, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1K

# Y1L - Same types in medulla: around serpentine; no prox/dist bias
Y1L <- dsplot(
  boot_new, "ME_R", "Notch Off_projection"
)
Y1L

Y1 <- list(Y1A, Y1B, Y1C, Y1D, Y1E, Y1F, Y1G, Y1H, Y1I, Y1J, Y1K, Y1L)

Y1_p <- wrap_plots(Y1, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y1.pdf", plot = Y1_p, width = 8.5, height = 11)

### Y2
boot_p1 <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "type_putative_1")
)

boot_p1$syn_type <- factor(
  boot_p1$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y2A <- dsplot(
  boot_p1, "LO_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y2A

boot_p2 <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "type_putative_2")
)

boot_p2$syn_type <- factor(
  boot_p2$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y2B <- dsplot(
  boot_p2, "LO_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y2B

Y2C <- dsplot(
  boot_p1, "ME_R", "all"
)
Y2C

Y2D <- dsplot(
  boot_p2, "ME_R", "all"
)
Y2D

boot_p3 <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "type_putative_3")
)

boot_p3$syn_type <- factor(
  boot_p3$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y2E <- dsplot(
  boot_p3, "ME_R", "all"
)
Y2E

Y2 <- list(Y2A, Y2B, Y2C, Y2D, Y2E, plot_spacer(), plot_spacer(), plot_spacer())

Y2_p <- wrap_plots(Y2, ncol = 2) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y2.pdf", plot = Y2_p, width = 8.5, height = 11)

### Y3
boot_fun <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "Subsystem Known")
)

boot_fun$syn_type <- factor(
  boot_fun$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y3A <- dsplot(
  boot_fun, "ME_R", "Notch On_intrinsic"
)
Y3A

Y3B <- dsplot(
  boot_fun, "ME_R", "Notch Off_intrinsic"
)
Y3B

Y3C <- dsplot(
  boot_fun, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3C

Y3D <- dsplot(
  boot_fun, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3D

Y3E <- dsplot(
  boot_fun, "ME_R", "Notch On_projection"
)
Y3E

Y3F <- dsplot(
  boot_fun, "ME_R", "Notch Off_projection"
)
Y3F

boot_fun_new <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "subsystem_new")
)

boot_fun_new$syn_type <- factor(
  boot_fun_new$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y3G <- dsplot(
  boot_fun_new, "ME_R", "Notch Off_intrinsic"
)
Y3G

Y3H <- dsplot(
  boot_fun_new, "ME_R", "Notch Off_projection"
)
Y3H

Y3I <- dsplot(
  boot_fun_new, "ME_R", "Notch On_projection"
)
Y3I

Y3J <- dsplot(
  boot_fun_new, "LO_R", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3J

Y3K <- dsplot(
  boot_fun, "LOP_R", c("Notch On_projection", "Notch Off_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)

Y3K

Y3 <- list(Y3A, Y3B, Y3C, Y3D, Y3E, Y3F, Y3G, Y3H, Y3I, Y3J, Y3K)

Y3_p <- wrap_plots(Y3, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y3.pdf", plot = Y3_p, width = 8.5, height = 11)

### Y4
boot_fun_putative <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "subsystem_putative")
)

boot_fun_putative$syn_type <- factor(
  boot_fun_putative$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y4A <- dsplot(
  boot_fun_putative, "ME_R", "all"
)
Y4A

Y4B <- dsplot(
  boot_fun_putative, "LO_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y4B

Y4C <-dsplot(
  boot_fun_putative, "LOP_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
) 
Y4C

Y4 <- list(
  Y4A, Y4B, Y4C,
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer()
)

Y4_p <- wrap_plots(Y4, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y4.pdf", plot = Y4_p, width = 8.5, height = 11)
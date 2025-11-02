anno <- read.csv("data/visual_neurons_anno.csv")
p15tf <- read.csv("data/P15_tf.csv", row.names = 1, check.names = FALSE)
p15cam <- read.csv("data/P15_CAM.csv", row.names = 1, check.names = FALSE)
sel <- read.csv("data/selectors.csv", row.names = 1, check.names = FALSE)
  

# Make conversion table for column names
lut <- paste0(
  anno$cell_type[!is.na(anno$ozel2021_cluster)],
  " (", anno$ozel2021_cluster[!is.na(anno$ozel2021_cluster)], ")"
)
names(lut) <- as.character(anno$ozel2021_cluster[!is.na(anno$ozel2021_cluster)])

colnames(p15tf)[colnames(p15tf) %in% names(lut)] <- lut[colnames(p15tf)[colnames(p15tf) %in% names(lut)]]
colnames(p15cam)[colnames(p15cam) %in% names(lut)] <- lut[colnames(p15cam)[colnames(p15cam) %in% names(lut)]]
colnames(sel)[colnames(sel) %in% names(lut)] <- lut[colnames(sel)[colnames(sel) %in% names(lut)]]

p15tfh <- p15tf
p15tfh[as.matrix(p15tf)] <- "On"
p15tfh[!as.matrix(p15tf)] <- "Off"

p15camh <- p15cam
p15camh[as.matrix(p15cam)] <- "On"
p15camh[!as.matrix(p15cam)] <- "Off"

selh <- sel
selh[as.matrix(sel)] <- "On"
selh[!as.matrix(sel)] <- "Off"

write.csv(p15tfh, "int/P15_TF_readable.csv", quote = FALSE)
write.csv(p15camh, "int/P15_CAM_readable.csv", quote = FALSE)
write.csv(selh, "int/Selector_readable.csv", quote = FALSE)

anno <- read.csv("data/visual_neurons_anno.csv")
anno <- anno[anno$Confident_annotation == "Y", ]

tws <- unique(anno$temporal_label)
funcs <- unique(anno$func)

out <- data.frame(
  temporal = character(56),
  subsystem = character(56),
  p.val = 1,
  or = 1
)

i <- 0
for (tw in tws) {
  for (func in funcs) {
    ta <- table(anno$temporal_label == tw, anno$func == func)
    i <- i + 1
    out$temporal[i] <- tw
    out$subsystem[i] <- func
    out$p.val[i] <- fisher.test(ta)$p.value
    out$or[i] <- fisher.test(ta)$estimate
  }
}
out$adj.p <- p.adjust(out$p.val, method = "BH")

write.csv(out, "int/Sup_tbl_3_per_window_func_enrichment.csv", row.names = FALSE)

# Path for mixture model objects
# The objects were derived from:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142787/suppl/GSE142787%5FMixture%5Fmodeling.xlsx
mmp <- "/scratch/cgsb/desplan/File_exchange/Ozel_and_Simon_2021_related/Mixture_Modeling/"
mmo <- list.files(mmp)

tf <- read.csv("data/FlyBase_TFs.csv")

# Drop a temporary file that I left there before.
mmo <- mmo[-1]

names(mmo) <- vapply(
  mmo, function(x) {
    return(strsplit(x, "_")[[1]][1])
  },
  FUN.VALUE = character(1)
)
obj_l <- lapply(mmo, function(x) {
  obj <- readRDS(file.path(mmp, x))
  return(obj)
})

# Get the genes that are present in all stages
consensus_genes <- lapply(obj_l, rownames)
consensus_genes <- Reduce(intersect, consensus_genes)
consensus_tfs <- intersect(consensus_genes, tf$SYMBOL) |>
  sort()

# Get the clusters that are present in all stages
consensus_clusters <- lapply(obj_l, colnames)
consensus_clusters <- Reduce(intersect, consensus_clusters) |>
  sort()

# Drop the last 9 columns that are QC info
consensus_clusters <- consensus_clusters[seq_len(length(consensus_clusters) - 9)]

# Format the probability matrices to be comparable across stages
obj_formatted <- lapply(obj_l, function(x) {
  x <- x[consensus_tfs, consensus_clusters]
  x <- as.matrix(x)
  class(x) <- "numeric"
  x <- x > 0.5
  return(x)
})

obj_bin <- Reduce(`&`, obj_formatted)
write.csv(obj_bin, "data/selector.csv")

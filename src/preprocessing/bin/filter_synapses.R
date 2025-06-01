#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValue = TRUE)

if (is.null(argvs$syn) || !file.exists(argvs$syn)) {
  stop("Cannot find the feather file containing synases.")
}

if (is.null(argvs$ann) || !file.exists(argvs$ann)) {
  stop("Cannot find the csv file containing proofread synapse counts.")
}

if (is.null(argvs$type) || !file.exists(argvs$type)) {
  stop("Cannot find the csv file containing type annotations.")
}

if (!is.null(argvs$npil)) {
  ol_npils <- strsplit(argvs$npil, ":")[[1]]
} else {
  message(
    "Default optic lobe neuropil labels are used.\n",
    "Please override the labels with --npil and use colon (:) as\n",
    "your delimiter (e.g., --npil ME:LO:LOP)."
  )
  ol_npils <- c("ME_L", "ME_R", "LO_L", "LO_R", "LOP_L", "LOP_R")
}

message(
  "The followings were considered as optic lobe neuropils:\n",
  paste(ol_npils, collapse = ", "), "."
)

message(
  "Filering the synapse table assuming neuropil membership is in\n",
  "a column named 'neuropil'."
)
syn_coord <- read_feather(argvs$syn)
syn_coord <- data.table(syn_coord)
syn_coord <- syn_coord[neuropil %in% ol_npils]


ann_filter <- fread(argvs$ann)

message(
  "Filtering synapse annotation to keep synapses with > 4 counts.\n",
  "Assuming the existence of a neuropil column containing neuropil labels\n",
  "and a syn_count column containing synapse counts to be filtered."
)
ann_filter <- ann_filter[neuropil %in% ol_npils & syn_count > 4]
ann_filter <- ann_filter[ , .(pre_pt_root_id, post_pt_root_id, neuropil)]

syn_coord <- merge(
  syn_coord, ann_filter, by = c("pre_pt_root_id", "post_pt_root_id", "neuropil")
)

message(
  "Annotating cell types involved in each synapse.\n",
  "Assuming the existence of a root_id column containing each root\n",
  "and a primary_type column containing the type labels."
)

type_ann <- fread(argvs$type)
type_ann <- type_ann[, .(root_id, primary_type)]

syn_coord <- merge(syn_coord, type_ann,
                   by.x = "pre_pt_root_id", by.y = "root_id", all.x = TRUE)
setnames(syn_coord, old = "primary_type", new = "pre_type")
syn_coord <- merge(syn_coord, type_ann,
                   by.x = "post_pt_root_id", by.y = "root_id", all.x = TRUE)
setnames(syn_coord, old = "primary_type", new = "post_type")

fwrite(syn_coord, "ol_passed_filter.csv.gz")
library(tidyverse)
library(plyranges)
library(nullranges)
library(furrr)

threads <- ifelse(exists("snakemake"),snakemake@threads,4)
plan(multisession, workers = threads)

lms_path <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
                   "upstream/final-models.collected-info.tsv.gz")

lms <- read_tsv(lms_path) %>% filter(significant_x)


anno_ins_path <- ifelse(exists("snakemake"),snakemake@input[["anno_ins"]],
                        "results/resources/annotated_fixed_insertions.gr.rds")

ins <- read_rds(anno_ins_path)

# 1. get fixed insertions to other script
# 2. For all TFs, and for all TEs
# 3. need bootstraps for sure. = maybe frame as bootstrap score rather than trying to do wilcox afterwards

# add1 is just the first TF, so this is the first index in our slice
tfs <- colnames(mcols(ins))[which(colnames(mcols(ins)) == "ADD1"):ncol(mcols(ins))]

tes <- unique(ins$repeat_element)
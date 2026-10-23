# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# If you use this code, please cite: Dashti-Gibson N., et al. Cross-cancer comparisons of multi-omics data for therapeutic selection.

# Figure 2e-g: evenness of sequencing coverage across the genome, per sample.
#
# Input: one deepTools bedGraph per sample (<sample>_deeptools_rmDupl.bedgraph; columns chr, start, end,
# coverage), the same files as for Supplementary Figure S3:
#   bamCoverage -b <sample>.sorted.bam -bs 1000000 -p 12 -of bedgraph --ignoreDuplicates \
#               -o <sample>_deeptools_rmDupl.bedgraph
#
# For each sample, using the bedGraph rows on chr1-chr22 (zero-coverage rows kept):
#   Gini coefficient = DescTools::Gini(coverage)     0 = perfectly even coverage, 1 = maximally uneven
#   CV               = sd(coverage) / mean(coverage)  coefficient of variation
# Each bedGraph row counts once. (deepTools writes adjacent bins with identical coverage as a single row.)
#
# Output: a table with the columns of Fig2e-g-Gini.tsv (sample, chr, cv, gini_coefficient).

library(data.table)
library(DescTools)

gini_cv <- function(bedgraph) {
  bg <- fread(bedgraph, header = FALSE, col.names = c("chr", "start", "end", "coverage"))
  cov <- bg[chr %in% paste0("chr", 1:22), coverage]
  data.table(sample           = sub("_deeptools_rmDupl\\.bedgraph$", "", basename(bedgraph)),
             chr              = "autosomes",
             cv               = sd(cov) / mean(cov),
             gini_coefficient = Gini(cov))
}

# Example: all bedGraphs in one folder
files <- list.files("bedgraphs", pattern = "_deeptools_rmDupl\\.bedgraph$", full.names = TRUE)
gini  <- rbindlist(lapply(files, gini_cv))
fwrite(gini, "Gini_autosomes.tsv", sep = "\t")

# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# If you use this code, please cite: Dashti-Gibson N., et al. Pan-cancer comparisons of multi-omics data for therapeutic selection.

# Figure 2b-c: high-confidence somatic variant density and cancer-associated genes (CAGs) with somatic variants, per model.
#
# Input: classified_results.2026-04-15.RData, the output of classifier ver0.1.0
# (FigS2/classifier_v0.1.0.R, classify_variants_mutect()) run on each model's Mutect2
# tumour-only calls. It holds one list element per model with a $somatic table (variants classified somatic,
# functional_uncertain_somatic or likely_somatic).
#
# Only variants with classification == "somatic" are counted (not functional_uncertain_somatic or likely_somatic).
#
# High-confidence somatic variant density (Figure 2b), in somatic variants per Mb
#   = number of somatic SNVs, insertions and deletions (Variant_Type SNP, INS, DEL) on chr1-chr22
#     / total length of chr1-chr22 in Mb (GRCh38; 2,875 Mb)
#   "Somatic" here is the classifier's highest-confidence tier: known hotspots, or variants in genes with > 100 COSMIC
#   alterations whose allele fraction is inconsistent with germline. The value is therefore not a genome-wide mutation
#   rate and not equivalent to clinical tumour mutational burden (TMB), which counts coding variants in a gene panel or exome.
#
# Cancer-associated genes (Figure 2c)
#   Genes with at least one somatic variant for which cancer_associated_gene == "YES" (gene in the OncoKB cancer
#   gene list, protein-altering or splice-site change, population frequency < 0.1 %), excluding intronic, UTR,
#   5' flank and silent variants. Each gene is counted once per model.
#
# Outputs: Fig2B-somatic_variant_density.tsv, Fig2C-CAGs.tsv

library(dplyr)
library(tibble)

load("classified_results.2026-04-15.RData")

# ---- High-confidence somatic variant density (Figure 2b) ---------------------------------------------
chrom_lengths <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555, chr5 = 181538259,
  chr6 = 170805979, chr7 = 159345973, chr8 = 145138636, chr9 = 138394717, chr10 = 133797422,
  chr11 = 135086622, chr12 = 133275309, chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
  chr16 = 90338345, chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468
)
autosomes <- names(chrom_lengths)
autosome_size_mb <- sum(chrom_lengths) / 1e6

density <- bind_rows(lapply(names(classified_results), function(sample) {
  v <- classified_results[[sample]]$somatic %>%
    filter(Chromosome %in% autosomes,
           classification == "somatic",
           Variant_Type %in% c("SNP", "INS", "DEL"))
  tibble(Sample = sample, n_somatic_variants = nrow(v), somatic_variants_per_Mb = nrow(v) / autosome_size_mb)
})) %>%
  mutate(Cancer_Type = sapply(strsplit(Sample, "-"), `[`, 2), .after = Sample)

write.table(density, "Fig2B-somatic_variant_density.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Cancer-associated genes (Figure 2c) --------------------------------------------------------------
cag <- lapply(names(classified_results), function(sample) {
  classified_results[[sample]]$somatic %>%
    filter(classification == "somatic") %>%
    filter(!Variant_Classification %in% c("Intron", "3'UTR", "5'Flank", "5'UTR", "Silent")) %>%
    filter(cancer_associated_gene == "YES") %>%
    pull(Hugo_Symbol) %>%
    unique()
})
names(cag) <- names(classified_results)

cag_table <- tibble(
  Sample = names(cag),
  Cancer_Type = sapply(strsplit(names(cag), "-"), `[`, 2),
  n_cancer_associated_genes = lengths(cag)
)

write.table(cag_table, "Fig2C-CAGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

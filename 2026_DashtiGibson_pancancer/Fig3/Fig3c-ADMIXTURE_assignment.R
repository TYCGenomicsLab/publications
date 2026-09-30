# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# If you use this code, please cite: Dashti-Gibson N., et al. Pan-cancer comparisons of multi-omics data for therapeutic selection.

# Figure 3c / Suppl Table: assign each PDX model to a 1000 Genomes Project (1KGP) superpopulation and population
# from its ADMIXTURE ancestry proportions.
#
# Inputs
#   pruned_merged_data.5.Q   ADMIXTURE output, K = 5: one row per sample, five ancestry proportions
#   pruned_merged_data.fam   PLINK sample file for the same samples, in the same order (1KGP samples + PDX models,
#                            listed by their model IDs)
#   20130606_g1k_3202_samples_ped_population.txt   1KGP sample labels (SampleID, Population, Superpopulation)
#
# Method
#   1. Reference profiles: for every 1KGP superpopulation (5) and population (26), the mean of the five
#      ancestry proportions over its 1KGP samples.
#   2. Each ADMIXTURE component is named after the superpopulation with the highest mean in it (AFR, AMR, EUR,
#      EAS, SAS).
#   3. Each PDX model is assigned the superpopulation, and separately the population, whose reference profile
#      has the highest Pearson correlation with the model's five proportions.
#   4. The reported superpopulation is the one the assigned population belongs to. Where this differs from
#      the superpopulation found in step 3, the model is flagged for manual review (PCA position, ADMIXTURE
#      bar). In this cohort only VCU-PC-065 is flagged: population PUR, so superpopulation AMR (confirmed).
#
# Outputs
#   Fig3c-ADMIXTURE.tsv              reference profiles, then the PDX models
#   ADMIXTURE_assignment_check.tsv   per model: both calls, best and second-best correlations, needs_review

library(data.table)

K <- paste0("Population", 1:5)

# ---- Load ---------------------------------------------------------------------------------------------
q   <- setnames(fread("pruned_merged_data.5.Q", header = FALSE), K)
fam <- setnames(fread("pruned_merged_data.fam", header = FALSE),
                c("FamilyID", "SampleID", "PaternalID", "MaternalID", "Sex", "Phenotype"))
kg  <- fread("20130606_g1k_3202_samples_ped_population.txt")[, .(SampleID, Population, Superpopulation)]
dat <- merge(cbind(fam[, .(FamilyID, SampleID)], q), kg, by = "SampleID", all.x = TRUE, sort = FALSE)

# ---- 1. Reference profiles from the 1KGP samples ------------------------------------------------------
sup_means <- dat[!is.na(Superpopulation), lapply(.SD, mean), by = Superpopulation, .SDcols = K]
pop_means <- dat[!is.na(Population),      lapply(.SD, mean), by = .(Population, Superpopulation), .SDcols = K]

# ---- 2. Name the ADMIXTURE components -----------------------------------------------------------------
comp_names <- sapply(K, function(k) sup_means$Superpopulation[which.max(sup_means[[k]])])
stopifnot(!anyDuplicated(comp_names))

# ---- 3. Best-correlated reference profile -------------------------------------------------------------
# Returns the best and second-best reference label and their correlations.
best_match <- function(x, means, label) {
  r <- apply(as.matrix(means[, ..K]), 1, function(profile) cor(x, profile, method = "pearson"))
  o <- order(r, decreasing = TRUE)[1:2]
  list(means[[label]][o[1]], r[o[1]], means[[label]][o[2]], r[o[2]])
}
pdx <- dat[grepl("^VCU-", FamilyID)]

# Three sequenced samples that are not part of the published cohort of 60 models.
not_in_cohort <- c("VCU-CO-078", "VCU-PC-070", "VCU-PC-095")
pdx <- pdx[!SampleID %in% not_in_cohort]
stopifnot(!anyDuplicated(pdx$SampleID))
for (i in seq_len(nrow(pdx))) {
  x <- unlist(pdx[i, ..K])
  s <- best_match(x, sup_means, "Superpopulation")
  p <- best_match(x, pop_means, "Population")
  set(pdx, i, c("sup_by_correlation", "sup_r", "sup_2nd", "sup_2nd_r"), s)
  set(pdx, i, c("AssignedPopulation", "pop_r", "pop_2nd", "pop_2nd_r"), p)
}

# ---- 4. Make the two levels agree; flag the models that need a manual check ---------------------------
# The superpopulation and the population are assigned independently, so they can disagree, most often for
# admixed genomes (e.g. a mainly European genome whose closest population is the admixed PUR). The reported
# superpopulation is the one the assigned population belongs to; every disagreement is listed for review.
pdx[, sup_of_population := pop_means$Superpopulation[match(AssignedPopulation, pop_means$Population)]]
pdx[, needs_review := sup_by_correlation != sup_of_population]
pdx[, AssignedSuperpopulation := sup_of_population]

check <- pdx[, .(SampleID, AssignedSuperpopulation, AssignedPopulation, sup_by_correlation, sup_r, sup_2nd,
                 sup_2nd_r, pop_r, pop_2nd, pop_2nd_r, needs_review)]
fwrite(check, "ADMIXTURE_assignment_check.tsv", sep = "\t", quote = FALSE)

if (any(pdx$needs_review)) {
  message(sum(pdx$needs_review), " model(s) where the superpopulation and population calls disagree.\n",
          "AssignedSuperpopulation was set from the assigned population; verify each one manually:\n",
          "  - PCA: does the model sit among the 1KGP samples of that superpopulation, or between clusters?\n",
          "  - ADMIXTURE bar: do its five proportions look like the reference profile of the assigned\n",
          "    population, or like a mixture of superpopulations?\n",
          "Correlations for all models are in ADMIXTURE_assignment_check.tsv.")
  print(check[needs_review == TRUE, .(SampleID, AssignedPopulation, population_superpop = AssignedSuperpopulation,
                                      superpop_by_correlation = sup_by_correlation, sup_r, pop_r)])
}

# ---- Output -------------------------------------------------------------------------------------------
out <- rbind(
  sup_means[, .(SampleID = paste0(Superpopulation, "_REF_PROFILE"), AssignedSuperpopulation = Superpopulation,
                AssignedPopulation = NA_character_, .SD), .SDcols = K],
  pop_means[, .(SampleID = paste0(Population, "_REF_PROFILE"), AssignedSuperpopulation = NA_character_,
                AssignedPopulation = Population, .SD), .SDcols = K],
  pdx[, c("SampleID", "AssignedSuperpopulation", "AssignedPopulation", K), with = FALSE],
  use.names = FALSE)
setnames(out, c("SampleID", "AssignedSuperpopulation", "AssignedPopulation", comp_names))
fwrite(out, "Fig3c-ADMIXTURE.tsv", sep = "\t", na = "NA", quote = FALSE)

# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# If you use this code, please cite: Dashti-Gibson N., et al. Cross-cancer comparisons of multi-omics data for therapeutic selection.

# Figure 3a-b: principal component analysis (PCA) of the VCU PDX models together with the 1000 Genomes Project (1KGP)
# samples, PC1 vs PC2.
#   3a  points coloured by 1KGP superpopulation; PDX models in dark grey
#   3b  points coloured by 1KGP population, one label per population
#
# Inputs (PLINK --pca on the jointly genotyped, LD-pruned 1KGP + PDX data)
#   pruned_pca_results.eigenvec   sample IDs and PC1-PC20 for 3,202 1KGP samples + 63 VCU PDX samples (model IDs)
#   pruned_pca_results.eigenval   eigenvalues of PC1-PC20; axis percentages = eigenvalue / sum of the 20 eigenvalues
#   20130606_g1k_3202_samples_ped_population.txt   1KGP sample labels (SampleID, Population, Superpopulation)
#
# Outputs: Fig3a_PCA_superpopulations.svg, Fig3b_PCA_populations.svg

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# Sequenced samples that are not part of the published cohort of 60 models.
EXCLUDE <- c("VCU-CO-078", "VCU-PC-070", "VCU-PC-095", "VCU-PC-067", "VCU-PC-133")
 
# ---- Read PCA and sample labels -----------------------------------------------------------------------
metadata <- read.table("20130606_g1k_3202_samples_ped_population.txt", header = TRUE, stringsAsFactors = FALSE)
pca_data <- read.table("pruned_pca_results.eigenvec", header = FALSE, comment.char = "")
colnames(pca_data) <- c("FamilyID", "SampleID", paste0("PC", 1:20))

pca_data <- pca_data[!pca_data$SampleID %in% EXCLUDE, ]
pca_data <- dplyr::select(pca_data, -FamilyID)

eigenval <- scan("pruned_pca_results.eigenval")
percent_var <- 100 * eigenval / sum(eigenval)

pca_data_join <- pca_data %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(FID = ifelse(!is.na(Superpopulation), Superpopulation, SampleID),
         Group = ifelse(grepl("^VCU-", FID), "VCU", FID),
         PopulationGroup = ifelse(grepl("^VCU-", FID), "VCU", Population),
         Alpha = ifelse(grepl("^VCU-", FID), 1.0, 0.1))   # PDX models opaque, 1KGP samples transparent

axis_labs <- labs(x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
                  y = paste0("PC2 (", round(percent_var[2], 1), "%)"))

# ---- 3a: superpopulations -----------------------------------------------------------------------------
non_vcu <- setdiff(unique(pca_data_join$Group), "VCU")
group_color_map <- setNames(c(colorRampPalette(brewer.pal(8, "Dark2"))(length(non_vcu)), "grey20"), c(non_vcu, "VCU"))

ggplot(pca_data_join, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group, alpha = Alpha), size = 2) +
  scale_color_manual(values = group_color_map) +
  scale_alpha_identity() +
  theme_minimal() +
  axis_labs + labs(color = "Population")
ggsave("Fig3a_PCA_superpopulations.svg", width = 6.5, height = 5, bg = "white")

# ---- 3b: populations, one label per population --------------------------------------------------------
non_vcu <- setdiff(unique(pca_data_join$PopulationGroup), "VCU")
color_map <- setNames(c(colorRampPalette(brewer.pal(8, "Dark2"))(length(non_vcu)), "grey20"), c(non_vcu, "VCU"))

highlight_points <- pca_data_join %>%
  filter(!is.na(Population)) %>%
  group_by(Population) %>%
  slice_head(n = 1) %>%          # the first sample of each population carries its label
  ungroup()

set.seed(1)                      # label placement by ggrepel is random; fixed for reproducibility
ggplot(pca_data_join, aes(x = PC1, y = PC2, color = PopulationGroup)) +
  geom_point(aes(alpha = Alpha), size = 2) +
  geom_text_repel(data = highlight_points, aes(label = Population, color = PopulationGroup), show.legend = FALSE,
                  size = 3, box.padding = 0.4, point.padding = 0.2, segment.size = 0.5, segment.curvature = 0.1,
                  segment.ncp = 3, segment.angle = 20, force = 2, force_pull = 1, min.segment.length = 0,
                  max.overlaps = Inf) +
  scale_color_manual(values = color_map) +
  scale_alpha_identity() +
  guides(alpha = "none") +
  theme_minimal() +
  axis_labs + labs(color = "Population") +
  theme(legend.position = "right")
ggsave("Fig3b_PCA_populations.svg", width = 7, height = 5, bg = "white")

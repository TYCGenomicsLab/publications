# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# If you use this code, please cite: Dashti-Gibson N., et al. Cross-cancer comparisons of multi-omics data for therapeutic selection.

# Supplementary Figure S3: genome-wide read coverage of one sample in 1-Mb bins.
# Usage: Rscript plot_genomewide_coverage.R <folder> <sample>_deeptools_rmDupl.bedgraph
# Output: <sample>.coverage.png in <folder>

# Load necessary libraries silently
suppressPackageStartupMessages({
  library(gtools)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to variables
WRD_ <- args[1]
file_name <- args[2]

# Keep everything until the first underscore
sample <- sub("_.*$", "", file_name)

print(sample)

# Set working directory
setwd(WRD_)

# Read the coverage data
coverage_data <- fread(file_name, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "coverage"))

# Keep main chromosomes only
coverage_data = coverage_data[grep("chr", coverage_data$chr),]
coverage_data = coverage_data[-grep("chrM", coverage_data$chr),]

# Prepare data for plotting
coverage_data[, midpoint := (start + end) / 2]

# Express coverage as approx. depth
##coverage_data[, depth := (coverage)*150 / 10^6]

# Convert chr to a factor to ensure alphanumeric ordering
coverage_data$chr <- factor(coverage_data$chr, levels = mixedsort(unique(coverage_data$chr)))

# Plotting all chromosomes together, using facet:
# Based on https://stackoverflow.com/questions/55717945/plotting-coverage-depth-in-1kb-windows

# Calculate mean coverage across all chromosomes
mean_coverage <- mean(coverage_data$coverage, na.rm = TRUE)

coverage_data$chr_num <- gsub("chr", "", coverage_data$chr)
coverage_data$chr_num[which(coverage_data$chr_num=="X")] <- 23
coverage_data$chr_num[which(coverage_data$chr_num=="Y")] <- 24
coverage_data$chr_num = as.numeric(coverage_data$chr_num)
coverage_data$fill <- ifelse(coverage_data$chr_num %% 2 == 0, "groupA", "groupB")
#coverage_data$fill <- (coverage_data$chr_num %% 2) + 1 # needs to be discrete values. This is red as continous.

coverage_plot <- ggplot(coverage_data, aes(x = midpoint, y = coverage)) +
  geom_rect( aes(fill = fill),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,show.legend =FALSE) +
  geom_point(size = 1, alpha = 0.5) +
  geom_hline(yintercept = mean_coverage, color = "red", linetype = "dashed") +  # Add red dashed line for mean coverage
  geom_hline(yintercept = 2*mean_coverage, color = "darkred", linetype = "dashed") +  # Add red dashed line for mean coverage
  geom_hline(yintercept = 0, color = "gold2", linetype = "dashed") +  # Add red dashed line for mean coverage
  facet_grid(~ chr_num, scales = "free_x", space = "free_x" ) +  # Use custom labeller
  scale_fill_manual(values = c("white", "grey95")) +
  scale_x_continuous(
    breaks = seq(0, max(coverage_data$end), by = 100000000),  # Set breaks at every 50 Mb
    labels = function(x) x / 10^6,  # Simplify labels to Mb
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 10^6),  # Set y-axis limits
    breaks = seq(0, 10^6, by = 250000),
    labels = function(x) x / 10^6  # Simplify labels to Mb
  ) +
  labs(
    x = "Genomic Position (Mb)",
    y = "Total Reads Mapped (x Million)",
    title = sample
  ) +
  theme(
    aspect.ratio = NULL,
    panel.spacing = unit(0, "mm"),
    #strip.text.x = element_text(face = "bold"),  # Bold the strip text
    panel.background = element_blank()
    #strip.background = element_rect(fill = rep(c("white", "gray90"), length.out = length(unique(coverage_data$chr))))  # Alternating colors
    )

png(paste0(sample,".coverage.png"), width = 15, height = 4, units = "in", res = 300)  # Open a PNG device
print(coverage_plot)  # Draw the ggplot object
dev.off()  # Close the device


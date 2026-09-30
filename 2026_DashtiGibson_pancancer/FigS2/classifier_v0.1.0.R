# Copyright (c) VCU 2026, CC BY-NC 4.0, Katarzyna M. Tyc.
# v0.1.0 2026-04-10
# If you use this code, please cite: Dashti-Gibson N., et al. Pan-cancer comparisons of multi-omics data for therapeutic selection.

library(dplyr)
library(stringr)
library(tidyr)
library(readxl)

oncogenes <- read.table("cancerGeneList.2026-03-27.tsv", header = T, sep = "\t", fill = T, comment.char = "#")
oncogene_set <- unique(oncogenes$Hugo_Symbol)
dim(oncogenes)
# [1] 1231   16

# ---- Variant classification (v0.1.0) ---- #

parse_popAF <- function(x) {
  # In most cases, the smallest non-missing, non-ref value is the one you want — that’s the population frequency of your variant allele
  # Split by comma
  alleles <- strsplit(x, ",")[[1]]
  # convert to numeric, ignoring "."
  nums <- suppressWarnings(as.numeric(alleles))
  nums <- nums[!is.na(nums)]
  # usually reference AF is the first, so take the MIN non-ref
  if (length(nums) == 0) return(NA_real_)
  return(min(nums, na.rm = TRUE))
}

hotspot_df <- read_excel("Hotspots_compiled.xlsx", skip = 1) 

hotspot_df <- hotspot_df %>%
  mutate(
    aa_change = paste0(
      Reference_Amino_Acid,
      Amino_Acid_Position,
      Variant_Amino_Acid
    ),
    hotspot_pos = paste0(
      Reference_Amino_Acid,
      Amino_Acid_Position
    )
  )

# Remove synonymous variants automatically
hotspot_df <- hotspot_df %>%
  filter(Reference_Amino_Acid != Variant_Amino_Acid)

escape_regex <- function(x) {
  stringr::str_replace_all(
    x,
    "([.|()\\^{}+$*?]|\\[|\\]|\\\\)",
    "\\\\\\1"
  )
}

hotspots <- hotspot_df %>%
  mutate(
    aa_change = paste0(
      Reference_Amino_Acid,
      Amino_Acid_Position,
      Variant_Amino_Acid
    ),
    aa_change = escape_regex(aa_change)
  ) %>%
  group_by(gene = Hugo_Symbol) %>%
  summarise(
    pattern = paste0(
      "^p\\.(",
      stringr::str_c(sort(aa_change), collapse = "|"),
      ")"
    ),
    .groups = "drop"
  )

is_hotspot <- function(gene, prot_change) {
  if (is.na(gene) || is.na(prot_change)) return(FALSE)
  patt <- hotspots$pattern[match(gene, hotspots$gene)]
  if (is.na(patt)) return(FALSE)
  str_detect(prot_change, regex(patt, ignore_case = TRUE))
}

# Remove: rsID is no longer reliable as a germline indicator
classify_variants_mutect <- function(
    maf,
    oncogene_set,
    POP_AF_COMMON = 0.001,
    MIN_DP = 10,
    MIN_ALT_COUNT = 3,
    VAF_HET_LOW = 0.35,
    VAF_HET_HIGH = 0.65,
    VAF_HOM_HIGH = 0.90,
    TLOD_WEAK = 6,
    TLOD_MODERATE = 10
) {
  
  library(dplyr)
  library(purrr)
  
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  coding_classes <- c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "Splice_Site",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "DE_NOVO_START_IN_FRAME",
    "DE_NOVO_START_OUT_FRAME"
  )
  
  maf2 <- maf %>%
    mutate(
      
      tumor_vaf =
        t_alt_count / (t_alt_count + t_ref_count),
      
      CAF_AF =
        if ("dbSNP_CAF" %in% names(.))
          vapply(dbSNP_CAF, parse_popAF, numeric(1))
      else NA_real_,
      
      TOPMED_AF =
        if ("dbSNP_TOPMED" %in% names(.))
          vapply(dbSNP_TOPMED, parse_popAF, numeric(1))
      else NA_real_,
      
      popmax_af =
        pmax(
          ifelse(is.na(CAF_AF), 0, CAF_AF),
          ifelse(is.na(TOPMED_AF), 0, TOPMED_AF),
          na.rm = TRUE
        ),
      
      pop_common =
        popmax_af >= POP_AF_COMMON,
      
      pop_rare =
        !pop_common,
      
      hotspot =
        purrr::map2_lgl(
          Hugo_Symbol,
          Protein_Change,
          is_hotspot
        ),
      
      is_coding_or_splice =
        Variant_Classification %in% coding_classes,
      
      in_oncogene_list =
        Hugo_Symbol %in% oncogene_set,
      
      evidence_strength =
        case_when(
          is.na(TLOD) ~ NA_character_,
          TLOD < TLOD_WEAK ~ "weak",
          TLOD < TLOD_MODERATE ~ "moderate",
          TRUE ~ "strong"
        ),
      
      vaf_het_like =
        tumor_vaf >= VAF_HET_LOW &
        tumor_vaf <= VAF_HET_HIGH,
      
      vaf_homozygous =
        tumor_vaf >= VAF_HOM_HIGH,
      
      support_flag = case_when(
        t_alt_count < 3 ~ "low_alt",
        (t_alt_count + t_ref_count) < 10 ~ "low_depth",
        !is.na(TLOD) & TLOD < 10 ~ "weak_tlod",
        TRUE ~ "ok"
      )
      
    ) %>%
    
    mutate(
      
      priority = case_when(
        
        hotspot ~ "hotspot",
        
        !is.na(COSMIC_total_alterations_in_gene) &
          COSMIC_total_alterations_in_gene > 100 ~ "cosmic_recurrent",
        
        TRUE ~ "none"
      )
      
    ) %>%
    
    mutate(
      
      cancer_associated_gene = case_when(
        
        in_oncogene_list &
          is_coding_or_splice &
          pop_rare ~ "YES",
        
        TRUE ~ "NO"
      )
      
    ) %>%
    
    mutate(
      
      classification = case_when(
        
        hotspot ~ "somatic",
        
        priority == "cosmic_recurrent" &
          !vaf_het_like & !vaf_homozygous ~ "somatic",
        
        pop_common & vaf_het_like ~ "germline",
        
        pop_common & vaf_homozygous ~ "germline",
        
        pop_rare & vaf_homozygous ~ "review_cn_loh",
        
        is_coding_or_splice & pop_rare &
          !vaf_het_like & !vaf_homozygous ~ "functional_uncertain_somatic",

        pop_rare & !vaf_het_like & !vaf_homozygous ~ "likely_somatic",
        
        is_coding_or_splice ~ "functional_uncertain",
        
        TRUE ~ "uncertain"
      )
      
    )
  
  maf_somatic <- maf2 %>%
    filter(
      classification %in% c(
        "somatic",
        "functional_uncertain_somatic",
        "likely_somatic"
      )
    )
  
  return(list(
    somatic = maf_somatic,
    all = maf2
  ))
  
}

# # Usage
# sample = "VCU-BC-052"
# classified_results <- classify_variants_mutect(
#   maf = maf_list[[sample]],
#   oncogene_set = oncogene_set
# )

# classified_results$`VCU-BC-052`$somatic %>%
#   filter(
#     classification %in% c(
#       "somatic",
#       "functional_uncertain_somatic",
#       "likely_somatic"
#     ),
#     Variant_Classification %in% c(
#       "Missense_Mutation",
#       "Nonsense_Mutation",
#       "Nonstop_Mutation",
#       "Translation_Start_Site",
#       "Splice_Site",
#       "Frame_Shift_Del",
#       "Frame_Shift_Ins",
#       "In_Frame_Del",
#       "In_Frame_Ins",
#       "DE_NOVO_START_IN_FRAME",
#       "DE_NOVO_START_OUT_FRAME"
#     )
#   ) %>%
#   count(classification)

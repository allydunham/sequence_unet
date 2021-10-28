#!/usr/bin/env Rscript
# Analyse proteomic measurements against Sequence UNET predictions
source("src/config.R")

import_abundance <- function() {
  df <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism"), skip = 1) %>%
    select(-row)
  
  protein_groups <- str_split(abundance$proteins, ";")
  group_counts <- map_int(proteins, length)

  tibble(protein = unlist(protein_groups),
         intensity = rep(abundance$intensity, times = group_counts),
         organism = rep(abundance$organism, times = group_counts))
}

abundance <- import_abundance()

preds <- read_tsv("data/abundance/muller_classifier_summary.tsv") %>%
  mutate(gene = str_match(gene, "[a-z]*\\|([A-Z0-9]*)\\|[A-Z0-9]*")[,2])

omics <- left_join(abundance, preds, by = c("protein"="gene")) %>%
  drop_na()

cor.test(omics$intensity, omics$mean_mut)

#!/usr/bin/env Rscript
# Initial proteome analysis to identify features to extract from main data
source("src/config.R")

# Take a random subset of proteins
# preds <- read_tsv("data/abundance/split/muller_classifier_0.tsv") %>%
#   mutate(gene = str_match(gene, "[a-z]*\\|([A-Z0-9]*)\\|[A-Z0-9]*")[,2])
# subset <- sample(unique(preds$gene), size = 1000, replace = FALSE)
# sub_preds <- filter(preds, gene %in% subset)
# write_tsv(sub_preds, "data/abundance/split/muller_classifier_subset.tsv")

import_abundance <- function() {
  df <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism"), skip = 1) %>%
    select(-row)
  
  protein_groups <- str_split(df$proteins, ";")
  group_counts <- map_int(protein_groups, length)
  
  tibble(protein = unlist(protein_groups),
         intensity = rep(df$intensity, times = group_counts),
         organism = rep(df$organism, times = group_counts))
}

abundance <- import_abundance() %>%
  group_by(organism) %>%
  mutate(intensity_fc = log2((intensity + 1) / median(intensity))) %>%
  ungroup()

preds <- inner_join(abundance, read_tsv("data/abundance/split/muller_classifier_subset.tsv"), by = c("protein"="gene"))

#### Analysis ####
# Average conservation

conserved_threshold <- 0.5
protein_summary <- pivot_longer(preds, A:Y, names_to = "mut", values_to = "pred") %>%
  group_by(protein, intensity, intensity_fc, organism, position, wt) %>%
  summarise(mean_mut_pred = mean(pred[wt != mut])) %>%
  group_by(protein, intensity, intensity_fc, organism) %>%
  summarise(mean_pos = mean(mean_mut_pred),
            percent_conserved = sum(mean_mut_pred > conserved_threshold) / n())



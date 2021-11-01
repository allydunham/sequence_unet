#!/usr/bin/env Rscript
# Prepare data for proteomics analysis
source("src/config.R")

### Preds ###
# Del thresholds used:
# classifier - 0.44
# pssm - 0.01
# clinvar - 0.93
# For clinvar and classifier maximise sqrt((1 - fpr)^2 + tpr^2)) to get closest to top left
# For pssm use the same cutoff as classifier training

### Abundance Intensity ###
abundance <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism"), skip = 1) %>%
  select(-row)

fasta <- Biostrings::readAAStringSet("data/abundance/muller.fa")
protein_lengths <- tibble(protein = str_match(names(fasta), "[a-z]{2}\\|([A-Z0-9]*)\\|.*")[,2],
                          length = Biostrings::width(fasta))

protein_groups <- str_split(abundance$proteins, ";")
group_counts <- map_int(protein_groups, length)
processed_abundance <- tibble(organism = rep(abundance$organism, times = group_counts),
                              protein = unlist(protein_groups),
                              intensity = rep(abundance$intensity, times = group_counts)) %>%
  left_join(protein_lengths, by = "protein") %>%
  group_by(organism) %>%
  mutate(intensity_per_len = intensity / length,
         intensity_fc = log2((intensity + 1) / median(intensity)),
         intensity_fc_per_len = log2((intensity_per_len + 1) / median(intensity_per_len, na.rm = TRUE))) %>%
  ungroup()

### SIFT ###


### Combine ###


#!/usr/bin/env Rscript
# Generate correlations between tools and DMS data
source('src/config.R')
source("src/analysis.R")

# Import DMS data
dms <- read_tsv("data/dms/long_combined_mutational_scans.tsv") %>%
  left_join(select(read_tsv("data/dms/gene_summary.tsv"), gene = Gene, uniprot = `Uniprot ID`), by = "gene") %>%
  filter(class == "Missense") %>%
  select(gene, uniprot, position, wt, mut, score, SIFT4G=sift, FoldX = total_energy)

# Import Sequence UNET data
unet <- bind_rows(
  `Baseline ClinVar` = read_tsv("data/dms/preds/baseline_clinvar.tsv"),
  `Baseline Frequency` = read_tsv("data/dms/preds/baseline_freq.tsv"),
  `UNET (Top)` = read_tsv("data/dms/preds/clinvar_top.tsv"),
  `UNET (Finetune)` = read_tsv("data/dms/preds/clinvar_finetune.tsv"),
  `UNET` = read_tsv("data/dms/preds/unet_freq.tsv"),
  .id = "model"
) %>% pivot_wider(names_from = model, values_from = pred)

# Import EVE data
uniprot_mapping <- read_tsv("data/dms/uniprot_mapping.tsv") %>% # Only Uniprot IDs in DMS included
  select(uniprot = Entry, uniprot_name = `Entry name`)

load_eve <- function(uniprot_name) {
  path <- str_c("data/eve_variants/", uniprot_name, ".csv")
  
  if (!file.exists(path)) {
    warning(str_c("Gene not found: ", uniprot_name))
    return()
  }
  
  eve_cols <- cols(.default = col_character(), position = col_integer(), EVE_scores_ASM = col_double())
  
  read_csv(path, col_types = eve_cols) %>%
    mutate(uniprot_name = uniprot_name) %>%
    left_join(uniprot_mapping, by = "uniprot_name") %>%
    select(uniprot, uniprot_name, position, wt = wt_aa, mut = mt_aa, `EVE` = EVE_scores_ASM) %>%
    return()
}

eve_scores <- unique(uniprot_mapping$uniprot_name) %>%
  map(~suppressMessages(load_eve(.))) %>%
  bind_rows() %>%
  drop_na()

# Generate Query for dbNPFS (for applicable human proteins)


# Import dbNPFS data


# Combine data 


# Calculate correlation

#! summarise number of studies with data too (often only a subset of human)


# Move ROC calculation here and make Jelier own script?
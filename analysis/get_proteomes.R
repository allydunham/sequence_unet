#!/usr/bin/env Rscript
# Identify proteomes to fetch from Uniprot
source("src/config.R")

uniprot <- read_tsv("data/abundance/uniprot_proteomes", col_names = c("uniprot_proteome", "name", "tax_id", "proteins"))

muller <- readxl::read_xlsx("data/abundance/muller_species.xlsx", range = "A2:C131", col_names = c("kingdom", "name", "tax_id"))

cross <- filter(uniprot, tax_id %in% muller$tax_id) %>% 
  group_by(tax_id) %>%
  filter(proteins == max(proteins)) %>%
  ungroup() %>%
  distinct(tax_id, .keep_all = TRUE)

proteomics <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism")) %>%
  select(-row)

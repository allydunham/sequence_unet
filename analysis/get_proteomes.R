#!/usr/bin/env Rscript
# Identify proteomes to fetch from Uniprot
source("src/config.R")

proteomics <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism"), skip = 1) %>%
  select(-row)

uniprot <- str_split(proteomics$proteins, ";") %>%
  flatten_chr() %>%
  unique() %>%
  sort()

cat(uniprot, file = "data/abundance/muller_uniprot_ids", sep = "\n")

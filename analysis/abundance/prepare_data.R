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

classifier <- read_tsv("data/abundance/muller_classifier_summary.tsv") %>%
  extract(gene, into = "gene", regex = "[a-z]{2}\\|([0-9A-Z]*)\\|.*")

pssm <- read_tsv("data/abundance/muller_pssm_summary.tsv") %>%
  extract(gene, into = "gene", regex = "[a-z]{2}\\|([0-9A-Z]*)\\|.*")

top_model <- read_tsv("data/abundance/muller_top_summary.tsv") %>%
  extract(gene, into = "gene", regex = "[a-z]{2}\\|([0-9A-Z]*)\\|.*")

### Abundance Intensity ###
abundance <- read_csv("data/abundance/muller_proteomics.csv", col_names = c("row", "proteins", "intensity", "organism"), skip = 1) %>%
  select(-row) %>%
  mutate(organism = ifelse(organism == "Dentriovibrio acetiphilus", "Denitrovibrio acetiphilus", organism))
  
fasta <- Biostrings::readAAStringSet("data/abundance/muller.fa")
protein_lengths <- tibble(protein = names(fasta), length = Biostrings::width(fasta)) %>%
  extract(protein, into = c("source", "protein"), regex = "([a-z]{2})\\|([0-9A-Z]*)\\|.*") %>%
  mutate(source = c(tr="TrEMBL", sp="SwissProt")[source])

species <- readxl::read_xlsx("data/abundance/muller_species.xlsx", sheet = 2)

protein_groups <- str_split(abundance$proteins, ";")
group_counts <- map_int(protein_groups, length)
processed_abundance <- tibble(organism = rep(abundance$organism, times = group_counts),
                              protein = unlist(protein_groups),
                              intensity = rep(abundance$intensity, times = group_counts)) %>%
  left_join(protein_lengths, by = "protein") %>%
  left_join(species, by = c("organism" = "species")) %>%
  group_by(organism) %>%
  mutate(intensity_per_len = intensity / length,
         intensity_fc = log2((intensity + 1) / median(intensity)),
         intensity_fc_per_len = log2((intensity_per_len + 1) / median(intensity_per_len, na.rm = TRUE))) %>%
  ungroup() %>%
  select(superkingdom, organism, taxid, source, protein, length, everything())

### SIFT ###
# Set E. coli SIFT4G uniprot IDs to ED1a equivalents strain rather than K12 since this is the strain with abundances
ed1a_k12 <- read_tsv("data/abundance/ed1a_to_k12.tsv")

sift <- read_tsv("data/abundance/mutfunc_sift_summary.tsv") %>%
  rename(protein = acc, mean_mut = mean_sift) %>%
  select(-organism) %>%
  left_join(select(ed1a_k12, ed1a, k12), by = c("protein"="k12")) %>%
  mutate(protein = ifelse(is.na(ed1a), protein, ed1a)) %>%
  select(-ed1a)

### Combine ###
comb <- bind_rows(
  SIFT4G = left_join(processed_abundance, sift, by = "protein") %>% drop_na(mean_mut),
  `UNET Classifier` = left_join(processed_abundance, classifier, by = c("protein" = "gene")) %>% drop_na(mean_mut),
  `UNET PSSM` = left_join(processed_abundance, pssm, by = c("protein" = "gene")) %>% drop_na(mean_mut),
  `UNET Top` = left_join(processed_abundance, top_model, by = c("protein" = "gene")) %>% drop_na(mean_mut),
  .id = "tool"
) %>%
  select(superkingdom, organism, taxid, source, protein, length, tool, everything()) %>%
  arrange(superkingdom, organism, protein, tool)
write_tsv(comb, "data/abundance/muller_proteomics_processed.tsv")

#!/usr/bin/env Rscript
# Initial proteome analysis to identify features to extract from main data
source("src/config.R")

# Take a random subset of proteins
# preds <- read_tsv("data/abundance/split/muller_classifier_0.tsv") %>%
#   mutate(gene = str_match(gene, "[a-z]*\\|([A-Z0-9]*)\\|[A-Z0-9]*")[,2])
# subset <- sample(unique(preds$gene), size = 1000, replace = FALSE)
# sub_preds <- filter(preds, gene %in% subset)
# write_tsv(sub_preds, "data/abundance/split/muller_classifier_subset.tsv")

preds <- inner_join(read_tsv("data/abundance/muller_proteomics_processed.tsv"),
                    read_tsv("data/abundance/split/muller_classifier_subset.tsv"), by = c("protein"="gene"))

#### Analysis ####
# Average conservation
conserved_threshold <- 0.5
pos_summary <- pivot_longer(preds, A:Y, names_to = "mut", values_to = "pred") %>%
  group_by(protein, intensity, intensity_per_len, intensity_fc, intensity_fc_per_len, organism, position, wt) %>%
  summarise(mean_mut_pred = mean(pred[wt != mut]),
            n_conserved = sum(pred > conserved_threshold), .groups = "drop")

protein_summary <- group_by(pos_summary, protein, intensity, intensity_per_len, intensity_fc, intensity_fc_per_len, organism) %>%
  summarise(mean_pos = mean(mean_mut_pred),
            percent_conserved = sum(mean_mut_pred > conserved_threshold) / n(),
            mean_conserved_muts = mean(n_conserved),
            percent_n_conserved = sum(n_conserved > 9) / n(),
            .groups = "drop")

cor_summary <- select(protein_summary, -intensity, -intensity_per_len, -intensity_fc) %>%
  pivot_longer(mean_pos:percent_n_conserved, names_to = "metric", values_to = "pred") %>%
  group_by(organism, metric) %>%
  group_modify(~tidy(cor.test(.$intensity_fc_per_len, .$pred))) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

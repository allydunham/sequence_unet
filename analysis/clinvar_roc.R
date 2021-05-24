#!/usr/bin/env Rscript
# Calculate ROC curves for ClinVar variants from the PSSM top models
source('src/analysis/config.R')
source("src/analysis/analysis.R")
data("BLOSUM62", package = "Biostrings")
plots <- list()

blosum <- as_tibble(BLOSUM62, rownames = 'wt') %>%
  pivot_longer(-wt, names_to = 'mut', values_to = 'BLOSUM62')

clinvar_stats <- read_tsv("data/clinvar/clinvar_test.tsv") %>% 
  mutate(pdb_id = str_to_upper(pdb_id)) %>%
  select(uniprot, position, wt, mut, pdb_id, chain, pdb_pos, clnsig, clnsig_patho, foldx_ddg, sift_score) %>%
  left_join(blosum, by = c("wt", "mut"))

preds <- bind_rows(
  # Trained on ProteinNet
  `Baseline Frequency` = read_tsv("data/clinvar/preds/baseline_freq.tsv"),
  `UNET DMS` = read_tsv("data/clinvar/preds/unet_dms.tsv"),
  `UNET Frequency (End to End)` = read_tsv("data/clinvar/preds/unet_end_to_end_freq.tsv"),
  `UNET Frequency (Top Model)` = read_tsv("data/clinvar/preds/unet_top_model_freq.tsv"),
  
  # Trained on ClinVar (incl. some ProteinNet pre-training for some models)
  `Baseline ClinVar` = read_tsv("data/clinvar/preds/baseline_clinvar.tsv"),
  `True PSSM` = read_tsv("data/clinvar/preds/true_pssm.tsv"),
  `Pred PSSM` = read_tsv("data/clinvar/preds/unet_pred_pssm.tsv"),
  `UNET Features (Sequence)` = read_tsv("data/clinvar/preds/unet_pred_features_seq.tsv"),
  `UNET Features (Structure)` = read_tsv("data/clinvar/preds/unet_pred_features_struct.tsv"),
  `UNET Features (Frequency)` = read_tsv("data/clinvar/preds/unet_pred_features_classifier.tsv"),
  .id = "model"
) %>%
  select(-wt) %>%
  pivot_wider(names_from = model, values_from = pred) %>%
  left_join(clinvar_stats, ., by = c("pdb_id", "chain", "pdb_pos", "mut")) %>%
  select(-pdb_id, -chain, -pdb_pos) %>%
  rename(SIFT4G = sift_score, FoldX = foldx_ddg)

# Scores where less than indicates deleterious
less = c("SIFT4G", "BLOSUM62")
roc <- pivot_longer(preds, c(-uniprot, -position, -wt, -mut, -clnsig, -clnsig_patho), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, clnsig_patho, pred, greater = !(.y %in% less))) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))

unet_models <- c("UNET Features (Frequency)", "UNET Features (Structure)", "UNET Frequency (End to End)", 
                 "UNET Features (Sequence)", "UNET DMS", "UNET Frequency (Top Model)")
plots$features_roc <- filter(roc, model %in% unet_models) %>%
  ggplot(aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

plots$features_pr <-  filter(roc, model %in% unet_models) %>%
  ggplot(aes(x = tpr, y = precision, colour = model)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  geom_line() +
  labs(x = 'Recall', y = 'Precision') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

comparisons <- c("UNET Features (Frequency)", "SIFT4G", "Baseline ClinVar", 
                 "Baseline Frequency", "FoldX", "BLOSUM62", "True PSSM", "Pred PSSM")
plots$comparisons_roc <- filter(roc, model %in% comparisons) %>%
  ggplot(aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Set1', name = '')

plots$comparisons_pr <- filter(roc, model %in% comparisons) %>%
  ggplot(aes(x = tpr, y = precision, colour = model)) + 
    geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
    geom_line() +
    labs(x = 'Recall', y = 'Precision') +
    scale_colour_brewer(type = 'qual', palette = 'Set1', name = '')
  
save_plotlist(plots, "figures/clinvar_unet/", overwrite = "all")

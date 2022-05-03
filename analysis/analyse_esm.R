#!/usr/bin/env Rscript
# Analyse ESM models
source('src/config.R')
source("src/analysis.R")

import_preds <- function() {
  true_freqs <- read_tsv('data/pssm/pn_casp12_testing.tsv') %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "freq")
  
  base_models <- read_tsv("data/esm/testing_logit_preds.tsv") %>%
    rename(protein = id) %>%
    pivot_longer(pred_A:pred_Y, names_to = "mut", names_prefix = "pred_", values_to = "pred") %>%
    pivot_wider(names_from = model, values_from = pred) %>%
    mutate(esm1v = rowMeans(across(starts_with("esm1v")))) %>%
    select(protein, position, wt, mut, esm1b, esm1v) %>%
    pivot_longer(c(esm1b, esm1v), names_to = "model", values_to = "pred") %>%
    mutate(pred = exp(pred))
  
  top_models <- bind_rows(
    freq_10 = read_tsv("data/esm/freq_preds_10_epoch.tsv"),
    freq_20 = read_tsv("data/esm/freq_preds_20_epoch.tsv"),
    pssm_10 = read_tsv("data/esm/pssm_preds_10_epoch.tsv"),
    pssm_20 = read_tsv("data/esm/pssm_preds_20_epoch.tsv"),
    .id = "model"
  ) %>%
    select(model, protein = id, position, wt, starts_with("pred")) %>%
    pivot_longer(pred_A:pred_Y, names_to = "mut", names_prefix = "pred_", values_to = "pred")
  
  bind_rows(full_join(true_freqs, base_models, by = c("protein", "position", "wt", "mut")),
            full_join(true_freqs, top_models, by = c("protein", "position", "wt", "mut")))
}

preds <- import_preds()

# PSSM correlation

pssm_cor <- filter(preds, model %in% c("esm1b", "esm1v", "pssm_10", "pssm_20")) %>%
  group_by(model) %>%
  group_modify(~broom::tidy(cor.test(.$pred, .$freq)))

freq_roc <- filter(preds, model %in% c("esm1b", "esm1v", "freq_10", "freq_20")) %>%
  mutate(deleterious = freq < 0.01) %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, deleterious, pred, greater = .y %in% c("freq_10", "freq_20"), max_steps = 6000)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))

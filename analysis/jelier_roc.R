#!/usr/bin/env Rscript
# Calculate ROC curves for Jelier and DMS predictions
source('src/config.R')
source("src/analysis.R")
plots <- list()

import_jelier <- function() {
  base <- read_tsv("data/jelier/jelier_variants.tsv") %>%
    mutate(class = effect == "deleterious") %>%
    left_join(read_tsv("data/jelier/preds/sift_foldx.tsv"), by = c("gene", "position", "wt", "mut")) %>%
    rename(SIFT4G = sift_score, FoldX = foldx)
  
  esm <- read_tsv("data/esm/1v_jelier_preds.tsv") %>%
    rename(gene = id) %>%
    semi_join(base, by = c("gene", "position", "wt")) %>%
    pivot_longer(starts_with("pred_"), names_to = "mut", names_prefix = "pred_", values_to = "s") %>%
    pivot_wider(names_from = "model", values_from = "s") %>%
    mutate(`ESM-1v` = rowMeans(across(c(`1`, `2`, `3`, `4`, `5`)))) %>%
    select(gene, position, wt, mut, `ESM-1v`)
    
  bind_rows(
    `Baseline ClinVar` = read_tsv("data/jelier/preds/baseline_clinvar.tsv"),
    `Baseline Frequency` = read_tsv("data/jelier/preds/baseline_freq.tsv"),
    `UNET (Top)` = read_tsv("data/jelier/preds/clinvar_top.tsv"),
    `UNET (Finetune)` = read_tsv("data/jelier/preds/clinvar_finetune.tsv"),
    `UNET` = read_tsv("data/jelier/preds/unet_freq.tsv"),
    .id = "model"
  ) %>%
    distinct(model, gene, position, wt, mut, .keep_all = TRUE) %>%
    pivot_wider(names_from = model, values_from = pred) %>% 
    left_join(base, ., by = c("gene", "position", "wt", "mut")) %>%
    left_join(esm, by = c("gene", "position", "wt", "mut"))
}

jelier_preds <- import_jelier()

jelier_roc <- pivot_longer(jelier_preds, c(-gene, -position, -wt, -mut, -effect, -class), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, class, pred, greater = !(.y %in% c("SIFT4G", "ESM-1v")), max_steps = 6000)) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))
write_tsv(jelier_roc, "data/jelier/roc.tsv")

plots$jelier_roc <- ggplot(jelier_roc, aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

plots$jelier_pr <- ggplot(jelier_roc, aes(x = tpr, y = precision, colour = model)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  geom_line() +
  labs(x = 'Recall', y = 'Precision') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

save_plotlist(plots, "figures/generalisation", overwrite = "all")

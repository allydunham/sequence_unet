#!/usr/bin/env Rscript
# Calculate ROC curves for Jelier and DMS predictions
source('src/config.R')
source("src/analysis.R")
plots <- list()

import_dms <- function() {
  long <- read_tsv("data/dms/long_combined_mutational_scans.tsv") %>%
    select(gene, position, wt, mut, SIFT4G=sift, FoldX = total_energy) %>%
    mutate(gene = str_to_lower(str_remove_all(gene, " ")))
  preds <- bind_rows(
    `Baseline Clinvar` = read_tsv("data/dms/preds/baseline_clinvar.tsv"),
    `Baseline Freq` = read_tsv("data/dms/preds/baseline_freq.tsv"),
    `UNET Top` = read_tsv("data/dms/preds/clinvar_top.tsv"),
    `UNET` = read_tsv("data/dms/preds/unet_freq.tsv"),
    .id = "model"
  ) %>% pivot_wider(names_from = model, values_from = pred)
  
  read_tsv("data/dms/test.tsv") %>%
    select(-study) %>%
    arrange(gene, position) %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "score") %>%
    left_join(long, by = c("gene", "position", "wt", "mut")) %>%
    left_join(preds, by = c("gene", "position", "wt", "mut")) %>%
    mutate(class = score < -0.5)
}

import_jelier <- function() {
  base <- read_tsv("data/jelier/jelier_variants.tsv") %>%
    mutate(class = effect == "deleterious") %>%
    left_join(read_tsv("data/jelier/preds/sift_foldx.tsv"), by = c("gene", "position", "wt", "mut")) %>%
    rename(SIFT4G = sift_score, FoldX = foldx)
  
  bind_rows(
    `Baseline Clinvar` = read_tsv("data/jelier/preds/baseline_clinvar.tsv"),
    `Baseline Freq` = read_tsv("data/jelier/preds/baseline_freq.tsv"),
    `UNET Top` = read_tsv("data/jelier/preds/clinvar_top.tsv"),
    `UNET` = read_tsv("data/jelier/preds/unet_freq.tsv"),
    .id = "model"
  ) %>%
    distinct(model, gene, position, wt, mut, .keep_all = TRUE) %>%
    pivot_wider(names_from = model, values_from = pred) %>% 
    left_join(base, ., by = c("gene", "position", "wt", "mut"))
}

dms_preds <- import_dms()
jelier_preds <- import_jelier()

jelier_roc <- pivot_longer(jelier_preds, c(-gene, -position, -wt, -mut, -effect, -class), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, class, pred, greater = !(.y == "SIFT4G"), max_steps = 6000)) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))
write_tsv(jelier_roc, "data/jelier/roc.tsv")

dms_roc <- pivot_longer(dms_preds, c(-gene, -position, -wt, -mut, -score, -class), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, class, pred, greater = !(.y == "SIFT4G"), max_steps = 6000)) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))
write_tsv(dms_roc, "data/dms/roc.tsv")

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

plots$dms_roc <- ggplot(dms_roc, aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

plots$dms_pr <- ggplot(dms_roc, aes(x = tpr, y = precision, colour = model)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  geom_line() +
  labs(x = 'Recall', y = 'Precision') +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', name = '')

save_plotlist(plots, "figures/generalisation", overwrite = "all")

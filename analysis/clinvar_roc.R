#!/usr/bin/env Rscript
# Calculate ROC curves for ClinVar variants models
source('src/config.R')
source("src/analysis.R")
plots <- list()

clinvar_stats <- read_tsv("data/clinvar/clinvar_test.tsv") %>% 
  mutate(pdb_id = str_to_upper(pdb_id)) %>%
  select(uniprot, position, wt, mut, pdb_id, chain, pdb_pos, clnsig, clnsig_patho)

fa_preds <- bind_rows(
  `UNET` = read_tsv("data/clinvar/preds/unet_freq_fa.tsv"),
  `UNET (Top)` = read_tsv("data/clinvar/preds/unet_freq_features_top_fa.tsv"),
  `UNET (Finetune)` = read_tsv("data/clinvar/preds/unet_freq_finetune_fa.tsv"),
  .id = "model"
) %>%
  pivot_wider(names_from = model, values_from = pred) %>%
  extract(gene, "uniprot", "[a-z]*\\|([A-Z0-9]*)\\|.*") %>%
  left_join(clinvar_stats, ., by = c("uniprot", "position", "wt", "mut")) %>%
  select(-pdb_id, -chain, -pdb_pos)

pdb_preds <- bind_rows(
  # Trained on ProteinNet
  `Baseline Frequency*` = read_tsv("data/clinvar/preds/baseline_freq.tsv"),
  `UNET*` = read_tsv("data/clinvar/preds/unet_freq.tsv"),
  `PreGraph UNET*` = read_tsv("data/clinvar/preds/unet_freq_structure.tsv"),
  
  # Freq Thresholds
  `UNET Thresh 0.1` = read_tsv("data/clinvar/preds/unet_freq_0.1.tsv"),
  `UNET Thresh 0.01` = read_tsv("data/clinvar/preds/unet_freq_0.01.tsv"),
  `UNET Thresh 0.001` = read_tsv("data/clinvar/preds/unet_freq_0.001.tsv"),
  `UNET Thresh 0.0001` = read_tsv("data/clinvar/preds/unet_freq_0.0001.tsv"),
  
  # Trained on ClinVar (incl. some ProteinNet pre-training for some models)
  `Baseline ClinVar*` = read_tsv("data/clinvar/preds/baseline_clinvar.tsv"),
  `UNET (Top)*` = read_tsv("data/clinvar/preds/unet_freq_features_top.tsv"),
  `PreGraph UNET (Top)*` = read_tsv("data/clinvar/preds/unet_freq_structure_features_top.tsv"),
  `UNET (Finetune)*` = read_tsv("data/clinvar/preds/unet_freq_finetune.tsv"),
  `PreGraph UNET (Finetune)*` = read_tsv("data/clinvar/preds/unet_freq_structure_finetune.tsv"),
  .id = "model"
) %>%
  select(-wt) %>%
  pivot_wider(names_from = model, values_from = pred) %>%
  left_join(clinvar_stats, ., by = c("pdb_id", "chain", "pdb_pos", "mut")) %>%
  select(-pdb_id, -chain, -pdb_pos)

preds <- full_join(fa_preds, pdb_preds, by = c("uniprot", "position", "wt", "mut", "clnsig", "clnsig_patho"))

# Scores where less than indicates deleterious
roc <- pivot_longer(preds, c(-uniprot, -position, -wt, -mut, -clnsig, -clnsig_patho), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, clnsig_patho, pred, greater = TRUE, max_steps = 6000)) %>%
  mutate(pr_auc = pr_auc(tpr, precision)) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))
write_tsv(roc, "data/clinvar/roc.tsv")

plots$roc <- ggplot(filter(roc, !str_detect(model, "Thresh")), aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Set3', name = '')

plots$pr <- ggplot(filter(roc, !str_detect(model, "Thresh")), aes(x = tpr, y = precision, colour = model)) + 
    geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
    geom_line() +
    labs(x = 'Recall', y = 'Precision') +
    scale_colour_brewer(type = 'qual', palette = 'Set3', name = '')

plots$roc_thresh <- ggplot(filter(roc, str_detect(model, "Thresh")), aes(x = fpr, y = tpr, colour = model_auc)) + 
  geom_abline(slope = 1, linetype = 'dashed', colour = 'black') +
  geom_step(direction = "hv") +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  scale_colour_brewer(type = 'qual', palette = 'Set3', name = '')

save_plotlist(plots, "figures/clinvar_unet/", overwrite = "all")


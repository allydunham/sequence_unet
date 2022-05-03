#!/usr/bin/env Rscript
# Analyse results of PSSM Profile predictions
source('src/config.R')
source("src/analysis.R")
data("BLOSUM62", package = "Biostrings")
plots <- list()

### Prepare Data ###
blosum <- as_tibble(BLOSUM62, rownames = 'wt') %>%
  pivot_longer(-wt, names_to = 'mut', values_to = 'BLOSUM62')

read_sift <- function(x) {
  protein <- str_split(basename(x), "\\.", simplify = TRUE)[1,1]
  col_names <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                 "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y", "Z", "*", "-")
  suppressMessages(read_table(x, skip = 5, comment = "//", col_names = col_names) %>%
    mutate(protein = protein, position = seq_along(protein)))
}

sift <- map(dir("data/freq/sift_test", full.names = TRUE), read_sift) %>%
  bind_rows() %>%
  pivot_longer(c(-protein, -position), names_to = "mut", values_to = "SIFT4G") %>%
  filter(mut %in% Biostrings::AA_STANDARD)

esm1b_top <- read_tsv('data/esm/freq_preds_10_epoch.tsv') %>%
  select(protein=id, position, wt, starts_with("pred")) %>%
  pivot_longer(starts_with("pred"), names_to = "mut", names_prefix = "pred_", values_to = "ESM-1b Top Model")

esm1b_logits <- read_tsv("data/esm/testing_logit_preds.tsv") %>%
  filter(model == "esm1b") %>%
  select(protein=id, position, wt, starts_with("pred")) %>%
  pivot_longer(pred_A:pred_Y, names_to = "mut", names_prefix = "pred_", values_to = "ESM-1b Logits")

models <- read_tsv('data/pssm/pn_casp12_testing.tsv') %>%
  select(protein, position, wt, A:Y) %>%
  pivot_longer(A:Y, names_to = "mut", values_to = "freq") %>%
  mutate(deleterious = freq < 0.01) %>%
  left_join(select(read_tsv('data/freq/unet_sequence_testing.tsv'), protein = pdb_id, position, wt, mut, UNET=pred),
            by = c("protein", "position", "wt", "mut")) %>%
  left_join(select(read_tsv('data/freq/unet_structure_testing.tsv'), protein = pdb_id, position, wt, mut, `PreGraph UNET`=pred),
            by = c("protein", "position", "wt", "mut")) %>%
  left_join(select(read_tsv('data/freq/baseline_testing.tsv'), protein = pdb_id, position, wt, mut, `Baseline CNN`=pred),
            by = c("protein", "position", "wt", "mut")) %>%
  left_join(esm1b_top, by = c("protein", "position", "wt", "mut")) %>%
  left_join(esm1b_logits, by = c("protein", "position", "wt", "mut")) %>%
  left_join(blosum, by = c("wt", "mut")) %>%
  left_join(sift, by = c("protein", "position", "mut")) %>%
  pivot_longer(c(UNET:SIFT4G), names_to = "model", values_to = "pred")
write_tsv(models, "data/freq/all_model_testing.tsv")

### Analyse ###
greater <- c(UNET = TRUE, `PreGraph UNET` = TRUE, SIFT4G = FALSE, `ESM-1b Logits`=FALSE, `ESM-1b Top Model`=TRUE, BLOSUM62 = FALSE, `Baseline CNN` = TRUE)
roc <- group_by(models, model) %>%
  group_modify(~calc_roc(., deleterious, pred, greater = greater[.y$model], max_steps = 6000)) %>%
  mutate(model_auc = auc_labeled_model(model, auc))
write_tsv(roc, "data/freq/roc.tsv")

plots$roc <- ggplot(roc, aes(x = fpr, y = tpr, colour = model_auc)) +
  geom_step() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_brewer(name = "", type = "qual", palette = "Dark2") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

plots$pr <- ggplot(roc, aes(x = tpr, y = precision, colour = model)) +
  geom_step() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_colour_brewer(name = "", type = "qual", palette = "Dark2") +
  labs(x = "Recall", y = "Precision")

### Save Plots ###
save_plotlist(plots, 'figures/freq_predictions', overwrite='all')


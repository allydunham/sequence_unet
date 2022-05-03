#!/usr/bin/env Rscript
# Analyse results of PSSM Profile predictions
source('src/config.R')
data("BLOSUM62", package = "Biostrings")

### Prepare Data ###
blosum <- as_tibble(BLOSUM62, rownames = 'wt') %>%
  pivot_longer(-wt, names_to = 'mut', values_to = 'BLOSUM62')

# from ExPASy Data Portal (https://web.expasy.org/docs/relnotes/relstat.html, 17/08/2020)
aa_freqs <- c('A'=8.25, 'Q'=3.93, 'L'=9.65, 'S'=6.63, 'R'=5.53, 'E'=6.72, 'K'=5.80, 'T'=5.35,
              'N'=4.05, 'G'=7.08, 'M'=2.41, 'W'=1.09, 'D'=5.46, 'H'=2.27, 'F'=3.86, 'Y'=2.92,
              'C'=1.38, 'I'=5.91, 'P'=4.73, 'V'=6.86) / 100

model_files <- c(ProteinNet='data/pssm/pn_casp12_testing.tsv',
                 SPBuild='data/pssm/spbuild_casp12_testing.tsv',
                 UNET='data/pssm/unet_sequence_testing.tsv',
                 `PreGraph UNET`='data/pssm/unet_structure_testing.tsv',
                 `Baseline CNN`='data/pssm/baseline_testing.tsv',
                 `ESM-1b Logits`='data/esm/testing_logit_preds.tsv',
                 `ESM-1b Top Model`='data/esm/pssm_preds_10_epoch.tsv')

# Convert output into the standard  log(O) / Log(E) format
models <- map(model_files, read_tsv)

models$ProteinNet <- pivot_longer(models$ProteinNet, A:Y, names_to = "mut", values_to = "pred") %>%
  mutate(pred = as.integer(log2((pred + 0.00001) / aa_freqs[mut]))) %>%
  drop_na()

models$SPBuild <- pivot_longer(models$SPBuild, A:V, names_to = "mut", values_to = "pred") %>%
  drop_na()

models$`ESM-1b Logits` <- filter(models$`ESM-1b Logits`, model == "esm1b") %>%
  select(protein=id, position, wt, starts_with("pred")) %>%
  pivot_longer(starts_with("pred"), names_prefix = "pred_", names_to = "mut", values_to = "pred") %>%
  mutate(pred = as.integer(log2((exp(pred) + 0.00001) / aa_freqs[mut]))) %>%
  drop_na()

models$`ESM-1b Top Model` <- select(models$`ESM-1b Top Model`, protein=id, position, wt, starts_with("pred")) %>%
  pivot_longer(starts_with("pred"), names_prefix = "pred_", names_to = "mut", values_to = "pred") %>%
  mutate(pred = as.integer(log2((exp(pred) + 0.00001) / aa_freqs[mut]))) %>%
  drop_na()

models$UNET <- mutate(models$UNET, pred = as.integer(log2((pred + 0.00001) / aa_freqs[mut]))) %>% 
  select(-chain) %>% 
  rename(protein = pdb_id) %>% 
  drop_na()

models$`PreGraph UNET` <- mutate(models$`PreGraph UNET`, pred = as.integer(log2((pred + 0.00001) / aa_freqs[mut]))) %>% 
  select(-chain) %>% 
  rename(protein = pdb_id) %>% 
  drop_na()

models$`Baseline CNN` <- mutate(models$`Baseline CNN`, pred = as.integer(log2((pred + 0.00001) / aa_freqs[mut]))) %>% 
  select(-chain) %>% 
  rename(protein = pdb_id) %>% 
  drop_na()

models <- bind_rows(models, .id = 'model') %>%
  pivot_wider(names_from = 'model', values_from = 'pred') %>%
  left_join(blosum, by = c('wt', 'mut')) %>%
  pivot_longer(c(SPBuild, UNET, BLOSUM62, `PreGraph UNET`, `Baseline CNN`, `ESM-1b Logits`, `ESM-1b Top Model`),
               names_to = 'model', values_to = 'pred') %>%
  rename(true = ProteinNet) %>%
  mutate(diff = pred - true)
write_tsv(models, "data/pssm/combined_preds.tsv")

### Analyse ###
plots <- list()
plots$overall <- (ggplot(models, aes(x = diff, y = ..prop.., fill = model)) +
  geom_bar(position = 'dodge') +
  labs(x = 'Pred - ProteinNet', y = 'Proportion') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill=guide_legend(title = ''))) %>%
  labeled_plot(units = 'cm', width = 15, height = 15)

plots$per_aa <- (ggplot(models, aes(x = diff, y = ..prop.., fill = model)) +
  facet_wrap(~wt, nrow = 4) +
  geom_bar(position = 'dodge') +
  labs(x = 'Pred - ProteinNet', y = 'Proportion') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill=guide_legend(title = ''))) %>%
  labeled_plot(units = 'cm', width = 25, height = 20)

plots$scatter <- (count(models, model, true, pred) %>%
  group_by(model) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = true, y = pred, size = prop, colour = model)) +
  facet_wrap(~model) +
  geom_point(shape = 15) +
  geom_abline(slope = 1) +
  coord_fixed() +
  labs(x = 'ProteinNet Log Odds', y = 'Predicted Log Odds') +
  scale_colour_brewer(name = 'Model', type = 'qual', palette = 'Dark2') +
  scale_size(name = 'Proportion of\nPredictions') + 
  theme(legend.position = 'top')) %>%
  labeled_plot(units = 'cm', width = 30, height = 10)

model_summary <- drop_na(models) %>%
  mutate(abs_diff = abs(diff)) %>%
  group_by(protein, position, wt, mut) %>%
  mutate(best_diff = min(abs_diff)) %>%
  ungroup() %>%
  mutate(good_model = abs_diff < 2,
         best_model = diff == best_diff) %>%
  group_by(model) %>%
  summarise(n_tot = n(),
            n_best = sum(best_model) / n_tot,
            n_good = sum(good_model) / n_tot,
            n_best_good = sum(best_model & good_model) / n_tot,
            .groups = 'drop') %>%
  select(-n_tot) %>%
  pivot_longer(-model, names_prefix = 'n_', names_to = 'type', values_to = 'prop')

plots$closest <- mutate(model_summary, type = factor(type, levels = c('good', 'best', 'best_good'))) %>%
  ggplot(aes(x = type, y = prop, fill = model)) +
  geom_col(position = 'dodge') +
  scale_fill_brewer(name = 'Model', type = 'qual', palette = 'Dark2') +
  labs(x = '', y = 'Proportion') +
  scale_x_discrete(labels = c(good='|Pred - True| <= 1', best='Best Model', best_good='Both')) +
  theme(axis.text.x = element_markdown(),
        axis.ticks.x = element_blank())

plots$combined <- ggarrange(plots$closest + guides(fill = "none"), plots$scatter,
                            ncol = 1, nrow = 2, common.legend = TRUE, labels = 'AUTO', heights = c(1, 2)) %>%
  labeled_plot(height = 20, width = 30, units = 'cm')

### Save Plots ###
save_plotlist(plots, 'figures/pssm_predictions', overwrite='all')


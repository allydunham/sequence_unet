#!/usr/bin/env Rscript
# Produce figure 2 - PSSM prediction performance
source('src/config.R')

### Panels - Example output ###
pssm_preds <- bind_rows(
  true = read_tsv('data/pssm/pn_casp12_validation.tsv') %>%
    extract(protein, into = c("pdb_id", "chain"), regex = "[0-9]*#([0-9A-Z]*)_[0-9]*_([A-Za-z0-9])") %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "pred") %>%
    drop_na(),
  pred = read_tsv('data/pssm/unet_sequence_validation.tsv'),
  .id = "model"
) %>%
  filter(pdb_id == "1WAZ") %>%
  pivot_wider(names_from = 'model', values_from = 'pred') %>%
  mutate(diff = pred - true)

p_preds <- select(pssm_preds, -diff) %>%
  pivot_longer(cols = c(true, pred)) %>%
  mutate(name = factor(name, levels = c("true", "pred"))) %>%
  ggplot(aes(x = position, fill = value)) +
  facet_wrap(~name, ncol = 1, scales = "free_x", labeller = labeller(name = c(true = "True", pred = "Predicted"))) +
  geom_tile(aes(y = mut), colour = "grey") +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_distiller(name = "Frequency", palette = "PuRd", direction = 1, limits = c(0, 1)) +
  scale_x_continuous(expand = expansion(0), breaks = c(1, 10, 20, 30, 40, 46)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), title.vjust = 1)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())

p_pred_diff <- ggplot(pssm_preds, aes(x = position, fill = diff)) +
  geom_tile(aes(y = mut), colour = "grey") +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "Difference", low = "#b2182b", high = "#2166ac", limits = c(-0.8, 0.8)) +
  scale_x_continuous(expand = expansion(0), breaks = c(1, 10, 20, 30, 40, 46)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), title.vjust = 1)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        axis.text.x = element_text())

### Panel - Correlation with true values ###
pssm_models <- read_tsv("data/pssm/combined_preds.tsv") %>%
  mutate(model = factor(model, levels = c("BLOSUM62", "SPBuild", "ESM-1b", "Baseline CNN", "UNET", "PreGraph UNET")))

pssm_cor <- group_by(pssm_models, model) %>%
  group_modify(~broom::tidy(cor.test(.$pred, .$true)))

p_pssm_cor <- ggplot(pssm_cor, aes(x = model, y = estimate, fill = model, ymin = conf.low, ymax = conf.high)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(width = 0.5) +
  scale_fill_manual(name = 'Model', values = TOOL_COLOURS) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.5)) +
  labs(y = expression("Pearson's"~rho)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank())

### Panel - Proportion best ###
pssm_summary <- drop_na(pssm_models) %>%
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
  pivot_longer(-model, names_prefix = 'n_', names_to = 'type', values_to = 'prop') %>%
  mutate(type = factor(type, levels = c('good', 'best', 'best_good')))

p_pssm_summary <- ggplot(pssm_summary, aes(x = type, y = prop, fill = model)) +
  geom_col(position = 'dodge', width = 0.75) +
  scale_fill_manual(name = '', values = TOOL_COLOURS, guide = "none") +
  labs(x = '', y = 'Proportion') +
  scale_x_discrete(labels = c(good='|Pred - True| ≤ 1', best='Best Model', best_good='Both')) +
  scale_y_continuous(expand = expansion(0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank())

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_preds + labs(tag = 'A') + size
p2 <- p_pred_diff + labs(tag = 'B') + size
p3 <- p_pssm_cor + labs(tag = 'C') + size
p4 <- p_pssm_summary + labs(tag = 'D') + size

figure2 <- multi_panel_figure(width = c(100, 50), height = rep(20,6), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1:4, column = 1) %>%
  fill_panel(p2, row = 5:6, column = 1) %>%
  fill_panel(p3, row = 1:3, column = 2) %>%
  fill_panel(p4, row = 4:6, column = 2)

ggsave('figures/figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2),
       units = 'mm', dpi = 600)
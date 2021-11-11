#!/usr/bin/env Rscript
# Produce figure 2 - PSSM prediction performance
source('src/config.R')

### Panel - Example output ###
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

p_preds <- ggplot(pssm_preds, aes(x = position, fill = diff)) +
  geom_tile(aes(y = mut), colour = "grey") +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "Pred - True", low = "#b2182b", high = "#2166ac") +
  scale_x_continuous(expand = expansion(0), breaks = c(1, 10, 20, 30, 40, 46)) +
  guides(fill = guide_colourbar(barwidth = unit(5, "mm"), title.vjust = 1)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 10))

### Panel - Correlation with true values ###
pssm_models <- read_tsv("data/pssm/combined_preds.tsv") %>%
  mutate(model = factor(model, levels = c("BLOSUM62", "SPBuild", "Baseline CNN", "UNET", "PreGraph UNET")))

pssm_cor <- group_by(pssm_models, model) %>%
  group_modify(~broom::tidy(cor.test(.$pred, .$true)))

p_pssm_cor <- ggplot(pssm_cor, aes(x = model, y = estimate, fill = model, ymin = conf.low, ymax = conf.high)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(width = 0.5) +
  coord_flip() +
  scale_fill_manual(name = 'Model', values = tool_colours) +
  labs(y = "Pearson Correlation Coefficient") +
  theme(axis.title.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank())

### Panel - Proportion best ###
pssm_summary <- drop_na(pssm_models) %>%
  mutate(abs_diff = abs(diff)) %>%
  group_by(pdb_id, chain, position, wt, mut) %>%
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
  coord_flip() +
  scale_fill_manual(name = '', values = tool_colours, guide = FALSE) +
  labs(x = '', y = 'Proportion of Predictions') +
  scale_x_discrete(labels = c(good='|Pred - True| ≤ 1', best='Best Model', best_good='Both')) +
  theme(axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        legend.key.size = unit(2, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"))

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure2 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2),
       units = 'mm')
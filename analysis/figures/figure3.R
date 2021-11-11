#!/usr/bin/env Rscript
# Produce figure 3 - Frequency Prediction
source('src/config.R')

### Panel - Example Output ###
preds <- bind_rows(
  true = read_tsv('data/pssm/pn_casp12_validation.tsv') %>%
    extract(protein, into = c("pdb_id", "chain"), regex = "[0-9]*#([0-9A-Z]*)_[0-9]*_([A-Za-z0-9])") %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "pred") %>%
    drop_na(),
  pred = read_tsv('data/freq/unet_sequence_validation.tsv'),
  .id = "model"
) %>%
  filter(pdb_id == "2EK0") %>%
  pivot_wider(names_from = 'model', values_from = 'pred') %>%
  mutate(del = true < 0.01,
         correct = (pred > 0.5) == del)

p_preds <- ggplot(preds, aes(x = position, fill = pred)) +
  geom_tile(aes(y = mut), colour = "grey") +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "P(Del)", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1), midpoint = 0.5, low = "#2166ac", high = "#b2182b") +
  scale_x_continuous(expand = expansion(0), breaks = c(1, seq(10, 90, 10))) +
  guides(fill = guide_colourbar(barwidth = unit(5, "mm"), title.vjust = 1),
         colour = guide_legend(override.aes = list(fill = NA))) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right")

p_preds_box <- ggplot(preds, aes(x = pred > 0.44, y = true, fill = pred > 0.44)) + 
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(method = "wilcox", comparisons = list(c("TRUE", "FALSE"))) +
  labs(x = "Prediction", y = "True Frequency") +
  scale_x_discrete(labels = c(`FALSE` = expression(f >= 0.01), `TRUE` = expression(f < 0.01))) +
  scale_fill_manual(values = c(`FALSE` = "#2166ac", `TRUE` = "#b2182b"))

### Panel - Thresholds ###
training_logs <- read_tsv("data/freq/training_logs.tsv") %>%
  extract(model, into = c("experiment", "model", "dataset"), regex = "models/classifier/([^/]*)/([^/]*)/([^/]*)")

### Panel - Overall performance ###

### Panel - ROC curve ###

### Panel - PR Curve ###

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure3 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure3.pdf', figure3, width = figure_width(figure3), height = figure_height(figure3),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure3.png', figure3, width = figure_width(figure3), height = figure_height(figure3),
       units = 'mm')
#!/usr/bin/env Rscript
# Produce combined figure 2/3 - PSSM and classification performance
source('src/config.R')
source("src/analysis.R")

### Panel - Correlation with true values ###
pssm_models <- read_tsv("data/pssm/combined_preds.tsv") %>%
  mutate(model = factor(model, levels = c("BLOSUM62", "SPBuild", "ESM-1b Logits", "ESM-1b Top Model", "Baseline CNN", "UNET", "PreGraph UNET")))

pssm_cor <- group_by(pssm_models, model) %>%
  group_modify(~broom::tidy(cor.test(.$pred, .$true)))

p_pssm_cor <- ggplot(pssm_cor, aes(x = model, y = estimate, fill = model, ymin = conf.low, ymax = conf.high)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(width = 0.5) +
  coord_flip() +
  scale_fill_manual(name = 'Model', values = TOOL_COLOURS) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.5)) +
  labs(y = expression("Pearson's"~rho)) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"))

### Panel - ROC curve ###
classifier_roc <- read_tsv("data/freq/roc.tsv") %>%
  group_by(model) %>%
  mutate(pr_auc = pr_auc(tpr, precision),
         n = tp + tn + fp + fn) %>%
  ungroup() %>%
  distinct(model, auc, pr_auc, n) %>%
  mutate(model = factor(model, levels = model[order(auc)])) %>%
  pivot_longer(c(auc, pr_auc), names_to = "metric", values_to = "value")

p_roc <- ggplot(classifier_roc, aes(x = model, y = value, fill = model)) +
  facet_wrap(~metric, labeller = labeller(metric = c(auc="ROC AUC", pr_auc="PR AUC")), strip.position = "bottom") +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(name = "", values = TOOL_COLOURS[names(TOOL_COLOURS) %in% classifier_roc$model]) +
  theme(strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"))

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_pssm_cor + labs(tag = 'A') + size
p2 <- p_roc + labs(tag = 'B') + size

figure2_3 <- multi_panel_figure(width = 120, height = c(60, 60), columns =  1,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 2, column = 1)

ggsave('figures/figures/figure2_3.pdf', figure2_3, width = figure_width(figure2_3), height = figure_height(figure2_3),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure2_3.png', figure2_3, width = figure_width(figure2_3), height = figure_height(figure2_3),
       units = 'mm', dpi = 600)

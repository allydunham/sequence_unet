#!/usr/bin/env Rscript
# Produce figure 4 - Generalisation
source('src/config.R')

### Panel - Example ###

### Panel - Performance overview ###

### Panel - ROC curve ###

### Panel - Other generalisation datasets ###

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure4 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm')
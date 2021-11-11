#!/usr/bin/env Rscript
# Produce figure 5 - Proteome scale predictions
source('src/config.R')

### Panel - Speed comparison (maybe earlier?) ###

### Panel - Abundance correlations plus phylogeny ###

### Panel - Correlation with length variance ###

### Panel - Correlation with abundance variance ###

### Panel - Mycoplasma SIFT4G summary ###

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure5 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure5.pdf', figure5, width = figure_width(figure5), height = figure_height(figure5),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure5.png', figure5, width = figure_width(figure5), height = figure_height(figure5),
       units = 'mm')
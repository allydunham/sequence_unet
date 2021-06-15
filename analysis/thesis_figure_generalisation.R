#!/usr/bin/env Rscript
# Generate thesis figure on model generalisation
source('src/config.R')
source("src/analysis.R")

#### Panel - Top model ####
p_model <- blank_plot("Top Model Strucutre")

#### Panel - ClinVar ROC ####
p_clinvar <- blank_plot("ClinVar ROC")

#### Panel - Jelier ROC ####
p_jelier <- blank_plot("Jelier ROC")

#### Panel - DMS ROC ####
p_dms <- blank_plot("DMS ROC")

#### Figure Assembly ####
size <- theme(text = element_text(size = 9))
p1 <- p_model + labs(tag = 'A') + size
p2 <- p_clinvar + labs(tag = 'A') + size
p3 <- p_jelier + labs(tag = 'A') + size
p4 <- p_dms + labs(tag = 'A') + size

figure <- multi_panel_figure(width = 180, height = 240, columns = 2, rows = 2,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)

ggsave('figures/thesis_figure_generalisation.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_generalisation.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
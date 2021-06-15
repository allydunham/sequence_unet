#!/usr/bin/env Rscript
# Generate thesis figure on the basic model
source('src/config.R')
source("src/analysis.R")

#### Panel - Model Structure ####
p_model <- blank_plot("Model Strucutre")

#### Panel - Training ####
p_train <- blank_plot("Training Trace")

#### Panel - Classifier Performance ####
p_classifier <- blank_plot("Freq Classification ROC")

#### Panel - PSSM Prediction Performance ####
p_pssm <- blank_plot("PSSM Prediction Performance")

#### Figure Assembly ####
size <- theme(text = element_text(size = 9))
p1 <- p_model + labs(tag = 'A') + size
p2 <- p_train + labs(tag = 'B') + size
p3 <- p_classifier + labs(tag = 'C') + size
p4 <- p_pssm + labs(tag = 'D') + size

figure <- multi_panel_figure(width = 180, height = 240, columns = 2, rows = 3,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 3, column = 1:2)

ggsave('figures/thesis_figure_model.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_model.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

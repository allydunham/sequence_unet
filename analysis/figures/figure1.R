#!/usr/bin/env Rscript
# Produce figure 1 - Model design and training
source('src/config.R')
library(png)

### Panel - Overall network architecture ###
p_schematic <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 890), clip = 'off') +
  annotation_raster(readPNG("figures/model_schematic.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 890)

### Panel - CNN details ###

### Panel - GraphCNN details ###

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure1 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure1.pdf', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure1.png', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm')
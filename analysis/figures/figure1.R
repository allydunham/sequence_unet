#!/usr/bin/env Rscript
# Produce figure 1 - Model design and training
source('src/config.R')
library(png)

### Panel - Overall network architecture ###
p_schematic <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 890), clip = 'off') +
  annotation_raster(readPNG("figures/model_schematic.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 890)

### Panel - CNN details ###
p_cnn <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 1000), clip = 'off') +
  annotation_raster(readPNG("figures/model_cnn.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 1000)

### Panel - Information flow ###
p_information <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 1000), clip = 'off') +
  annotation_raster(readPNG("figures/model_receptive_field.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 1000)

### Panel - GraphCNN details ###
p_graphcnn <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 2000), clip = 'off') +
  annotation_raster(readPNG("figures/model_graph.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)

### Panel - Speed Comparison ###
tool_speed <- tibble(tool = c("FoldX", "SIFT4G", "UNET (inc. initialisation)", "UNET (prediction only)"),
                     tool_id = c("FoldX", "SIFT4G", "UNET", "UNET"),
                     time_lab = c("232h", "22m", "85s", "0.42s"), 
                     speed = c(835200000, 1320000, 85000, 420))

p_speed <- ggplot(tool_speed, aes(x = tool, y = speed, fill = tool_id, label = time_lab)) + 
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_y_continuous(breaks = c(1, 1000, 60*1000, 60*60*1000, 24*60*60*1000, 7*24*60*60*1000),
                     labels = c("0", "1 Second", "1 Minute", "1 Hour", "1 Day", "1 Week"),
                     expand = expansion(mult = c(0, 0.12)), trans = scales::pseudo_log_trans(sigma = 100, base = 10)) +
  labs(x = "", y = "SARS-CoV-2 Spike Prediction Time") +
  scale_fill_manual(values = TOOL_COLOURS) +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_schematic + labs(tag = 'A') + size
p2 <- p_cnn + labs(tag = 'B') + size
p3 <- p_information + labs(tag = 'C') + size
p4 <- p_graphcnn + labs(tag = 'D') + size
p5 <- p_speed + labs(tag = 'E') + size

figure1 <- multi_panel_figure(width = c(75, 75), height = c(65, 75/2, 75/2, 30), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 3, column = 1) %>%
  fill_panel(p4, row = 2:3, column = 2) %>%
  fill_panel(p5, row = 4, column = 1:2)

ggsave('figures/figures/figure1.pdf', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure1.png', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm')


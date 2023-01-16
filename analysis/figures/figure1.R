#!/usr/bin/env Rscript
# Produce figure 1 - Model design and training
source('src/config.R')
library(png)

### Panel - Overall network architecture ###
p_schematic <- ggplot() +
  coord_fixed(expand = FALSE, xlim = c(0, 2000), ylim = c(0, 890), clip = 'off') +
  annotation_raster(readPNG("figures/model_schematic.png"), xmin = 0, xmax = 2000, ymin = 0, ymax = 890)

### Panel - Speed Comparison ###
tool_speed <- tibble(class = c("FoldX", "SIFT4G", "Initialisation", "Prediction"),
                     tool = c("FoldX", "SIFT4G", "UNET", "UNET"),
                     time_lab = c("232h", "22m", "85s", "0.42s"), 
                     speed = c(835200000, 1320000, 85000, 420))

tool_cols <- c(Initialisation="#cab2d6",
               Prediction="#6a3d9a",
               SIFT4G="#e6ab02",
               FoldX="#66c2a5")
p_speed <- ggplot(mapping = aes(x = tool, y = speed, fill = class, label = time_lab)) + 
  geom_col(data = filter(tool_speed, tool == "UNET"), width = 0.4, show.legend = TRUE, position = position_identity()) +
  geom_col(data = filter(tool_speed, tool != "UNET"), width = 0.4, show.legend = FALSE, position = position_identity()) +
  geom_text(data = tool_speed, hjust = -0.2, size = 2.5) +
  coord_flip() +
  scale_y_continuous(breaks = c(1, 1000, 60*1000, 60*60*1000, 24*60*60*1000),
                     labels = c("0s", "1s", "1m", "1h", "24h"),
                     expand = expansion(mult = c(0, 0.12)), trans = scales::pseudo_log_trans(sigma = 100, base = 10)) +
  labs(x = "", y = "") +
  scale_fill_manual(name = "", values = tool_cols, breaks = c("Initialisation", "Prediction")) +
  guides(fill = guide_legend(keywidth = unit(3, "mm"), keyheight = unit(3, "mm"), reverse = TRUE)) +
  theme(axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-10, 0, -5, 0),
        legend.margin = margin())

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_schematic + labs(tag = 'A') + size
p2 <- p_speed + labs(tag = 'B') + size

figure1 <- multi_panel_figure(width = 140, columns = 1, height = c(65, 30), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 2, column = 1)

ggsave('figures/figures/figure1.pdf', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure1.png', figure1, width = figure_width(figure1), height = figure_height(figure1),
       units = 'mm', dpi = 600)


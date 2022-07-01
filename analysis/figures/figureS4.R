#!/usr/bin/env Rscript
# Generate Figure S2 testing how much padding varies predictions
source('src/config.R')
source("src/analysis.R")

padding <- bind_rows(
  freq = read_tsv("data/freq/padding_variation.tsv"),
  pssm = read_tsv("data/pssm/padding_variation.tsv"),
  .id = "model"
) %>% 
  pivot_longer(pad_32:pad_1024, names_to = "pad", values_to = "pred", names_prefix = "pad_", names_transform = list(pad = as.integer)) %>%
  group_by(model, id) %>%
  mutate(length = max(position)) %>%
  ungroup() %>%
  mutate(prop_pad = pad / length,
         diff = abs(pred - pad_0))

# Correlation
correlation <- group_by(padding, model, id, length, pad) %>%
  group_modify(~broom::tidy(cor.test(.x$pad_0, .x$pred))) %>%
  ungroup()

p_correlation <- ggplot(correlation, aes(x = length, y = estimate, group = id, colour = as.factor(pad))) +
  facet_wrap(~model, labeller = labeller(model = c(freq="Frequency Classifier", pssm="PSSM Predictor"))) +
  geom_path(colour = "black") +
  geom_point(shape = 20) +
  labs(x = "Protein Length", y = expression("Pearson's"~rho)) +
  scale_colour_manual(name = "Additional\nPadding", values = c("#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026"))

# Class Swaps
class_change <- filter(padding, model == "freq") %>% 
  mutate(class_0 = pad_0 > 0.5, class_pad = pred > 0.5) %>%
  group_by(id, length, pad) %>%
  summarise(n = n(), n_swap = sum(class_0 != class_pad), p = n_swap / n, .groups = "drop")

p_class_change <- ggplot(class_change, aes(x = length, y = p, group = id, colour = as.factor(pad))) +
  geom_path(colour = "black") +
  geom_point(shape = 20) +
  labs(x = "Protein Length", y = "Proportion of variants\nchanging class") +
  scale_colour_manual(name = "Additional\nPadding", values = c("#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026"))

# Difference
mean_diff <- group_by(padding, model, id, length, pad) %>%
  summarise(mean_diff = mean(diff), .groups = "drop")

p_diff <- ggplot(mean_diff, aes(x = length, y = mean_diff, group = id, colour = as.factor(pad))) +
  facet_wrap(~model, labeller = labeller(model = c(freq="Frequency Classifier", pssm="PSSM Predictor")), scales = "free_y") +
  geom_path(colour = "black") +
  geom_point(shape = 20) +
  labs(x = "Protein Length", y = "Mean Prediction\nDifference") +
  scale_colour_manual(name = "Additional\nPadding", values = c("#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026"))

# Accuracy
accuracy <- filter(padding, model == "freq") %>%
  select(id, position, wt, mut, length, freq, pad_0, pad, pred) %>%
  pivot_wider(names_from = pad, values_from = pred, names_prefix = "pad_") %>%
  pivot_longer(starts_with("pad"), names_to = "pad", values_to = "pred", names_prefix = "pad_", names_transform = list(pad = as.integer)) %>%
  mutate(end_dist = length - position,
         pos_cat = case_when(end_dist < 10 ~ "< 10", end_dist < 20 ~ "< 20", end_dist < 30 ~ "< 30", TRUE ~ "≥ 30")) %>%
  group_by(id, pos_cat, pad) %>%
  summarise(acc = sum((freq < 0.01 & pred > 0.5) | (freq >= 0.01 & pred <= 0.5)) / n(), .groups = "drop")

p_acc <- ggplot(accuracy, aes(x = pos_cat, y = acc, fill = as.factor(pad))) +
  geom_boxplot(position = position_dodge(), outlier.shape = 20) +
  geom_point(aes(colour = as.factor(pad)), alpha = 0, shape = 15, size = 2) +
  labs(x = "Distance to end of Protein", y = "Accuracy") +
  scale_fill_manual(name = "Additional\nPadding", values = c("#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026")) +
  scale_colour_manual(name = "Additional\nPadding", values = c("#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026"))

### Assemble Figure ###
size <- theme(text = element_text(size = 10))
p_legend <- get_legend(p_acc + guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1)), fill = "none") + 
                         theme(legend.position = "top") + size)
p1 <- p_correlation + guides(colour = "none") + labs(tag = 'A') + size
p2 <- p_class_change + guides(colour = "none") + labs(tag = 'B') + size
p3 <- p_diff + guides(colour = "none") + labs(tag = 'C') + size
p4 <- p_acc + guides(fill = "none", colour = "none") + labs(tag = 'D') + size

figure <- multi_panel_figure(width = 150, height = c(50, 50, 50, 50, 10), columns = 1, panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 3, column = 1) %>%
  fill_panel(p4, row = 4, column = 1) %>%
  fill_panel(p_legend, row = 5, column = 1)

ggsave('figures/figures/figureS2.pdf', figure, width = figure_width(figure), height = figure_height(figure),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figureS2.png', figure, width = figure_width(figure), height = figure_height(figure),
       units = 'mm', dpi = 600)


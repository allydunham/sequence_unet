#!/usr/bin/env Rscript
# Generate Figure S2 testing how much padding varies predictions
source('src/config.R')
source("src/analysis.R")

padding <- bind_rows(
  freq = read_tsv("data/freq/padding_variation.tsv"),
  pssm = read_tsv("data/pssm/padding_variation.tsv"),
  .id = "model"
) %>% 
  pivot_longer(pad_32:pad_1024, names_to = "pad", values_to = "pred", names_prefix = "pad_", names_transform = list(pad = as.integer))

correlation <- group_by(padding, model, id, pad) %>%
  group_modify(~broom::tidy(cor.test(.x$pad_0, .x$pred))) %>%
  ungroup()

p_correlation <- ggplot(correlation, aes(x = as.factor(pad), y = estimate, fill = model)) +
  geom_boxplot()

class_change <- filter(padding, model == "freq") %>% 
  mutate(class_0 = pad_0 > 0.5, class_pad = pred > 0.5) %>%
  group_by(id, pad) %>%
  summarise(n = n(), n_swap = sum(class_0 != class_pad), p = n_swap / n, .groups = "drop")

p_class_change <- ggplot(class_change, aes(x = as.factor(pad), y = p)) +
  geom_boxplot()

diff <- group_by(padding, model, id, pad) %>%
  summarise(mean_diff = mean(abs(pred - pad_0)), .groups = "drop")

p_diff <- ggplot(diff, aes(x = as.factor(pad), y = mean_diff)) +
  geom_boxplot()

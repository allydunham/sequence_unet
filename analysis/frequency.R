#!/usr/bin/env Rscript
# Analyse PSSM frequency distribution from ProteinNet
source("src/config.R")
plots <- list()

freqs <- read_tsv("data/validation_freqs.tsv") %>%
  pivot_longer(A:Y, names_to = "mut", values_to = "freq")

# Overall histogram
plots$histogram <- ggplot(freqs, aes(x = freq)) +
  geom_histogram(binwidth = 0.05) +
  labs(x = "PSSM Frequency", y = "Count")

# Threshold counts
counts <- tibble(threshold = c(0.1, 0.01, 0.001, 0.0001)) %>%
  mutate(Neutral = map_int(threshold, ~sum(freqs$freq >= .)),
         Deleterious = map_int(threshold, ~sum(freqs$freq < .)),
         ratio = Deleterious / Neutral)

plots$categories <- pivot_longer(counts, c(-threshold, -ratio), names_to = "type", values_to = "count") %>%
  ggplot(aes(x = threshold, y = count)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_x_log10(breaks = c(0.1, 0.01, 0.001, 0.0001), labels = c("0.1", "0.01", "0.001", "0.0001"),
                sec.axis = dup_axis(name = '', labels = signif(counts$ratio, digits = 3))) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set1") +
  labs(x = "Threshold", y = "Count") +
  theme(axis.ticks.x.top = element_blank())

save_plotlist(plots, "figures/proteinnet_frequency", overwrite = "all")

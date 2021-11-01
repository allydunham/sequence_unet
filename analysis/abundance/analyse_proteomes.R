#!/usr/bin/env Rscript
# Analyse proteomic measurements against Sequence UNET predictions
source("src/config.R")

# Import combined data

### Analysis ###
plots <- list()

cor_summary <- select(omics, organism, intensity, intensity_per_len, log2_fc, log2_fc_per_len, mean_sift, mean_pred, mean_wt, mean_mut) %>%
  bind_rows(., mutate(., organism = "Overall")) %>%
  pivot_longer(intensity:log2_fc_per_len, names_to = "abundance_type", values_to = "abundance") %>%
  pivot_longer(mean_sift:mean_mut, names_to = "pred_type", values_to = "pred") %>%
  filter(!(organism == "Overall" & pred_type == "mean_sift")) %>%
  mutate(pred = ifelse(pred_type == "mean_sift", -pred, pred)) %>%
  drop_na() %>%
  group_by(organism, abundance_type, pred_type) %>%
  filter(n() > 1) %>%
  group_modify(~tidy(cor.test(.$abundance, .$pred))) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

plots$intensity_type <- filter(cor_summary, pred_type == "mean_mut") %>%
  ggplot(aes(x = abundance_type, y = estimate)) +
  geom_boxplot() +
  labs(x = "Abundance Measurement", y = "Correlation estimate")

cor_summary_reduced <- filter(cor_summary, abundance_type == "log2_fc_per_len", pred_type %in% c("mean_mut", "mean_sift")) %>%
  mutate(tool = ifelse(pred_type == "mean_sift", "SIFT4G", "UNET"),
         organism = ifelse(tool == "SIFT4G", str_c("SIFT4G (", organism,")"), organism))

organism_levels <- c(sort(unique(cor_summary_reduced$organism) %>% magrittr::extract(!. %in% c("SIFT4G (Saccharomyces cerevisiae)", "Overall"))),
                     "Overall",
                     "SIFT4G (Saccharomyces cerevisiae)")

plots$correlations <- (ggplot(cor_summary_reduced, aes(x = factor(organism, levels = organism_levels), y = estimate,
                                                      ymin = conf.low, ymax = conf.high, fill = tool)) +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c(UNET="#377eb8", SIFT4G="#4daf4a")) +
  labs(x = "", y = expression("Pearsons's"~rho)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank(),
        text = element_text(size = 11))) %>%
  labeled_plot(unit = "cm", height = 15, width = 8)

### Save plots ###
save_plotlist(plots, "figures/abundance", overwrite = "all")

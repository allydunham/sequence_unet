#!/usr/bin/env Rscript
# Analyse proteomic measurements against Sequence UNET predictions
source("src/config.R")
library(ggtree)

# Import combined data
phyl <- ggtree::read.tree("data/abundance/muller_species.phy")
species <- readxl::read_xlsx("data/abundance/muller_species.xlsx", sheet = 2) %>%
  mutate(taxid = factor(as.character(taxid), levels = phyl$tip.label)) %>%
  arrange(taxid)

taxid_order <- levels(species$taxid)
species_order <- species$species

col_types <- cols(.default = col_double(), superkingdom = col_character(), organism = col_character(), taxid = col_character(),
                  source = col_character(), protein = col_character(), length = col_integer(), tool = col_character())
omics <- read_tsv("data/abundance/muller_proteomics_processed.tsv", col_types = col_types) %>%
  mutate(organism = factor(organism, levels = species_order),
         taxid = factor(taxid, levels = taxid_order),
         superkingdom = factor(superkingdom, levels = c("Eukaryote", "Archea", "Bacteria")))

### Analysis ###
plots <- list()

plots$length_intensity <- (ggplot(distinct(omics, protein, length, intensity), aes(x = cut_number(length, 10), y = intensity)) +
  geom_boxplot() +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x) +
  labs(x = "Length", y = "Intensity")) %>%
  labeled_plot(file_format = "png")

plots$kingdom_length <- ggplot(distinct(omics, superkingdom, organism, protein, length), aes(x = superkingdom, y = length, fill = superkingdom)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  stat_compare_means(method = "t.test", comparisons = list(c("Eukaryote", "Archea"), c("Archea", "Bacteria"), c("Eukaryote", "Bacteria"))) +
  labs(x = "", y = "Length")

get_cor <- function(x, ...) {
  if (nrow(x) < 3) {
    return(tibble(estimate=NA, statistic=NA, p.value=NA, parameter=NA, method=NA, alternative=NA))
  }
  return(tidy(cor.test(x$abundance, x$pred)))
}

# cor_summary <- select(omics, tool, organism, source, intensity:intensity_fc_per_len, mean_mut:mean_wt) %>%
#   bind_rows(., mutate(., organism = "Overall")) %>%
#   pivot_longer(intensity:intensity_fc_per_len, names_to = "abundance_type", values_to = "abundance") %>%
#   pivot_longer(mean_mut:mean_wt, names_to = "pred_type", values_to = "pred") %>%
#   filter(!(organism == "Overall" & tool == "SIFT4G")) %>%
#   drop_na() %>%
#   group_by(tool, organism, source, abundance_type, pred_type) %>%
#   group_modify(get_cor) %>%
#   ungroup() %>%
#   mutate(p.adj = p.adjust(p.value, method = "fdr")) %T>%
#   write_tsv("data/abundance/pred_correlations.tsv")
cor_summary <- read_tsv("data/abundance/pred_correlations.tsv")

plots$metric_cor_summary <- (ggplot(cor_summary, aes(x = tool, y = estimate, fill = source)) +
  facet_grid(rows = vars(abundance_type), cols = vars(pred_type)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5) +
  labs(x = "", y = expression("Pearson's"~rho)) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))) %>%
  labeled_plot(unit = "cm", height = 20, width = 20)

cor_summary_reduced <- filter(cor_summary, abundance_type == "log2_fc_per_len", pred_type %in% c("mean_mut", "mean_sift")) %>%
  mutate(tool = ifelse(pred_type == "mean_sift", "SIFT4G", "UNET"),
         organism = ifelse(tool == "SIFT4G", str_c("SIFT4G (", organism,")"), organism))

plots$metric_cor_summary_reduced <- (filter(cor_summary, organism %in% c("Homo sapiens", "Saccharomyces cerevisiae", "Escherichia coli")) %>%
  ggplot(aes(x = tool, y = estimate, fill = tool)) +
  facet_grid(rows = vars(abundance_type), cols = vars(pred_type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(subtitle = "Filtered to H. sapiens, S. cerevisiae, & E. coli") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))) %>%
  labeled_plot(unit = "cm", height = 10, width = 10)

# TODO compare TrEMBL vs SwissProt
omics_summary <- select(omics, superkingdom, organism, taxid, source, protein, tool, length, intensity_fc_per_len,
                        mean_conserved, percent_n_conserved, mean_mut) %>%
  filter(tool %in% c("SIFT4G", "UNET PSSM")) %>%
  group_by(superkingdom, organism, taxid, tool) %>%
  summarise(tidy(cor.test(intensity_fc_per_len, mean_conserved)),
            mean_length = mean(length),
            mean_mut = mean(mean_mut),
            mean_percent_conserved = mean(percent_n_conserved),
            .groups = "drop")

plots$correlations <- ggplot(omics_summary, aes(x = organism, y = estimate, ymin = pmax(conf.low, 0), ymax = conf.high, fill = superkingdom)) +
  facet_grid(cols = vars(tool), scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.7) +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(expand = expansion(0)) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "", y = expression("Pearsons's"~rho)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        legend.position = "bottom")
plots$correlations <- labeled_plot(plots$correlations, unit = "cm", height = 10, width = 20)

plots$conservation <- ggplot(omics, aes(x = organism, y = mean_conserved, fill = superkingdom)) +
  facet_wrap(~tool, ncol = 1) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  labs(x = "", y = "Mean Deleterious Variants") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plots$conservation <- labeled_plot(plots$conservation, unit = "cm", height = 20, width = 10)

plots$kingdom_conservation <- ggplot(omics_summary, aes(x = superkingdom, fill = superkingdom, y = mean_percent_conserved)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.5) +
  stat_compare_means(method = "t.test", comparisons = list(c("Eukaryote", "Archea"), c("Archea", "Bacteria"), c("Eukaryote", "Bacteria"))) +
  labs(x = "", y = "Mean % Conserved")

plots$kingdom_mean_length <- ggplot(omics_summary, aes(x = superkingdom, y = mean_length, fill = superkingdom)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  stat_compare_means(method = "t.test", comparisons = list(c("Eukaryote", "Archea"), c("Archea", "Bacteria"), c("Eukaryote", "Bacteria"))) +
  labs(x = "", y = "Length")

plots$kingdom_mean_mut <- ggplot(omics, aes(x = superkingdom, y = mean_mut, fill = superkingdom)) +
  facet_wrap(~tool, nrow = 2, scales = "free_y") +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  stat_compare_means(method = "t.test", comparisons = list(c("Eukaryote", "Archea"), c("Archea", "Bacteria"), c("Eukaryote", "Bacteria"))) +
  labs(x = "", y = "Mean predicted score")
plots$kingdom_mean_mut <- labeled_plot(plots$kingdom_mean_mut, file_format = "png")
  
### Save plots ###
save_plotlist(plots, "figures/abundance", overwrite = "all")

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

comp_time <- bind_rows(
  read_tsv("data/abundance/times_other_tools.tsv"),
  read_tsv("data/abundance/times_cpu.tsv") %>%
    extract(id, "gene", "[a-z]*\\|([A-Z0-9]*)\\|.*") %>%
    mutate(tool = "UNET (CPU)", variants = length * 19),
  read_tsv("data/abundance/times_gpu.tsv") %>%
    extract(id, "gene", "[a-z]*\\|([A-Z0-9]*)\\|.*") %>%
    mutate(tool = "UNET (GPU)", variants = length * 19)
)

### Analysis ###
plots <- list()

# Tool efficiency
common_time_axis <- dup_axis(name = "", breaks = c(0.01, 1, 60, 60*60, 60*60*24),
                             labels = c("1 Millisecond", "1 Second", "1 Minute", "1 Hour", "1 Day"))

plots$comp_time <- ggplot(comp_time, aes(x = variants, y = time, colour = tool, size = tool, alpha = tool)) +
  geom_point(shape = 20) +
  scale_y_log10(labels = c("0.01", "0.1", "1", "10", "100", "10,000", "1,000,000"), breaks = c(0.01, 0.1, 1, 10, 100, 10000, 1000000),
                limits = c(0.01, 1000000), sec.axis = common_time_axis) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_brewer(name = "", palette = "Dark2") +
  scale_size_manual(values = c(SIFT4G = 1, FoldX = 1, `UNET (GPU)` = 0.1, `UNET (CPU)` = 0.1)) +
  scale_alpha_manual(values = c(SIFT4G = 1, FoldX = 1, `UNET (GPU)` = 0.5, `UNET (CPU)` = 0.5)) +
  guides(alpha = "none", size = "none") +
  labs(x = "Variants Computed", y = "Computation Time (s)") +
  theme(axis.ticks.y.right = element_blank())
plots$comp_time <- labeled_plot(plots$comp_time, file_format = "png", unit = "cm", height = 10, width = 10)

# Omics
plots$length_intensity <- distinct(omics, superkingdom, protein, length, intensity_fc) %>%
  {ggplot(., aes(x = length, y = intensity_fc, colour = superkingdom)) +
      facet_wrap(~superkingdom, scales = "free_x", ncol = 1) +
      geom_point(show.legend = FALSE, shape = 20, size = 0.2, alpha = 0.5) +
      geom_smooth(method = 'lm', formula = y ~ x, colour = "black") +
      labs(x = "Length", y = expression("log"[2]~"Intensity FC"))} %>%
  labeled_plot(file_format = "png")

plots$kingdom_length <- ggplot(distinct(omics, superkingdom, organism, protein, length), aes(x = superkingdom, y = length, fill = superkingdom)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  stat_compare_means(method = "t.test", comparisons = list(c("Eukaryote", "Archea"), c("Archea", "Bacteria"), c("Eukaryote", "Bacteria"))) +
  labs(x = "", y = "Length")

# Abundance vs preds
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

# Protein Source
source_summary <- select(omics, superkingdom, organism, source, tool, intensity_fc_per_len, mean_conserved) %>%
  filter(tool %in% c("SIFT4G", "UNET PSSM")) %>%
  drop_na() %>%
  group_by(superkingdom, source, organism, tool) %>%
  filter(n() > 2) %>%
  summarise(tidy(cor.test(intensity_fc_per_len, mean_conserved)), .groups = "drop")

plots$protein_source <- ggplot(source_summary, aes(x = source, y = estimate, fill = source)) +
  facet_wrap(~superkingdom, nrow = 1) +
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(method = "t.test", comparisons = list(c("SwissProt", "TrEMBL"))) +
  labs(x = "", y = expression("Pearson's"~rho)) +
  lims(y = c(-0.4, 1))

# Reduced 
omics_summary <- select(omics, superkingdom, organism, taxid, source, protein, tool, length, intensity, intensity_fc_per_len,
                        mean_conserved, percent_n_conserved, mean_mut) %>%
  filter(tool %in% c("SIFT4G", "UNET PSSM")) %>%
  group_by(superkingdom, organism, taxid, tool) %>%
  summarise(tidy(cor.test(intensity_fc_per_len, mean_conserved)),
            mean_intensity = mean(intensity),
            sd_intensity = sd(intensity),
            mean_length = mean(length),
            sd_length = sd(length),
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

plots$correlation_vs_length <- ggplot(filter(omics_summary, tool == "UNET PSSM"), aes(x = estimate, y = mean_length, colour = superkingdom)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  scale_colour_brewer(name = "", palette = "Dark2") +
  labs(x = expression("Pearson's"~rho), y = "Mean Protein Length")

plots$correlation_vs_abundance <- ggplot(filter(omics_summary, tool == "UNET PSSM"),
                                         aes(x = estimate, y = mean_intensity, colour = superkingdom)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  scale_y_log10() +
  scale_colour_brewer(name = "", palette = "Dark2") +
  labs(x = expression("Pearson's"~rho), y = "Mean Abundance")

plots$correlation_vs_length_sd <- ggplot(filter(omics_summary, tool == "UNET PSSM"), aes(x = estimate, y = sd_length)) +
  geom_point(aes(colour = superkingdom)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = superkingdom)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", size = 1) +
  scale_colour_brewer(name = "", palette = "Dark2") +
  labs(x = expression("Pearson's"~rho), y = "Protein Length SD")

plots$correlation_vs_abundance_sd <- ggplot(filter(omics_summary, tool == "UNET PSSM"),
                                         aes(x = estimate, y = sd_intensity)) +
  geom_point(aes(colour = superkingdom)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = superkingdom)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", size = 1) +
  scale_y_log10() +
  scale_colour_brewer(name = "", palette = "Dark2") +
  labs(x = expression("Pearson's"~rho), y = "Abundance SD")

### Save plots ###
save_plotlist(plots, "figures/abundance")

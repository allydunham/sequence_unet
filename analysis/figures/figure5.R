#!/usr/bin/env Rscript
# Produce figure 5 - Proteome scale predictions
source('src/config.R')

phyl <- phylogram::read.dendrogram("data/abundance/muller_species.phy")
phyl_data <- ggdendro::dendro_data(phyl, type = "rectangle")

species <- readxl::read_xlsx("data/abundance/muller_species.xlsx", sheet = 2) %>%
  mutate(taxid = factor(as.character(taxid), levels = phyl_data$labels$label)) %>%
  arrange(taxid) %>%
  drop_na(taxid)

taxid_order <- levels(species$taxid)
species_order <- species$species

col_types <- cols(.default = col_double(), superkingdom = col_character(), organism = col_character(), taxid = col_character(),
                  source = col_character(), protein = col_character(), length = col_integer(), tool = col_character())
omics <- read_tsv("data/abundance/muller_proteomics_processed.tsv", col_types = col_types) %>%
  mutate(organism = factor(organism, levels = species_order),
         taxid = factor(taxid, levels = taxid_order),
         superkingdom = factor(superkingdom, levels = c("Eukaryote", "Archea", "Bacteria")))

superkingdom <- distinct(omics, superkingdom, organism, taxid) %>%
  mutate(across(.fns = as.character))

### Panel - Speed comparison ###
comp_time <- bind_rows(
  read_tsv("data/abundance/times_other_tools.tsv"),
  read_tsv("data/abundance/times_cpu.tsv") %>%
    extract(id, "gene", "[a-z]*\\|([A-Z0-9]*)\\|.*") %>%
    mutate(tool = "UNET (CPU)", variants = length * 19),
  read_tsv("data/abundance/times_gpu.tsv") %>%
    extract(id, "gene", "[a-z]*\\|([A-Z0-9]*)\\|.*") %>%
    mutate(tool = "UNET (GPU)", variants = length * 19)
)

common_time_axis <- dup_axis(name = "", breaks = c(0.01, 1, 60, 60*60, 60*60*24),
                             labels = c("1 Millisecond", "1 Second", "1 Minute", "1 Hour", "1 Day"))

comp_time_colours <- c(SIFT4G = unname(TOOL_COLOURS["SIFT4G"]), FoldX = unname(TOOL_COLOURS["FoldX"]),
                       `UNET (CPU)` = "#6a3d9a", `UNET (GPU)` = "#e31a1c")

p_comp_time <- ggplot(comp_time, aes(x = variants, y = time, colour = tool, size = tool, alpha = tool)) +
  geom_point(shape = 20) +
  scale_y_log10(labels = c("0.01", "0.1", "1", "10", "100", "10,000", "1,000,000"), breaks = c(0.01, 0.1, 1, 10, 100, 10000, 1000000),
                limits = c(0.01, 1000000), sec.axis = common_time_axis) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_manual(name = "", values = comp_time_colours) +
  scale_size_manual(values = c(SIFT4G = 1, FoldX = 1, `UNET (GPU)` = 0.1, `UNET (CPU)` = 0.1)) +
  scale_alpha_manual(values = c(SIFT4G = 1, FoldX = 1, `UNET (GPU)` = 0.5, `UNET (CPU)` = 0.5)) +
  guides(alpha = "none", size = "none", colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Variants Computed", y = "Time (s)") +
  theme(axis.ticks.y.right = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-10, 0, 0, 0),
        legend.margin = margin(),
        legend.key.size = unit(2, "mm"))

### Panel - Abundance correlations plus phylogeny ###
omics_summary <- select(omics, superkingdom, organism, taxid, source, protein, tool, length, intensity, intensity_fc_per_len,
                        mean_conserved, percent_n_conserved, mean_mut) %>%
  filter(tool %in% c("SIFT4G", "UNET PSSM")) %>%
  group_by(superkingdom, organism, taxid, tool) %>%
  summarise(tidy(cor.test(intensity_fc_per_len, mean_conserved)),
            sd_intensity = sd(intensity, na.rm = TRUE),
            mean_intensity = mean(intensity, na.rm = TRUE),
            sd_length = sd(length, na.rm = TRUE),
            mean_length = mean(length, na.rm = TRUE),
            sd_mut = sd(mean_mut, na.rm = TRUE),
            mean_mut = mean(mean_mut, na.rm = TRUE),
            sd_percent_conserved = sd(percent_n_conserved, na.rm = TRUE),
            mean_percent_conserved = mean(percent_n_conserved, na.rm = TRUE),
            .groups = "drop")

phyl_labels <- left_join(phyl_data$labels,
                         filter(omics_summary, tool == "UNET PSSM") %>% 
                           select(label = taxid, superkingdom, organism, estimate, sd_intensity, mean_intensity, sd_length,
                                  mean_length, sd_mut, mean_mut, sd_percent_conserved, mean_percent_conserved),
                         by = "label") %>%
  as_tibble()

phyl_segments <- as_tibble(phyl_data$segments) %>%
  mutate(superkingdom = ifelse(xend <= max(phyl_labels[phyl_labels$superkingdom == "Eukaryote", "x"]),
                               "Eukaryote",
                               ifelse(xend <= max(phyl_labels[phyl_labels$superkingdom == "Archea", "x"]),
                                      "Archea",
                                      "Bacteria"))) %>%
  left_join(select(phyl_labels, -superkingdom), by = c("xend" = "x", "yend" = "y"))

omics_unet <- filter(omics_summary, tool == "UNET PSSM") %>%
  mutate(organism_int = as.integer(organism))
omics_sift <- filter(omics_summary, tool == "SIFT4G") %>%
  mutate(organism_int = 105 + order(as.integer(organism)))

org_labels <- c(levels(omics_unet$organism), "", "", as.character(omics_sift$organism))

p_correlation <- ggplot() +
  annotate("text", x = c(42.5, 107), y = 0.55 + 0.05*7, label = c("UNET", "SIFT4G"), size = 3, vjust = -1.5) +
  
  # UNET correlation columns
  geom_col(data = omics_unet, mapping = aes(x = organism_int, y = estimate, fill = superkingdom), width = 0.7) +
  geom_errorbar(data = omics_unet, mapping = aes(x = organism_int, ymin = pmax(conf.low, 0), ymax = conf.high), width = 0.5) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(name = "", values = KINGDOM_COLOURS) +
  
  # SIFT4G correlation columns
  geom_col(data = omics_sift, mapping = aes(x = organism_int, y = estimate, fill = superkingdom), width = 0.7) +
  geom_errorbar(data = omics_sift, mapping = aes(x = organism_int, ymin = pmax(conf.low, 0), ymax = conf.high), width = 0.5) +
  
  # Phylogeny
  geom_segment(data = phyl_segments, aes(x = x, y = 0.55 + 0.05*y, xend = xend, yend = 0.55 + 0.05*yend, colour = superkingdom)) +
  annotate("text", x = c(18, 42.5, 79), y = 0.55 + 0.05*7, label = c("Eukaryote", "Archea", "Bacteria"), colour = KINGDOM_COLOURS,
           size = 2.5, vjust = -0.5) +
  scale_colour_manual(name = "", values = KINGDOM_COLOURS) +
  
  # Setup
  scale_x_continuous(breaks = 1:length(org_labels), labels = org_labels) +
  scale_y_continuous(expand = expansion(c(0, 0.13)), breaks = c(0, 0.25, 0.5)) +
  guides(fill = "none", colour = "none") +
  labs(x = "", y = expression("Pearsons's"~rho)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(0,0,0,0)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.2),
        legend.position = "bottom")

### Panel - Correlation with length variance ###
p_correlation_length_sd <- ggplot(filter(omics_summary, tool == "UNET PSSM"), aes(x = estimate, y = sd_length)) +
  geom_point(aes(colour = superkingdom), shape = 20, size = 1.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = superkingdom), size = 0.5) +
  scale_colour_manual(name = "", values = KINGDOM_COLOURS) +
  labs(x = expression("Pearson's"~rho), y = expression(sigma~"Length"))

### Panel - Correlation with abundance variance ###
p_correlation_abundance_sd <- ggplot(filter(omics_summary, tool == "UNET PSSM"), aes(x = estimate, y = sd_intensity)) +
  geom_point(aes(colour = superkingdom), shape = 20, size = 1.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = superkingdom), size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", size = 0.5) +
  scale_y_log10() +
  scale_colour_manual(name = "", values = KINGDOM_COLOURS) +
  labs(x = expression("Pearson's"~rho), y = expression(sigma~"Abunance"))

### Figure Assembly ###
size <- theme(text = element_text(size = 9))
p1 <- p_comp_time + labs(tag = 'A') + size
p2 <- p_correlation + labs(tag = 'B') + size
p3 <- p_correlation_abundance_sd + guides(colour = "none") + labs(tag = 'C') + size
p4 <- p_correlation_length_sd + guides(colour = "none") + labs(tag = 'D') + size
pleg <- {p_correlation_abundance_sd + 
    theme(legend.box.margin = margin(-10, 0, -10, 0), legend.margin = margin()) + 
    size} %>%
  get_legend("bottom") %>%
  as_ggplot()

figure5 <- multi_panel_figure(width = c(90, 90), height = c(40, 120, 40, 10), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1:2) %>%
  fill_panel(p3, row = 3, column = 1) %>%
  fill_panel(p4, row = 3, column = 2) %>%
  fill_panel(pleg, row = 4, column = 1:2)

#ggsave('figures/figures/figure5.pdf', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure5.png', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')

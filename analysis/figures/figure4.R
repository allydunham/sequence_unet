#!/usr/bin/env Rscript
# Produce figure 4 - Generalisation
source('src/config.R')
source("src/analysis.R")

model_type <- c(
  `UNET (Top)` = "UNET (All)", `UNET (Finetune)` = "UNET (All)", `UNET` = "UNET (All)",
  `UNET (Top)*` = "UNET (Structured)", `PreGraph UNET (Top)*` = "UNET (Structured)", `UNET (Finetune)*` = "UNET (Structured)",
  `PreGraph UNET (Finetune)*` = "UNET (Structured)", `UNET*` = "UNET (Structured)", `PreGraph UNET*` = "UNET (Structured)",
  `Baseline ClinVar*` = "Neural Network", `Baseline Frequency*` = "Neural Network",
  `Baseline ClinVar` = "Neural Network", `Baseline Frequency` = "Neural Network",
  `ESM-1v`="Neural Network", `EVE` = "Neural Network",
  `PrimateAI` = "Neural Network", `DANN` = "Neural Network", 
  `SIFT4G` = "MSA", `BLOSUM62` = "MSA", `PROVEAN` = "MSA", `MutationAssessor` = "MSA", `LIST S2` = "MSA", `GERP++` = "MSA", 
  `FoldX` = "Structure", 
  `ClinPred` = "Ensemble", `REVEL` = "Ensemble", `M-CAP` = "Ensemble",`MPC` = "Ensemble", `MVP` = "Ensemble",
  `PolyPhen2` = "Other ML", `VEST4` = "Other ML", `CADD` = "Other ML", `DEOGEN2` = "Other ML", 
  `FATHMM-XF` = "Other ML", `MutPred` = "Other ML", `FATHMM` = "Other ML", `FATHMM-MKL` = "Other ML",
  `GenoCanyon` = "Other ML"
)

type_colours <- c(
  `UNET (Structured)`="#6a3d9a",
  `UNET (All)`="#cab2d6",
  `Neural Network` = "#1f78b4",
  `Other ML`="#a6cee3",
  `MSA`="#33a02c",
  `Structure`="#a65628",
  Ensemble="#e6ab02"
)

### Panel - Example ###
preds <- read_tsv("data/clinvar/preds/pn_testing_features.tsv") %>%
  filter(pdb_id == "TBM#T0865")

p_preds <- ggplot(preds, aes(x = position, fill = pred)) +
  geom_raster(aes(y = mut)) +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "P(del)", low = "#2166ac", high = "#b2182b", midpoint = 0.5) +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

### Panel - ClinVar ROC/PR AUC Summary ###
clinvar_roc <- bind_rows(read_tsv("data/clinvar/roc.tsv"), read_tsv("data/clinvar/other_tools_roc.tsv")) %>%
  filter(!str_detect(model, "UNET Thresh"))

clinvar_auc <- mutate(clinvar_roc, n = tp + tn + fp + fn) %>%
  distinct(model, auc, pr_auc, n) %>%
  arrange(desc(auc)) %>%
  mutate(type = model_type[model],
         model = factor(model, levels = rev(model))) %>%
  arrange(model)

p_clinvar <- ggplot(clinvar_auc, aes(x = as.integer(model), y = auc, fill = type)) +
  geom_col(width = 0.5) +
  scale_fill_manual(name="", values = type_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 1)) +
  scale_x_continuous(labels = levels(clinvar_auc$model), breaks = 1:nrow(clinvar_auc),
                     sec.axis = dup_axis(labels = clinvar_auc$n, breaks = 1:nrow(clinvar_auc))) +
  coord_flip() +
  labs(x = "", y = "ClinVar AUC") +
  guides(fill = guide_legend(keywidth = unit(3, "mm"), keyheight = unit(3, "mm"))) +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-10, 0, -5, 0),
        legend.margin = margin())

### Panel - DMS Summary ###
dms_cor <- read_tsv("data/dms/correlation.tsv") %>%
  drop_na() %>%
  arrange(desc(mean_spearman)) %>%
  mutate(type = model_type[tool],
         tool = factor(tool, levels = rev(tool)))

p_dms <- ggplot(dms_cor, aes(x = as.integer(tool), y = mean_spearman, fill = type)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  geom_errorbar(aes(ymin = mean_spearman - stderr_spearman, ymax = mean_spearman + stderr_spearman), width = 0.5) +
  scale_fill_manual(values = type_colours) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.5)) +
  scale_x_continuous(labels = levels(dms_cor$tool), breaks = 1:nrow(dms_cor),
                     sec.axis = dup_axis(labels = dms_cor$variants, breaks = 1:nrow(dms_cor))) +
  labs(x = "", y = expression("Spearman's"~rho)) +
  coord_flip() +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

### Alt Panel - Specieswise DMS Summary ###
dms_cor_species <- read_tsv("data/dms/correlation_species.tsv") %>%
  drop_na() %>%
  filter(str_starts(tool, "UNET"))

p_dms_species <- ggplot(dms_cor_species, aes(x = pearson, y = spearman, colour = species, shape = tool)) +
  geom_point() +
  scale_color_brewer(name = "Species", palette = "Dark2") +
  scale_shape(name = "") +
  coord_fixed() +
  labs(x = expression("Pearson's"~rho), y = expression("Spearman's"~rho))
ggsave("figures/figures/figure_S3.pdf", p_dms_species, units = "cm", height = 12, width = 17)
ggsave("figures/figures/figure_S3.png", p_dms_species, units = "cm", height = 12, width = 17)

### Panel - Jelier Yeast Summary ###
jelier_roc <- read_tsv("data/jelier/roc.tsv") %>%
  mutate(n = tp + tn + fp + fn) %>%
  distinct(model, auc, n) %>%
  arrange(desc(auc)) %>%
  mutate(type = model_type[model],
         model = factor(model, levels = rev(model))) %>%
  arrange(model)

p_jelier <- ggplot(jelier_roc, aes(x = as.integer(model), y = auc, fill = type)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  scale_x_continuous(labels = levels(jelier_roc$model), breaks = 1:nrow(jelier_roc),
                     sec.axis = dup_axis(labels = jelier_roc$n, breaks = 1:nrow(jelier_roc))) +
  coord_flip() +
  scale_fill_manual(values = type_colours) +
  labs(y = "Jelier AUC") +
  lims(y = c(0,1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_preds + labs(tag = 'A') + size
p2 <- p_clinvar + guides(fill = "none") + labs(tag = 'B') + size
p3 <- p_dms + labs(tag = 'C') + size
p4 <- p_jelier + labs(tag = 'D') + size
p_legend <- as_ggplot(get_legend(p_clinvar))

figure4 <- multi_panel_figure(width = c(90, 90), height = c(45, 95, 50, 10), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2:3, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 3, column = 2) %>%
  fill_panel(p_legend, row = 4, column = 1:2)
  
ggsave('figures/figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', dpi = 600)

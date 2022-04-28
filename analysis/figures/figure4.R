#!/usr/bin/env Rscript
# Produce figure 4 - Generalisation
source('src/config.R')
source("src/analysis.R")

model_colours <- c(
  `UNET (Top)` = "#1f78b4", `PreGraph UNET (Top)` = "#a6cee3",
  `UNET (Finetune)` = "#e31a1c", `PreGraph UNET (Finetune)` = "#fb9a99",
  `UNET` = "#6a3d9a", `PreGraph UNET` = "#cab2d6",
  `Baseline ClinVar` = "#33a02c", `Baseline Frequency` = "#33a02c", 
  `ClinPred` = "#989898", `REVEL` = "#989898", `MVP` = "#989898", `VEST4` = "#989898", `M-CAP` = "#989898",
  `DEOGEN2` = "#989898", `CADD` = "#989898", `EVE` = "#e6ab02", `FATHMM-XF` = "#989898", `PolyPhen2` = "#989898",
  `PROVEAN` = "#989898", `MutationAssessor` = "#989898", `MutPred` = "#989898", `FATHMM` = "#989898", `ESM-1v`="#e6ab02",
  `SIFT4G` = "#e6ab02", `PrimateAI` = "#989898", `LIST S2` = "#989898", `MPC` = "#989898", `FATHMM-MKL` = "#989898",
  `FoldX` = "#e6ab02", `DANN` = "#989898", `GenoCanyon` = "#989898", `GERP++` = "#989898", `BLOSUM62` = "#989898"
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

clinvar_auc <- distinct(clinvar_roc, model, auc, pr_auc) %>%
  arrange(desc(auc)) %>%
  mutate(model = factor(model, levels = rev(model)))

p_clinvar <- ggplot(clinvar_auc, aes(x = model, y = auc, fill = model)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 1)) +
  coord_flip() +
  labs(x = "", y = "ClinVar AUC") +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

### Panel - DMS Summary ###
dms_cor <- read_tsv("data/dms/correlation.tsv") %>%
  drop_na() %>%
  arrange(desc(mean_spearman)) %>%
  mutate(tool = factor(tool, levels = rev(tool)))

p_dms <- ggplot(dms_cor, aes(x = tool, y = mean_spearman, fill = tool)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  geom_errorbar(aes(ymin = mean_spearman - stderr_spearman, ymax = mean_spearman + stderr_spearman), width = 0.5) +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.5)) +
  labs(x = "", y = expression("Spearman's"~rho)) +
  coord_flip() +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

### Panel - Jelier Yeast Summary ###
jelier_roc <- read_tsv("data/jelier/roc.tsv") %>%
  select(model, auc) %>%
  distinct() %>%
  arrange(desc(auc)) %>%
  mutate(model = factor(model, levels = rev(model)))

p_jelier <- ggplot(jelier_roc, aes(x = model, y = auc, fill = model)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = model_colours) +
  labs(y = "Jelier AUC") +
  lims(y = c(0,1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        axis.title.y = element_blank())

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_preds + labs(tag = 'A') + size
p2 <- p_clinvar + labs(tag = 'B') + size
p3 <- p_dms + labs(tag = 'C') + size
p4 <- p_jelier + labs(tag = 'D') + size

figure4 <- multi_panel_figure(width = c(90, 90), height = c(45, 95, 50), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2:3, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 3, column = 2)

ggsave('figures/figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', dpi = 600)

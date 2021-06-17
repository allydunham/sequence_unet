#!/usr/bin/env Rscript
# Generate thesis figure on model generalisation
source('src/config.R')
source("src/analysis.R")

#### Panel - Pathogenicity predictions ####
preds <- read_tsv("data/clinvar/preds/unet_freq_structure_features_top.tsv")

p_preds <- blank_plot("Pathogenicity predictions")

#### Panel - Performance? ####
p_performance <- blank_plot("General Performance")

#### Panel - ClinVar ROC ####
clinvar_roc <- read_tsv("data/clinvar/roc.tsv") %>%
  mutate(model_auc = auc_labeled_model(model, auc))

clinvar_auc <- distinct(clinvar_roc, model_auc, auc) %>%
  arrange(desc(auc)) %>%
  mutate(fpr = 1, tpr = rev(seq(from = 0.03, by = 0.05, length.out = n())))

p_clinvar <- ggplot(clinvar_roc, aes(x = fpr, y = tpr, colour = model_auc, label = model_auc)) +
  geom_step(show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = clinvar_auc, hjust = 1, show.legend = FALSE, size = 2.5) +
  coord_fixed() +
  scale_colour_brewer(name = "", type = "qual", palette = "Set1") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

#### Panel - Generalisation ####
gen_roc <- bind_rows(DMS = read_tsv("data/dms/roc.tsv"), Jelier = read_tsv("data/jelier/roc.tsv"), .id = "type") %>%
  select(type, model, auc) %>%
  distinct() %>%
  mutate(model = factor(model, levels = reorder(model[type == "Jelier"], -auc[type == "Jelier"])))

p_generalisation <- ggplot(gen_roc, aes(x = model, y = auc, fill = model)) +
  facet_wrap(~type, ncol = 1) +
  geom_col(show.legend = FALSE, width = 0.5) +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "AUC") +
  lims(y = c(0,1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        axis.title.y = element_blank())

#### Figure Assembly ####
size <- theme(text = element_text(size = 9))
p1 <- p_preds + labs(tag = 'A') + size
p2 <- p_performance + labs(tag = 'B') + size
p3 <- p_clinvar + labs(tag = 'C') + size
p4 <- p_generalisation + labs(tag = 'D') + size

figure <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)

ggsave('figures/thesis_figure_generalisation.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_generalisation.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
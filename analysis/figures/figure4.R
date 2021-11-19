#!/usr/bin/env Rscript
# Produce figure 4 - Generalisation
source('src/config.R')
source("src/analysis.R")

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

### Panel - Performance overview ###
clinvar_acc <- bind_rows(
  # Trained on ProteinNet
  `Baseline Frequency` = read_tsv("data/clinvar/preds/baseline_freq.tsv"),
  `UNET` = read_tsv("data/clinvar/preds/unet_freq.tsv"),
  `PreGraph UNET` = read_tsv("data/clinvar/preds/unet_freq_structure.tsv"),
  
  # Trained on ClinVar (incl. some ProteinNet pre-training for some models)
  `Baseline ClinVar` = read_tsv("data/clinvar/preds/baseline_clinvar.tsv"),
  `UNET (Top)` = read_tsv("data/clinvar/preds/unet_freq_features_top.tsv"),
  `PreGraph UNET (Top)` = read_tsv("data/clinvar/preds/unet_freq_structure_features_top.tsv"),
  `UNET (Finetune)` = read_tsv("data/clinvar/preds/unet_freq_finetune.tsv"),
  `PreGraph UNET (Finetune)` = read_tsv("data/clinvar/preds/unet_freq_structure_finetune.tsv"),
  .id = "model"
) %>%
  select(-wt) %>%
  pivot_wider(names_from = model, values_from = pred) %>% 
  {
    t <- read_tsv("data/clinvar/clinvar_test.tsv") %>% 
      mutate(pdb_id = str_to_upper(pdb_id)) %>%
      select(uniprot, position, wt, mut, pdb_id, chain, pdb_pos, clnsig, clnsig_patho, foldx_ddg, sift_score)
    left_join(t, ., by = c("pdb_id", "chain", "pdb_pos", "mut"))
  } %>%
  select(-pdb_id, -chain, -pdb_pos) %>%
  rename(SIFT4G = sift_score, FoldX = foldx_ddg) %>%
  mutate( FoldX = FoldX > 1, SIFT4G = SIFT4G < 0.05, `Baseline Frequency` = `Baseline Frequency` > 0.5,
          UNET = UNET > 0.5, `PreGraph UNET` = `PreGraph UNET` > 0.5, `Baseline ClinVar` = `Baseline ClinVar` > 0.5,
          `UNET (Top)` = `UNET (Top)` > 0.5, `PreGraph UNET (Top)` = `PreGraph UNET (Top)` > 0.5,
          `UNET (Finetune)` = `UNET (Finetune)` > 0.5, `PreGraph UNET (Finetune)` = `PreGraph UNET (Finetune)` > 0.5) %>%
  select(true = clnsig_patho, FoldX:`PreGraph UNET (Finetune)`) %>%
  pivot_longer(-true, names_to = "model", values_to = "pred") %>%
  drop_na() %>%
  group_by(model) %>%
  summarise(tp = sum(pred & true),
            tn = sum(!pred & !true),
            fp = sum(pred & !true),
            fn = sum(!pred & true)) %>%
  mutate(accuracy = (tp + tn) / (tp + tn + fp + fn),
         f1 = 2 * tp / (2 * tp + fp + fn),
         kappa = 2 * (tp * tn - fn * fp) / ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn)),
         model = factor(model, levels = rev(names(TOOL_COLOURS)))) %>%
  select(model, accuracy, f1, kappa) %>%
  pivot_longer(-model, names_to = "metric", values_to = "value")

p_performance <- ggplot(clinvar_acc, aes(x = model, y = value, fill = model)) +
  facet_wrap(~metric, scales = "free_x", strip.position = "bottom", nrow = 1,
             labeller = as_labeller(c(accuracy = "Accuracy", f1 = "F1 Score", kappa = "Cohen's &kappa;"))) +
  coord_flip() +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = TOOL_COLOURS) +
  labs(x = "", y = "") + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(2, "mm"),
        strip.placement = "outside",
        strip.text = element_markdown(margin = margin(0, 0, 0, 0, unit = "mm")))

### Panel - ROC curve ###
clinvar_roc <- read_tsv("data/clinvar/roc.tsv") %>%
  mutate(model_auc = auc_labeled_model(model, auc)) %>%
  filter(!str_detect(model, "Thresh"))

clinvar_auc <- distinct(clinvar_roc, model, model_auc, auc) %>%
  arrange(desc(auc)) %>%
  mutate(fpr = 1, tpr = rev(seq(from = 0.03, by = 0.05, length.out = n())))

p_clinvar <- ggplot(clinvar_roc, aes(x = fpr, y = tpr, colour = model, label = model_auc)) +
  geom_step(show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = clinvar_auc, hjust = 1, show.legend = FALSE, size = 2.5) +
  coord_fixed() +
  scale_colour_manual(values = TOOL_COLOURS) +
  labs(x = "False Positive Rate", y = "True Positive Rate")

### Panel - Other generalisation datasets ###
gen_roc <- bind_rows(DMS = read_tsv("data/dms/roc.tsv"), Jelier = read_tsv("data/jelier/roc.tsv"), .id = "type") %>%
  select(type, model, auc) %>%
  distinct() %>%
  mutate(model = factor(model, levels = reorder(model[type == "Jelier"], -auc[type == "Jelier"])))

p_generalisation <- ggplot(gen_roc, aes(x = model, y = auc, fill = model)) +
  facet_wrap(~type, ncol = 1) +
  geom_col(show.legend = FALSE, width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = TOOL_COLOURS) +
  labs(y = "AUC") +
  lims(y = c(0,1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        axis.title.y = element_blank())

### Figure Assembly ###
size <- theme(text = element_text(size = 11))
p1 <- p_preds + labs(tag = 'A') + size
p2 <- p_performance + labs(tag = 'B') + size
p3 <- p_clinvar + labs(tag = 'C') + size
p4 <- p_generalisation + labs(tag = 'D') + size

figure4 <- multi_panel_figure(width = c(90, 90), height = c(45, 55, 90), panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1:2) %>%
  fill_panel(p3, row = 3, column = 1) %>%
  fill_panel(p4, row = 3, column = 2)

ggsave('figures/figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4),
       units = 'mm', dpi = 600)
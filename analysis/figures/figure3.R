#!/usr/bin/env Rscript
# Produce figure 3 - Frequency Prediction
source('src/config.R')
source("src/analysis.R")

### Panel - Example Output ###
preds <- bind_rows(
  true = read_tsv('data/pssm/pn_casp12_validation.tsv') %>%
    extract(protein, into = c("pdb_id", "chain"), regex = "[0-9]*#([0-9A-Z]*)_[0-9]*_([A-Za-z0-9])") %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "pred") %>%
    drop_na(),
  pred = read_tsv('data/freq/unet_sequence_validation.tsv'),
  .id = "model"
) %>%
  filter(pdb_id == "2EK0") %>%
  pivot_wider(names_from = 'model', values_from = 'pred') %>%
  mutate(del = true < 0.01,
         correct = (pred > 0.5) == del)

p_preds <- ggplot(preds, aes(x = position, fill = pred)) +
  geom_tile(aes(y = mut), colour = "grey") +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "P(Del)", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1), midpoint = 0.5, low = "#2166ac", high = "#b2182b") +
  scale_x_continuous(expand = expansion(0), breaks = c(1, seq(10, 90, 10))) +
  guides(fill = guide_colourbar(barwidth = unit(5, "mm"), title.vjust = 1),
         colour = guide_legend(override.aes = list(fill = NA))) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right")

p_preds_box <- ggplot(preds, aes(x = pred > 0.44, y = true, fill = pred > 0.44)) + 
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(method = "wilcox", comparisons = list(c("TRUE", "FALSE"))) +
  labs(x = "Prediction", y = "True Frequency") +
  scale_x_discrete(labels = c(`FALSE` = expression(f >= 0.01), `TRUE` = expression(f < 0.01))) +
  scale_fill_manual(values = c(`FALSE` = "#2166ac", `TRUE` = "#b2182b"))

### Panel - Thresholds ###
training_logs <- read_tsv("data/freq/training_logs.tsv") %>%
  extract(model, into = c("experiment", "model", "dataset"), regex = "models/classifier/([^/]*)/([^/]*)/([^/]*)")

p_train_threshold <- filter(training_logs, experiment == "threshold", metric == "epoch_masked_accuracy", dataset == "validation") %>% 
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  mutate(model = str_remove(model, "t")) %>%
  ggplot(aes(x = model, y = value)) +
  geom_point(shape = 20, colour = "#e41a1c", size = 4) +
  scale_y_continuous(limits = c(0.75, 0.9), labels = function(x) str_remove(x, "0*^")) +
  labs(x = "Frequency Threshold", y = "Accuracy")

### Panel - Overall performance ###
classifier_performance <- read_tsv("data/freq/all_model_testing.tsv") %>%
  mutate(threshold = c(`UNET`=0.5, `PreGraph UNET`=0.5, `Baseline CNN`=0.5, `BLOSUM62`=-2, `SIFT4G`=0.05)[model],
         pred_del = ifelse(model %in% c("SIFT4G"), pred < threshold, pred > threshold)) %>%
  group_by(model, protein) %>%
  summarise(tp = sum(pred_del & deleterious),
            tn = sum(!pred_del & !deleterious),
            fp = sum(pred_del & !deleterious),
            fn = sum(!pred_del & deleterious),
            .groups = "drop") %>%
  mutate(accuracy = (tp + tn) / (tp + tn + fp + fn),
         recall = tp / (tp + fn),
         precision = tp / (tp + fp),
         f1 = 2 * tp / (2 * tp + fp + fn),
         kappa =  2 * (tp * tn - fn * fp) / ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn))) %>%
  select(model, protein, accuracy, recall, precision, f1, kappa) %>%
  pivot_longer(c(-model, -protein), names_to = "metric", values_to = "value") %>%
  mutate(model = factor(model, levels = c("BLOSUM62", "SIFT4G", "Baseline CNN", "UNET", "PreGraph UNET")),
         metric = factor(metric, levels = c("accuracy", "precision", "recall", "f1", "kappa")))

metric_labs <- c(accuracy="Accuracy", recall="Recall", precision="Precision", f1="F1 Score", kappa="Cohen's &kappa;")
p_performance <- ggplot(classifier_performance, aes(x = model, y = value, fill = model)) + 
  facet_wrap(~metric, nrow = 1, scales = "free_y", labeller = labeller(metric = metric_labs)) + 
  geom_boxplot(outlier.shape = 20, size = 0.25) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(name = "", values = TOOL_COLOURS) +
  scale_y_continuous(limits = function(x) if (any(x < 0)) c(-0.5, 1) else c(0, 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,-10,0),
        strip.text = element_markdown())

### Panel - ROC curve ###
classifier_roc <- read_tsv("data/freq/roc.tsv") %>%
  group_by(model) %>%
  mutate(pr_auc = pr_auc(tpr, precision)) %>%
  ungroup() %>%
  mutate(model_auc = auc_labeled_model(model, auc),
         model_pr_auc = auc_labeled_model(model, pr_auc))

classifier_auc <- distinct(classifier_roc, model, model_auc, auc) %>%
  arrange(desc(auc)) %>%
  mutate(fpr = 1, tpr = c(0.23, 0.18, 0.13, 0.08, 0.03),
         col = TOOL_COLOURS[model])

p_roc <- ggplot(classifier_roc, aes(x = fpr, y = tpr, colour = model_auc, label = model_auc)) +
  geom_step(show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = classifier_auc, hjust = 1, show.legend = FALSE, size = 3.5) +
  coord_fixed() +
  scale_colour_manual(name = "", values = structure(classifier_auc$col, names = as.character(classifier_auc$model_auc))) +
  labs(x = "False Positive Rate", y = "True Positive Rate")

### Panel - PR Curve ###
classifier_pr_auc <- distinct(classifier_roc, model, model_pr_auc, pr_auc) %>%
  arrange(desc(pr_auc)) %>%
  mutate(tpr = 1, precision = c(1, 0.95, 0.9, 0.85, 0.8),
         col = TOOL_COLOURS[model])

p_pr <- ggplot(classifier_roc, aes(x = tpr, y = precision, colour = model_pr_auc, label = model_pr_auc)) +
  geom_step(show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(data = classifier_pr_auc, hjust = 1, vjust = 1, show.legend = FALSE, size = 3.5) +
  coord_fixed() +
  scale_colour_manual(name = "", values = structure(classifier_pr_auc$col, names = as.character(classifier_pr_auc$model_pr_auc))) +
  labs(x = "Recall", y = "Precision")

### Figure Assembly ###
size <- theme(text = element_text(size = 12))
p1 <- p_rep_ubi + labs(tag = 'A') + sizet

figure3 <- multi_panel_figure(width = 150, height = 150, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1)

ggsave('figures/figures/figure3.pdf', figure3, width = figure_width(figure3), height = figure_height(figure3),
       units = 'mm', device = cairo_pdf())
ggsave('figures/figures/figure3.png', figure3, width = figure_width(figure3), height = figure_height(figure3),
       units = 'mm')
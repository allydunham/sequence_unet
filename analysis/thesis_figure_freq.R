#!/usr/bin/env Rscript
# Generate thesis figure on the basic model
source('src/config.R')
source("src/analysis.R")

tool_colours <- c(`PreGraph UNET`="#39f0f5",
                  UNET="#1f78b4",
                  `Baseline CNN`="#ff7f00",
                  BLOSUM62="#717171",
                  SIFT4G="#33a02c",
                  SPBuild="#f300ff")

#### Panel - Classifier Performance ####
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
  scale_fill_manual(name = "", values = tool_colours) +
  scale_y_continuous(limits = function(x) if (any(x < 0)) c(-0.5, 1) else c(0, 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,-10,0),
        strip.text = element_markdown())

#### Panel - Classifier ROC ####
classifier_roc <- read_tsv("data/freq/roc.tsv") %>%
  mutate(model_auc = auc_labeled_model(model, auc))

classifier_auc <- distinct(classifier_roc, model, model_auc, auc) %>%
  arrange(desc(auc)) %>%
  mutate(fpr = 1, tpr = c(0.23, 0.18, 0.13, 0.08, 0.03),
         col = tool_colours[model])

p_roc <- ggplot(classifier_roc, aes(x = fpr, y = tpr, colour = model_auc, label = model_auc)) +
  geom_step(show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = classifier_auc, hjust = 1, show.legend = FALSE, size = 3.5) +
  coord_fixed() +
  scale_colour_manual(name = "", values = structure(classifier_auc$col, names = as.character(classifier_auc$model_auc))) +
  labs(x = "False Positive Rate", y = "True Positive Rate")

#### Figure Assembly ####
size <- theme(text = element_text(size = 11))
p1 <- p_performance + labs(tag = 'A') + size
p2 <- p_roc + labs(tag = 'B') + size

figure <- multi_panel_figure(width = 150, height = c(50, 150), columns = 1,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 2, column = 1)

ggsave('figures/thesis_figure_freq.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_freq.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

#!/usr/bin/env Rscript
# Generate thesis figure on the basic model
source('src/config.R')
source("src/analysis.R")
library(png)

#### Panel - Model Structure ####
model_png <- readPNG("figures/model_schematic.png")
p_model <- ggplot() +
  geom_blank() +
  lims(x = c(0, 180), y = c(0, 80)) +
  coord_fixed() +
  annotation_raster(model_png, interpolate = TRUE, xmin = 0, xmax = 180, ymin = 0, ymax = 80) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())

#### Panel - Training ####
p_train <- blank_plot("Training Trace Over\nDifferent Factors")

#### Panel - Classifier Performance ####
classifier_roc <- read_tsv("data/freq/roc.tsv") %>%
  mutate(model_auc = auc_labeled_model(model, auc))

classifier_auc <- distinct(classifier_roc, model_auc, auc) %>%
  arrange(desc(auc)) %>%
  mutate(fpr = 1, tpr = c(0.2, 0.15, 0.1, 0.05))

p_classifier <- ggplot(classifier_roc, aes(x = fpr, y = tpr, colour = model_auc, label = model_auc)) +
  geom_step(show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = classifier_auc, hjust = 1, show.legend = FALSE, size = 2.5) +
  coord_fixed() +
  scale_colour_brewer(name = "", type = "qual", palette = "Dark2") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

#### Panel - PSSM Prediction ####
pssm_preds <- bind_rows(
  true = read_tsv('data/pssm/pn_casp12_validation.tsv') %>%
    extract(protein, into = c("pdb_id", "chain"), regex = "[0-9]*#([0-9A-Z]*)_[0-9]*_([A-Za-z0-9])") %>%
    pivot_longer(A:Y, names_to = "mut", values_to = "pred") %>%
    drop_na(),
  pred = read_tsv('data/pssm/unet_sequence_validation.tsv'),
  .id = "model"
) %>%
  filter(pdb_id == "1WAZ") %>%
  pivot_wider(names_from = 'model', values_from = 'pred') %>%
  mutate(diff = pred - true)

p_pssm_pred <- ggplot(pssm_preds, aes(x = position, fill = diff)) +
  geom_raster(aes(y = mut)) +
  geom_tile(aes(y = wt), fill = NA, colour = "black") +
  scale_fill_gradient2(name = "Pred - True", low = "#b2182b", high = "#2166ac") +
  scale_x_continuous(expand = expansion(0)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

#### Panel - PSSM Prediction Performance ####
pssm_models <- read_tsv("data/pssm/combined_preds.tsv")

pssm_cor <- group_by(pssm_models, model) %>%
  group_modify(~broom::tidy(cor.test(.$pred, .$true)))

p_pssm_cor <- ggplot(pssm_cor, aes(x = model, y = estimate, fill = model, ymin = conf.low, ymax = conf.high)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_errorbar(width = 0.6) +
  coord_flip() +
  scale_fill_brewer(name = 'Model', type = 'qual', palette = 'Dark2') +
  labs(y = "Pearson Correlation Coefficient") +
  theme(axis.title.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank())

pssm_summary <- drop_na(pssm_models) %>%
  mutate(abs_diff = abs(diff)) %>%
  group_by(pdb_id, chain, position, wt, mut) %>%
  mutate(best_diff = min(abs_diff)) %>%
  ungroup() %>%
  mutate(good_model = abs_diff < 2,
         best_model = diff == best_diff) %>%
  group_by(model) %>%
  summarise(n_tot = n(),
            n_best = sum(best_model) / n_tot,
            n_good = sum(good_model) / n_tot,
            n_best_good = sum(best_model & good_model) / n_tot,
            .groups = 'drop') %>%
  select(-n_tot) %>%
  pivot_longer(-model, names_prefix = 'n_', names_to = 'type', values_to = 'prop') %>%
  mutate(type = factor(type, levels = c('good', 'best', 'best_good')))

p_pssm_summary <- ggplot(pssm_summary, aes(x = type, y = prop, fill = model)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  scale_fill_brewer(name = '', type = 'qual', palette = 'Dark2', guide = guide_legend(reverse = TRUE)) +
  labs(x = '', y = 'Proportion of Predictions') +
  scale_x_discrete(labels = c(good='|Pred - True| <= 1', best='Best Model', best_good='Both')) +
  theme(axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        legend.key.size = unit(2, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"))

#### Figure Assembly ####
size <- theme(text = element_text(size = 9))
p1 <- p_model + labs(tag = 'A') + size
p2 <- p_train + labs(tag = 'B') + size
p3 <- p_classifier + labs(tag = 'C') + size
p4 <- p_pssm_pred + labs(tag = 'D') + size
p5 <- p_pssm_cor + labs(tag = 'E') + size
p6 <- p_pssm_summary + labs(tag = 'F') + size

figure <- multi_panel_figure(width = 180, height = 240, columns = 6, rows = 6,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1:2, column = 1:6) %>%
  fill_panel(p2, row = 3:4, column = 1:3) %>%
  fill_panel(p3, row = 3:4, column = 4:6) %>%
  fill_panel(p4, row = 5:6, column = 1:3) %>%
  fill_panel(p5, row = 5, column = 4:6) %>%
  fill_panel(p6, row = 6, column = 4:6)

ggsave('figures/thesis_figure_model.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_model.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

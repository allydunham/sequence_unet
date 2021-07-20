#!/usr/bin/env Rscript
# Generate thesis figure on the basic model
source('src/config.R')
source("src/analysis.R")

#### Panel - Training ####
training_logs <- read_tsv("data/freq/training_logs.tsv") %>%
  extract(model, into = c("experiment", "model", "dataset"), regex = "models/classifier/([^/]*)/([^/]*)/([^/]*)")

# Calculate variance over replicates for text
# filter(training_logs, experiment == "replicates") %>% 
#   filter(step == max(step)) %>% 
#   group_by(dataset, metric) %>% 
#   summarise(mean = mean(value), var = var(value))

p_train_threshold <- filter(training_logs, experiment == "threshold", metric == "epoch_masked_accuracy", dataset == "validation") %>% 
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  mutate(model = str_remove(model, "t")) %>%
  ggplot(aes(x = model, y = value)) +
  geom_point(shape = 20, colour = "#e41a1c") +
  scale_y_continuous(limits = c(0.75, 0.9), labels = function(x) str_remove(x, "0*^")) +
  labs(x = "Frequency Threshold", y = "Accuracy")

p_train_activation <- filter(training_logs, experiment == "activation", metric == "epoch_masked_accuracy", dataset == "validation") %>%
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(model, value), y = value)) +
  geom_point(shape = 20, colour = "#377eb8") +
  coord_flip(clip = "off") +
  scale_x_discrete(labels = c(elu = "ELU", hard_sigmoid = "Hard Sigmoid", relu = "ReLU", swish = "Swish", tanh = "Tanh")) +
  scale_y_continuous(limits = c(0.73, 0.75), labels = function(x) str_remove(x, "0*^")) +
  labs(x = "", y = "Accuracy") +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank())

training_log_size <- filter(training_logs, experiment == "size", dataset == "validation") %>%
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  extract(model, c("filters", "kernel_size", "layers"), "f([0-9]*)_k([0-9]*)_l([0-9]*)", convert = TRUE)

p_train_filters <- filter(training_log_size, layers == 6, kernel_size == 9, metric == "epoch_masked_accuracy") %>%
  ggplot(aes(x = filters, y = value)) +
  geom_point(shape = 20, colour = "#4daf4a") +
  coord_cartesian(clip = "off") +labs(x = "Base CNN Filters", y = "Accuracy")

p_train_kernel <- filter(training_log_size, layers == 6, filters == 32, metric == "epoch_masked_accuracy") %>%
  ggplot(aes(x = kernel_size, y = value)) +
  geom_point(shape = 20, colour = "#984ea3") +
  scale_x_continuous(breaks = c(3,5,7,9)) +
  coord_cartesian(clip = "off") +labs(x = "Kernel Size", y = "Accuracy")

p_train_layers <- filter(training_log_size, filters == 32, kernel_size == 5, metric == "epoch_masked_accuracy") %>%
  ggplot(aes(x = layers, y = value)) +
  geom_point(shape = 20, colour = "#ff7f00") +
  scale_x_continuous(breaks = c(2,4,6)) +
  coord_cartesian(clip = "off") +labs(x = "Layers", y = "Accuracy")

p_train_structure <- filter(training_logs, experiment == "structure", metric == "epoch_masked_accuracy", dataset == "validation") %>%
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  separate(model, c("activation", "layers"), extra = "merge") %>%
  mutate(nlayers = ifelse(layers == "none", 0, str_count(layers, "_") + 1)) %>%
  filter(activation != "relu", str_starts(layers, "none|32")) %>%
  ggplot(aes(x = nlayers, y = value)) +
  geom_point(shape = 20, colour = "#a65628") +
  coord_cartesian(clip = "off") +
  labs(x = "Graph Layers", y = "Accuracy")

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
  guides(fill = guide_colourbar(barheight = unit(3, "mm"), title.vjust = 1)) +
  labs(x = "Position", y = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

#### Panel - PSSM Prediction Performance ####
pssm_models <- read_tsv("data/pssm/combined_preds.tsv") %>%
  mutate(model = factor(model, levels = c("BLOSUM62", "SPBuild", "UNET", "PreGraph UNET")))

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
p1 <- p_train_threshold + labs(tag = 'A') + size
p2 <- p_train_activation + labs(tag = 'B') + size
p3 <- p_train_filters + labs(tag = 'C') + size
p4 <- p_train_kernel + labs(tag = 'D') + size
p5 <- p_train_layers + labs(tag = 'E') + size
p6 <- p_train_structure + labs(tag = 'F') + size

p7 <- p_classifier + labs(tag = 'G') + size
p8 <- p_pssm_pred + labs(tag = 'H') + size
p9 <- p_pssm_cor + labs(tag = 'I') + size
p10 <- p_pssm_summary + labs(tag = 'J') + size

figure <- multi_panel_figure(width = c(30, 30, 30, 30, 30, 30), height = c(40, 40, 80, 40),
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 1, column = 3:4) %>%
  fill_panel(p3, row = 1, column = 5:6) %>%
  fill_panel(p4, row = 2, column = 1:2) %>%
  fill_panel(p5, row = 2, column = 3:4) %>%
  fill_panel(p6, row = 2, column = 5:6) %>%
  fill_panel(p7, row = 3, column = 1:3) %>%
  fill_panel(p8, row = 3, column = 4:6) %>%
  fill_panel(p9, row = 4, column = 1:3) %>%
  fill_panel(p10, row = 4, column = 4:6)

ggsave('figures/thesis_figure_model.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_model.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

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
  geom_point(shape = 20, colour = "#e41a1c", size = 4) +
  scale_y_continuous(limits = c(0.75, 0.9), labels = function(x) str_remove(x, "0*^")) +
  labs(x = "Frequency Threshold", y = "Accuracy")

p_train_activation <- filter(training_logs, experiment == "activation", metric == "epoch_masked_accuracy", dataset == "validation") %>%
  group_by(model) %>%
  filter(step == step[which.max(value)]) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(model, value), y = value)) +
  geom_point(shape = 20, colour = "#377eb8", size = 4) +
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
  geom_point(shape = 20, colour = "#4daf4a", size = 4) +
  scale_x_continuous(breaks = c(32, 48, 64, 80, 96)) +
  ylim(0.75, 0.76) +
  coord_cartesian(clip = "off") +labs(x = "Base CNN Filters", y = "Accuracy")

p_train_kernel <- filter(training_log_size, layers == 6, filters == 32, metric == "epoch_masked_accuracy") %>%
  ggplot(aes(x = kernel_size, y = value)) +
  geom_point(shape = 20, colour = "#984ea3", size = 4) +
  scale_x_continuous(breaks = c(3,5,7,9)) +
  coord_cartesian(clip = "off") +labs(x = "Kernel Size", y = "Accuracy")

p_train_layers <- filter(training_log_size, filters == 32, kernel_size == 5, metric == "epoch_masked_accuracy") %>%
  ggplot(aes(x = layers, y = value)) +
  geom_point(shape = 20, colour = "#ff7f00", size = 4) +
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
  geom_point(shape = 20, colour = "#a65628", size = 4) +
  coord_cartesian(clip = "off", ylim = c(0.75, 0.775)) +
  labs(x = "Graph Layers", y = "Accuracy")

#### Figure Assembly ####
size <- theme(text = element_text(size = 11))
p1 <- p_train_threshold + labs(tag = 'A') + size
p2 <- p_train_activation + labs(tag = 'B') + size
p3 <- p_train_filters + labs(tag = 'C') + size
p4 <- p_train_kernel + labs(tag = 'D') + size
p5 <- p_train_layers + labs(tag = 'E') + size
p6 <- p_train_structure + labs(tag = 'F') + size

figure <- multi_panel_figure(width = c(90, 90), height = c(50, 50, 50),
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2) %>%
  fill_panel(p5, row = 3, column = 1) %>%
  fill_panel(p6, row = 3, column = 2)

ggsave('figures/thesis_figures/thesis_figure_training.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm',
       device = cairo_pdf)
ggsave('figures/thesis_figures/thesis_figure_training.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

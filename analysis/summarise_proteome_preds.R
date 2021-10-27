#!/usr/bin/env Rscript
# Summarise wide Sequence UNET prediction files
library(tidyverse)

input_files <- commandArgs(TRUE)

process_gene <- function(x, ...) {
  long <- pivot_longer(x, A:Y, names_to = "mut", values_to = "pred")
  
  mat <- group_by(long, wt, mut) %>%
    summarise(m = mean(pred), .groups = "drop") %>%
    mutate(pair = str_c(wt, mut, sep = "_")) %>%
    select(pair, m) %>%
    pivot_wider(names_from = pair, values_from = m)
  
  tibble(length = max(x$position),
         mean_pred = mean(long$pred),
         mean_wt = mean(long$pred[long$mut == long$wt]),
         mean_mut = mean(long$pred[!long$mut == long$wt])) %>%
    bind_cols(mat)
}

process_file <- function(x, print_col_names=FALSE) {
  read_tsv(x) %>%
    group_by(gene) %>%
    group_modify(process_gene) %>%
    format_tsv(col_names = print_col_names) %>%
    cat()
}

walk(input_files, ~process_file(., print_col_names = . == input_files[1]))

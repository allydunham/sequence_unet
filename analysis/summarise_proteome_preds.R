#!/usr/bin/env Rscript
# Summarise wide Sequence UNET prediction files
library(tidyverse)

input_files <- commandArgs(TRUE)

process_file <- function(x) {
  df <- read_tsv(x)
  group_by(df, gene) %>% summarise(across(A:Y, mean)) %>% mutate(overall = rowMeans(select(., A:Y)))
}

map_df(input_files, process_file) %>%
  bind_rows() %>%
  format_tsv() %>%
  cat()

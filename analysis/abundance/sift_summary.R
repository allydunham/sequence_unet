#!/usr/bin/env Rscript
# Summarise SIFT4G results
library(tidyverse)

bind_rows(
  `Homo sapiens`=read_tsv("/nfs/research1/beltrao/ally/mutations/data/mutfunc/human/conservation/sift_parsed_all.tab", comment = "#"),
  `Saccharomyces cerevisiae`=read_tsv("/nfs/research1/beltrao/ally/mutations/data/mutfunc/yeast/conservation/sift_parsed_all.tab", comment = "#"),
  `Escherichia coli`=read_tsv("/nfs/research1/beltrao/ally/mutations/data/mutfunc/ecoli/conservation/sift.tab", comment = "#"),
  .id = "organism"
) %>%
  group_by(acc, pos) %>%
  summarise(mean_score = mean(score[ref != alt]),
            n_conserved = sum(score[ref != alt] < 0.05),
            .groups = "drop_last") %>%
  summarise(mean_sift = mean(mean_score),
            percent_avg_conserved = sum(mean_score < 0.05) / n(),
            mean_conserved = mean(n_conserved),
            percent_n_conserved = sum(n_conserved > 9) / n()) %>%
  write_tsv("data/abundance/mutfunc_sift_summary.tsv")
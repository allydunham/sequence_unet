#!/usr/bin/env Rscript
# DMS Data Split
library(tidyverse)

dms <- read_tsv("data/dms/long_combined_mutational_scans.tsv") %>%
  select(study, gene, position, wt, mut, score=imputed_score) %>%
  filter(!gene == "UBI" | study == "roscoe_2014_ubi") %>% # Take the later UBI study
  filter(!(gene == "HSP90" & study == "hietpas_2011_hsp90")) %>% # Take the later of Hietpas 2011 & Jiang 2013 which cover the same part of hsp90
  filter(!mut == "*") %>% # Ignore nonsense
  pivot_wider(names_from = mut, values_from = score) %>%
  mutate(gene = str_to_lower(str_remove_all(gene, " ")))

test_study <- filter(dms, study == "matreyek_2018_pten")
dms <- filter(dms, !study == "matreyek_2018_pten")

row_train <- sample(seq_len(nrow(dms)), 5500)
dms_train <- dms[row_train,]
dms_test <- dms[-row_train,]
  
row_test <- sample(seq_len(nrow(dms_test)), 400)
dms_val <- dms_test[-row_test,]
dms_test <- bind_rows(test_study, dms_test[row_test,])

write_tsv(dms_train, "data/dms/train.tsv")
write_tsv(dms_val, "data/dms/val.tsv")
write_tsv(dms_test, "data/dms/test.tsv")

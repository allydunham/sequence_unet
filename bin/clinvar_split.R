#!/usr/bin/env Rscript
# ClinVar Data Split
library(tidyverse)
library(caret)

clinvar <- read_tsv("data/clinvar/clinvar_variants.tsv") %>%
  mutate(clnsig_patho = clnsig %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")) %>%
  distinct(uniprot, position, wt, mut, .keep_all = TRUE)

split <- createDataPartition(clinvar$clnsig_patho, p = 0.95)

clinvar_train <- clinvar[split$Resample1,]
clinvar_test <- clinvar[-split$Resample1,]

# Split val / test
val_split <- createDataPartition(clinvar_test$clnsig_patho, p = 0.9)
clinvar_val <- clinvar_test[-val_split$Resample1,]
clinvar_test <- clinvar_test[val_split$Resample1,]

# Write tables
write_tsv(clinvar_train, "data/clinvar/clinvar_train.tsv")
write_tsv(clinvar_test, "data/clinvar/clinvar_test.tsv")
write_tsv(clinvar_val, "data/clinvar/clinvar_val.tsv")

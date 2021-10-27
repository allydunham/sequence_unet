#!/usr/bin/env Rscript
# Summarise wide Sequence UNET prediction files
library(data.table)
library(purrr)

amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

exp_pairs <- gtools::permutations(20, 2, v = amino_acids, repeats.allowed = TRUE)
exp_pairs <- paste(exp_pairs[,1], exp_pairs[,2], sep = "_")

exp_mat <- data.table(pair = exp_pairs)

input_files <- commandArgs(TRUE)

process_gene <- function(x, ...) {
  long <- melt(x, measure.vars = amino_acids, variable.name = "mut", value.name = "pred")
  mat <- long[, .(pair = paste0(wt, "_", mut), mean = mean(pred)), keyby = .(wt, mut)][,.(pair, mean)]
  mat <- merge(mat, exp_mat, all = TRUE)   
  
  cbind(data.table(length = nrow(x),
                   mean_pred = mean(long$pred),
                   mean_wt = mean(long$pred[long$mut == long$wt]),
                   mean_mut = mean(long$pred[!long$mut == long$wt])),
         dcast(mat, formula = . ~ pair, value.var = "mean", drop = FALSE)[,-1])
}

process_file <- function(x, print_col_names=FALSE) {
  dt <- fread(x, sep = "\t")
  s <- dt[, process_gene(.SD), by = gene]
  fwrite(s, quote = FALSE, sep = "\t", col.names = print_col_names)
}

walk(input_files, ~process_file(., print_col_names = . == input_files[1]))

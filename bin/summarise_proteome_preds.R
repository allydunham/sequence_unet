#!/usr/bin/env Rscript
# Summarise wide Sequence UNET prediction files
library(argparser)
library(data.table)
library(purrr)

amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

parser <- arg_parser("Summarise wide Sequence UNET prediction files")
parser <- add_argument(parser, "files", help="Files to summarise", type="character")
parser <- add_argument(parser, "--mut_thresh", help="Number of deleterious variants required for a position to be conserved", default=9)
parser <- add_argument(parser, "--del_thresh", help="Deleteriousness threshold", default=0.5)
args <- parse_args(parser)

input_files <- commandArgs(TRUE)

# Guess thresholds from file name, defaulting

process_gene <- function(x, ...) {
  long <- melt(x, measure.vars = amino_acids, variable.name = "mut", value.name = "pred")
  pos_summary <- long[!wt == mut, .(mean_pred = mean(pred), n_conserved = sum(pred > args$del_thresh)), keyby = position]
  
  data.table(mean_pred = mean(long$pred),
             mean_wt = mean(long$pred[long$mut == long$wt]),
             mean_mut = mean(long$pred[!long$mut == long$wt]),
             percent_avg_conserved = sum(pos_summary$mean_pred > args$del_thresh) / nrow(x),
             percent_n_conserved = sum(pos_summary$n_conserved > args$mut_thresh) / nrow(x))
}

process_file <- function(x, print_col_names=FALSE) {
  dt <- fread(x, sep = "\t")
  s <- dt[, process_gene(.SD), by = gene]
  fwrite(s, quote = FALSE, sep = "\t", col.names = print_col_names)
}

walk(args$files, ~process_file(., print_col_names = . == input_files[1]))

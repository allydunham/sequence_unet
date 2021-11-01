#!/usr/bin/env Rscript
# Summarise wide Sequence UNET prediction files
library(argparser)
library(data.table)
library(purrr)

amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

parser <- arg_parser("Summarise wide Sequence UNET prediction files")
parser <- add_argument(parser, "--files", help="Files to summarise", nargs = Inf, type="character")
parser <- add_argument(parser, "--mut_thresh", help="Number of deleterious variants required for a position to be conserved", default=9)
parser <- add_argument(parser, "--del_thresh", help="Deleteriousness threshold", default=0.5)
parser <- add_argument(parser, "--del_less", short="-l", help="Deleteriousness threshold", flag=TRUE)
args <- parse_args(parser)

if (args$del_less) {
  mut_del <- function(x) x < args$del_thresh
} else {
  mut_del <- function(x) x > args$del_thresh
}


process_gene <- function(x, ...) {
  long <- melt(x, measure.vars = amino_acids, variable.name = "mut", value.name = "pred")
  pos_summary <- long[!wt == mut, .(mean_pred = mean(pred), n_conserved = sum(mut_del(pred))), keyby = position]
  
  data.table(mean_pred = mean(long$pred),
             mean_wt = mean(long$pred[long$mut == long$wt]),
             mean_mut = mean(long$pred[!long$mut == long$wt]),
             percent_avg_conserved = sum(mut_del(pos_summary$mean_pred)) / nrow(x),
             percent_n_conserved = sum(pos_summary$n_conserved > args$mut_thresh) / nrow(x))
}

process_file <- function(x, print_col_names=FALSE) {
  dt <- fread(x, sep = "\t")
  s <- dt[, process_gene(.SD), by = gene]
  fwrite(s, quote = FALSE, sep = "\t", col.names = print_col_names)
}

walk(args$files, ~process_file(., print_col_names = . == input_files[1]))

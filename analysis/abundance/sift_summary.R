#!/usr/bin/env Rscript
# Summarise SIFT4G results
library(data.table)

sift <- rbind(
  `Homo sapiens`=fread("/nfs/research1/beltrao/ally/mutations/data/mutfunc/human/conservation/sift_parsed_all.tab", sep="\t", skip="acc"),
  `Saccharomyces cerevisiae`=fread("/nfs/research1/beltrao/ally/mutations/data/mutfunc/yeast/conservation/sift_parsed_all.tab", sep="\t", skip="acc"),
  `Escherichia coli`=fread("/nfs/research1/beltrao/ally/mutations/data/mutfunc/ecoli/conservation/sift.tab", sep="\t", skip="acc"),
  .idcol = "organism"
)

pos_summary <- sift[ref != alt, .(mean_score = mean(score), n_conserved = sum(score < 0.05)), by = .(acc, pos)]
overall_summary <- pos_summary[, .(mean_sift = mean(mean_score),
                                   percent_avg_conserved = sum(mean_score < 0.05) / .N,
                                   mean_conserved = mean(n_conserved),
                                   percent_n_conserved = sum(n_conserved > 9) / .N), acc]
fwrite(overall_summary, "data/abundance/mutfunc_sift_summary.tsv", quote = FALSE, sep = "\t", col.names = print_col_names)

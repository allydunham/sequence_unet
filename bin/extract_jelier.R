#!/usr/bin/env Rscript
# Format Jelier data for analysis
source('src/analysis/config.R')

variants <- readxl::read_xls("data/jelier/jelier_variants.xls", skip = 3) %>%
  select(gene = `gene ID`, variation = Variation, effect = Effect) %>%
  extract(variation, c("wt", "position", "mut"), "([A-Z])([0-9]*)([A-Z])") %>%
  select(gene, position, wt, mut, effect) %>%
  mutate(effect = ifelse(effect == 0, "neutral", "deleterious"))
write_tsv(variants, "data/jelier/jelier_variants.tsv")

seqs <- Biostrings::readAAStringSet("data/jelier/cerevisiae.fasta")
seqid <- str_split_fixed(names(seqs), " ", 2)[,1]
names(seqs) <- seqid
seqs <- seqs[seqid %in% variants$gene]
seqs <- Biostrings::subseq(seqs, end = -2)
Biostrings::writeXStringSet(seqs, "data/jelier/variant_genes.fa")


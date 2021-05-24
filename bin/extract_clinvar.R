#!/usr/bin/env Rscript
# Extract gene SIFT/FoldX/PDB IDs from ClinVar VCF
library(tidyverse)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg38)
library(siftr)

id_mapping <- read_tsv('meta/HUMAN_9606_idmapping_selected.tab',
                       col_names = c('uniprot', 'uniprot_id', 'entrez_id', 'refseq', 'gi',
                                     'pdb', 'go', 'uniref100', 'uniref90', 'uniref50', 'uniparc', 'pir', 'ncbi_taxon',
                                     'mim', 'unigene', 'pubmed', 'embl', 'embl_cds', 'ensembl', 'ensembl_trs', 'ensembl_pro', 'other_pubmed'),
                       col_types = cols(.default = col_character())) %>%
  dplyr::select(entrez_id, uniprot) %>%
  distinct()

# Read ClinVar VCF
chroms <- c(1:22, 'X', 'Y', 'MT')
chr_filter <- GRanges(seqnames = chroms,
                      ranges = IRanges(start = rep(0, length(chroms)),
                                       end = rep(536870912, length(chroms)),))

tbi <- TabixFile('data/clinvar/clinvar_20210424.vcf.gz')
clinvar <- readVcf(tbi, "hg38", param = chr_filter)

# Change chr names to match hg19 packages
seqlevels(rowRanges(clinvar)) <- structure(str_c('chr', str_replace(chroms, 'MT', 'M')), names=chroms)

# Identify Coding Variants
clinvar_coding <- predictCoding(clinvar, TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg38)
clinvar_coding <- clinvar_coding[(sapply(clinvar_coding$PROTEINLOC, length) == 1) & (clinvar_coding$CONSEQUENCE == 'nonsynonymous')]

# Generate clean table of variants
clnsig <- tibble(chr = as.vector(seqnames(rowRanges(clinvar))),
                 clnsig = as.matrix(info(clinvar)$CLNSIG)[,1]) %>%
  bind_cols(as_tibble(ranges(rowRanges(clinvar)))) %>%
  dplyr::select(names, chr, genomic_start = start, genomic_end = end, width, clnsig)

protein_variants <- bind_cols(as_tibble(ranges(clinvar_coding)),
          tibble(chr = as.vector(seqnames(clinvar_coding)),
                 query = clinvar_coding$QUERYID,
                 entrez_id = clinvar_coding$GENEID,
                 position = drop(clinvar_coding$PROTEINLOC),
                 wt = as.character(clinvar_coding$REFAA),
                 mut = as.character(clinvar_coding$VARAA))) %>%
  dplyr::select(names, chr, genomic_start = start, genomic_end = end, width, query, entrez_id, position, wt, mut) %>%
  left_join(clnsig, by = c('names', 'chr', 'genomic_start', 'genomic_end', 'width')) %>%
  distinct() %>%
  filter(clnsig %in% c('Benign', 'Likely_benign', 'Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic'))

# Add Sift and FoldX scores
foldx <- read_tsv('data/foldx.tsv.gz')

id_mapping <- filter(id_mapping, uniprot %in% unique(foldx$uniprot_id)) %>%
  distinct(entrez_id, .keep_all = TRUE)

protein_variants <- dplyr::select(protein_variants, entrez_id, position, wt, mut, clnsig) %>%
  left_join(id_mapping, by = 'entrez_id') %>%
  drop_na()

protein_variants <- filter(foldx, uniprot_id %in% unique(protein_variants$uniprot)) %>%
  group_by(pdb_id) %>%
  mutate(first_pdb_pos = min(pdb_pos)) %>%
  dplyr::select(uniprot=uniprot_id, position=uniprot_pos, wt=aa_wt, mut=aa_mt, foldx_ddg=ddG, pdb_id, chain, pdb_pos, first_pdb_pos) %>%
  left_join(protein_variants, ., by = c('uniprot', 'position', 'wt', 'mut')) %>%
  drop_na()

load('data/sift_mat_all.Rdata')
ids <- unique(protein_variants$uniprot)
sift <- map(sift_mat_all[ids[ids %in% names(sift_mat_all)]], ~as_tibble(filterPredictions(., score_thresh = 1, ic_thresh = 100, residue_thresh = 2))) %>%
  bind_rows(.id = 'uniprot')

protein_variants <- left_join(protein_variants, dplyr::select(sift, uniprot, position=pos, wt=ref, mut=alt, sift_score=score),
                              by = c('uniprot', 'position', 'wt', 'mut')) %>%
  drop_na()

# Write Table
write_tsv(protein_variants, 'data/clinvar/clinvar_variants.tsv')


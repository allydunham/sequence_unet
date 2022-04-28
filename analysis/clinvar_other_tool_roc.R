#!/usr/bin/env Rscript
# dbNDSFP and EVE tool clinar roc
source('src/config.R')
source("src/analysis.R")
data("BLOSUM62", package = "Biostrings")
plots <- list()

# Load ClinVar (inc FoldX, SIFT4G and BLOSUM62)
blosum <- as_tibble(BLOSUM62, rownames = 'wt') %>%
  pivot_longer(-wt, names_to = 'mut', values_to = 'BLOSUM62')

clinvar_stats <- read_tsv("data/clinvar/clinvar_test.tsv") %>% 
  select(uniprot, pdb_id, chain, pdb_pos, position, wt, mut, clnsig, clnsig_patho, FoldX = foldx_ddg, SIFT4G = sift_score) %>%
  left_join(blosum, by = c("wt", "mut"))

# Load EVE results
# Requires a download of the EVE csv results from evemodel.org placed in data/eve_variants
uniprot_mapping <- read_tsv("data/clinvar/uniprot_name_mapping.tsv") %>% # Only Uniprot IDs in ClinVar included
  select(uniprot = Entry, uniprot_name = `Entry name`)

load_eve <- function(uniprot_name) {
  path <- str_c("data/eve_variants/", uniprot_name, ".csv")
  
  if (!file.exists(path)) {
    warning(str_c("Gene not found: ", uniprot_name))
    return()
  }
  
  eve_cols <- cols(.default = col_character(), position = col_integer(), EVE_scores_ASM = col_double())
  
  read_csv(path, col_types = eve_cols) %>%
    mutate(uniprot_name = uniprot_name) %>%
    left_join(uniprot_mapping, by = "uniprot_name") %>%
    select(uniprot, uniprot_name, position, wt = wt_aa, mut = mt_aa, `EVE` = EVE_scores_ASM) %>%
    return()
}

eve_scores <- unique(uniprot_mapping$uniprot_name) %>%
  map(~suppressMessages(load_eve(.))) %>%
  bind_rows() %>%
  drop_na()

# Import ESM-1v results
esm_scores <- read_tsv("data/esm/1v_clinvar_preds.tsv") %>%
  pivot_longer(starts_with("pred_"), names_to = "mut", names_prefix = "pred_", values_to = "s") %>%
  pivot_wider(names_from = "model", values_from = "s") %>%
  mutate(`ESM-1v` = rowMeans(across(c(`1`, `2`, `3`, `4`, `5`)))) %>%
  select(id, position, wt, mut, `ESM-1v`) %>%
  separate(id, c("pdb_id", "pdb_model", "chain"))

# Generate dbNSFP query
ensembl_cannonical <- read_tsv("data/clinvar/human_grch38_105_ensembl_cannonical.tsv") %>%
  pull(transcript_stable_id) %>%
  str_split("\\.", simplify = TRUE) %>%
  magrittr::extract(1:nrow(.), 1) %>%
  unique()

dbnsfp_query <- read_tsv("data/clinvar/uniprot_name_mapping.tsv") %>%
  select(uniprot = Entry, ensembl = `Ensembl transcript`) %>%
  left_join(read_tsv("data/clinvar/clinvar_variants.tsv"), ., by = "uniprot") %>%
  separate(ensembl, into = str_c("ensembl", 1:41), sep = ";", fill = "right") %>%
  pivot_longer(ensembl1:ensembl41, values_to = "ensembl") %>%
  select(-name) %>%
  filter(!is.na(ensembl), !ensembl == "", !str_count(ensembl, "\\[") > 0, ensembl %in% ensembl_cannonical) %>%
  distinct() %>%
  mutate(query = str_c("Ensembl:", ensembl, ":", wt, position, mut))
# write_lines(dbnsfp_query$query, "data/clinvar/dbnsfp_ensebml_query.txt")

# Load dbNSFP results
# Requires a local search of dbNSFP using the query generated above:
# java -Xmx10g search_dbNSFP43a -i query.txt -o results.out -w 1-18,29-31,38-172,631-639
# The results must be placed in data/clinvar/dbnsfp_clinvar.tsv
get_cannonical_index <- function(x) {
  x <- str_split(x, ";", simplify = TRUE)
  cannonical <- matrix(x %in% ensembl_cannonical, nrow = nrow(x), ncol = ncol(x))
  inds <- as_tibble(which(cannonical, arr.ind = TRUE)) %>%
    distinct(row, .keep_all = TRUE) %>% # Only use first case for each row
    arrange(row)
  return(inds$col)
}

select_cannonical <- function(x, ind) {
  if (!is.character(x)) {
    return(x)
  }
  
  if (!any(str_detect(x, ";"))) {
    return(x)
  }
  
  inds <- matrix(c(1:length(x), ind), byrow = FALSE, nrow = length(x), ncol = 2)
  out <- str_split(x, ";", simplify = TRUE)[inds]
  out[out == "."] <- NA
  return(out)
}

dbnsfp <- read_tsv("data/clinvar/dbnsfp_clinvar.tsv", na = c("NA", "na", ".", "-")) %>%
  select(wt = aaref, mut = aaalt, position = aapos, gene_name = genename, ensembl_transcript = Ensembl_transcriptid, uniprot = Uniprot_acc,
         uniprot_entry = Uniprot_entry, PolyPhen2 = Polyphen2_HVAR_score, MutationAssessor = MutationAssessor_score,
         FATHMM = FATHMM_score, PROVEAN = PROVEAN_score, VEST4 = VEST4_score, `M-CAP` = `M-CAP_score`, REVEL = REVEL_score,
         MutPred = MutPred_score, MVP = MVP_score, MPC = MPC_score, PrimateAI = PrimateAI_score, DEOGEN2 = DEOGEN2_score,
         ClinPred = ClinPred_score, `LIST S2` = `LIST-S2_score`, CADD = CADD_raw, DANN = DANN_score, `FATHMM-MKL` = `fathmm-MKL_coding_score`,
         `FATHMM-XF` = `fathmm-XF_coding_score`, GenoCanyon = GenoCanyon_score, `GERP++` = `GERP++_RS`, clinvar_id, clinvar_clnsig) %>%
  mutate(cannonical = get_cannonical_index(ensembl_transcript),
         across(position:`GERP++`, ~select_cannonical(., cannonical)),
         across(c(position, PolyPhen2:`GERP++`), as.numeric))

# Combine data
preds <- left_join(clinvar_stats, select(eve_scores, uniprot, position, wt, mut, EVE), by = c("uniprot", "position", "wt", "mut")) %>%
  left_join(select(dbnsfp, uniprot, position, wt, mut, PolyPhen2:`GERP++`), by = c("uniprot", "position", "wt", "mut")) %>%
  left_join(mutate(select(esm_scores, -pdb_model), pdb_id = str_to_lower(pdb_id)), by = c("pdb_id", "chain", "pdb_pos"="position", "wt", "mut"))

# Calculate ROC
less = c("SIFT4G", "BLOSUM62", "FATHMM", "PROVEAN", "ESM-1v")
roc <- pivot_longer(preds, c(-pdb_id, -chain, -pdb_pos, -uniprot, -position, -wt, -mut, -clnsig, -clnsig_patho),
                    names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  group_modify(~calc_roc(.x, clnsig_patho, pred, greater = !(.y %in% less), max_steps = 5000)) %>%
  mutate(pr_auc = pr_auc(tpr, precision)) %>%
  ungroup() %>%
  arrange(model)
write_tsv(roc, "data/clinvar/other_tools_roc.tsv")

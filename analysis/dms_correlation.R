#!/usr/bin/env Rscript
# Generate correlations between tools and DMS data
source('src/config.R')
source("src/analysis.R")
plots <- list()

# Import DMS data
# Uses the data from Dunham and Beltrao 2021 (MSB)
dms <- read_tsv("data/dms/long_combined_mutational_scans.tsv") %>%
  left_join(select(read_tsv("data/dms/gene_summary.tsv"), gene = Gene, uniprot = `Uniprot ID`, species = Species), by = "gene") %>%
  filter(class == "Missense") %>%
  select(gene, uniprot, species, position, wt, mut, score, SIFT4G=sift, FoldX = total_energy)

# Import Sequence UNET data
unet <- bind_rows(
  `Baseline ClinVar` = read_tsv("data/dms/preds/baseline_clinvar.tsv"),
  `Baseline Frequency` = read_tsv("data/dms/preds/baseline_freq.tsv"),
  `UNET (Top)` = read_tsv("data/dms/preds/clinvar_top.tsv"),
  `UNET (Finetune)` = read_tsv("data/dms/preds/clinvar_finetune.tsv"),
  `UNET` = read_tsv("data/dms/preds/unet_freq.tsv"),
  .id = "model"
) %>% pivot_wider(names_from = model, values_from = pred)

# Import EVE data
# Uniprot mapping found using Uniprot search/retrieve IDs service based on DMS Uniprot IDs
uniprot_mapping <- read_tsv("data/dms/uniprot_mapping.tsv") %>%
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
esm_scores <- read_tsv("data/esm/1v_dms_preds.tsv") %>%
  pivot_longer(starts_with("pred_"), names_to = "mut", names_prefix = "pred_", values_to = "s") %>%
  pivot_wider(names_from = "model", values_from = "s") %>%
  mutate(`ESM-1v` = rowMeans(across(c(`1`, `2`, `3`, `4`, `5`)))) %>%
  select(gene=id, position, wt, mut, `ESM-1v`)

# Generate Query for dbNPFS (for applicable human proteins)
# Uses list of cannonical transcripts available from Ensembl
ensembl_cannonical <- read_tsv("data/clinvar/human_grch38_105_ensembl_cannonical.tsv") %>%
  pull(transcript_stable_id) %>%
  str_split("\\.", simplify = TRUE) %>%
  magrittr::extract(1:nrow(.), 1) %>%
  unique()

dbnsfp_query <- read_tsv("data/dms/uniprot_mapping.tsv") %>%
  select(uniprot = Entry, ensembl = `Ensembl transcript`) %>%
  left_join(select(dms, uniprot, position, wt, mut), ., by = "uniprot") %>%
  separate(ensembl, into = str_c("ensembl", 1:13), sep = ";", fill = "right") %>%
  pivot_longer(ensembl1:ensembl13, values_to = "ensembl") %>%
  select(-name) %>%
  filter(!is.na(ensembl), !ensembl == "", !str_count(ensembl, "\\[") > 0, ensembl %in% ensembl_cannonical) %>%
  distinct() %>%
  mutate(query = str_c("Ensembl:", ensembl, ":", wt, position, mut))
# write_lines(dbnsfp_query$query, "data/dms/dbnsfp_ensebml_query.txt")

# Import dbNPFS data
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

dbnsfp <- read_tsv("data/dms/dbnsfp_results.tsv", na = c("NA", "na", ".", "-")) %>%
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
preds <- mutate(dms, gene = str_to_lower(gene)) %>%
  left_join(unet, by = c("gene", "position", "wt", "mut")) %>%
  left_join(esm_scores, by = c("gene", "position", "wt", "mut")) %>%
  left_join(select(eve_scores, uniprot, position, wt, mut, EVE), by = c("uniprot", "position", "wt", "mut")) %>%
  left_join(select(dbnsfp, uniprot, position, wt, mut, PolyPhen2:`GERP++`), by = c("uniprot", "position", "wt", "mut"))

# Calculate correlation
less = c("SIFT4G", "BLOSUM62", "FATHMM", "PROVEAN", "ESM-1v")
correlations <- pivot_longer(preds, SIFT4G:`GERP++`, names_to = "tool", values_to = "pred") %>%
  drop_na() %>%
  mutate(pred = ifelse(tool %in% less, pred, -pred)) %>% # Change preds so lower score means more deleterious like score
  group_by(tool, uniprot) %>%
  summarise(pearson = cor(score, pred, method = "pearson"), spearman = cor(score, pred, method = "spearman"), n = n(),.groups = "drop_last") %>%
  summarise(studies = n(), variants = sum(n), mean_pearson = mean(pearson), sd_pearson = sd(pearson),
            mean_spearman = mean(spearman), sd_spearman = sd(spearman)) %>%
  mutate(stderr_pearson = sd_pearson/sqrt(studies), stderr_spearman = sd_spearman/sqrt(studies))
write_tsv(correlations, "data/dms/correlation.tsv")

correlations_species <- pivot_longer(preds, SIFT4G:`GERP++`, names_to = "tool", values_to = "pred") %>%
  drop_na() %>%
  mutate(pred = ifelse(tool %in% less, pred, -pred)) %>% # Change preds so lower score means more deleterious like score
  group_by(tool, species, uniprot) %>%
  summarise(pearson = cor(score, pred, method = "pearson"), spearman = cor(score, pred, method = "spearman"), n = n(), .groups = "drop")
write_tsv(correlations_species, "data/dms/correlation_species.tsv")

roc <- pivot_longer(preds, SIFT4G:`GERP++`, names_to = "tool", values_to = "pred") %>%
  mutate(class = score < -0.5) %>%
  group_by(tool) %>%
  group_modify(~calc_roc(.x, class, pred, greater = !(.y %in% less), max_steps = 1000)) %>%
  ungroup() %>%
  arrange(desc(auc))
write_tsv(roc, "data/dms/roc.tsv")


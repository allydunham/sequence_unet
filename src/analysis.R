#!/usr/bin/env Rscript
# Shared analysis functions

# Blank ggplot with label
blank_plot <- function(text = ''){
  ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    annotate(geom = 'text', x = 0.5, y = 0.5, label = text) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

# Area under a stepwise curve
step_auc <- function(x, y){
  sum(diff(x) * y[-length(y)])
}

# Calculate ROC curve
calc_roc <- function(tbl, true_col, var_col, greater = TRUE, max_steps = 500, max_matrix_size = 100000){
  true_col <- enquo(true_col)
  var_col <- enquo(var_col)
  
  tbl <- select(tbl, !!true_col, !!var_col) %>%
    drop_na()
  if (nrow(tbl) == 0){
    return(tibble(TP=NA, TN=NA, FP=NA, FN=NA))
  }
  
  true <- pull(tbl, !!true_col)
  var <- pull(tbl, !!var_col)
  
  unique_values <- unique(var)
  if (length(unique_values) > max_steps) {
    steps <- c(-Inf, seq(from = min(unique_values), to = max(unique_values), length.out = max_steps), Inf)
  } else {
    steps <- c(-Inf, sort(unique_values), Inf)
  }
  
  if (length(true) * length(steps) < max_matrix_size) {
    tbl <- calc_roc_matrix(true, var, steps, greater = greater)
  } else {
    tbl <- calc_roc_loop(true, var, steps, greater = greater)
  }
  
  if (greater){
    tbl <- arrange(tbl, desc(thresh))
  } else {
    tbl <- arrange(tbl, thresh)
  }
  
  tbl$auc <- step_auc(tbl$fpr, tbl$tpr)
  return(tbl)
}

calc_roc_matrix <- function(true, var, steps, greater = TRUE) {
  true_mat <- matrix(true, nrow = length(true), ncol = length(steps))
  var_mat <- matrix(var, nrow = length(var), ncol = length(steps))
  thresh_mat <- matrix(steps, nrow = length(var), ncol = length(steps), byrow = TRUE)
  
  if (greater){
    preds <- var_mat >= thresh_mat
  } else {
    preds <- var_mat <= thresh_mat
  }
  
  tp <- colSums(preds & true_mat)
  tn <- colSums(!preds & !true_mat)
  fp <- colSums(preds & !true_mat)
  fn <- colSums(!preds & true_mat)
  tbl <- tibble(thresh = steps, tp = tp, tn = tn, fp = fp, fn = fn,
                tpr = tp / (tp + fn),
                tnr = tn / (tn + fp),
                fpr = fp / (tn + fp),
                precision = tp / (tp + fp))
  return(tbl)
}

calc_roc_loop <- function(true, var, steps, greater = TRUE) {
  if (greater){
    tp <- map_int(steps, ~sum((var >= .) & true))
    tn <- map_int(steps, ~sum(!(var >= .) & !true))
    fp <- map_int(steps, ~sum((var >= .) & !true))
    fn <- map_int(steps, ~sum(!(var >= .) & true))
  } else {
    tp <- map_int(steps, ~sum((var <= .) & true))
    tn <- map_int(steps, ~sum(!(var <= .) & !true))
    fp <- map_int(steps, ~sum((var <= .) & !true))
    fn <- map_int(steps, ~sum(!(var <= .) & true))
  }
  
  tbl <- tibble(thresh = steps, tp = tp, tn = tn, fp = fp, fn = fn,
                tpr = tp / (tp + fn),
                tnr = tn / (tn + fp),
                fpr = fp / (tn + fp),
                precision = tp / (tp + fp))
  return(tbl)
}


auc_labeled_model <- function(model, auc){
  out <- str_c(model, " (AUC = ", signif(auc, 2), ")")
  
  ord <- order(auc, decreasing = TRUE)
  levs <- unique(out[ord])
  
  return(factor(out, levels = levs))
}
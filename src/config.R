#!/usr/bin/env Rscript
# Config for R part of predictor project, containing shared functions, packages and setup

# Import packages
library(tidyverse)
library(broom)
library(ggpubr)
library(ggtext)
library(multipanelfigure)

# Custom packages - available at github.com/allydunham
library(plotlistr) 

### GGPlot theme ###
# clean with centered title by default
theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

### Constants ###
TOOL_COLOURS <- c(`PreGraph UNET`="#cab2d6",
                  UNET="#6a3d9a",
                  `UNET (Finetune)` = "#e31a1c",
                  `PreGraph UNET (Finetune)` = "#fb9a99", 
                  `UNET (Top)` = "#1f78b4",
                  `PreGraph UNET (Top)` = "#a6cee3",
                  `Baseline CNN`="#33a02c",
                  `Baseline Frequency` = "#b2df8a",
                  `Baseline ClinVar` = "#33a02c",
                  BLOSUM62="#a65628",
                  SIFT4G="#e6ab02",
                  SPBuild="#e78ac3",
                  `FoldX`="#a65628",
                  `ESM-1b`="#377eb8",
                  `ESM-1b Logits`="#91ceff",
                  `ESM-1b Top Model`="#377eb8")

KINGDOM_COLOURS <- c(Eukaryote = "#7570b3", Archea = "#d95f02", Bacteria = "#1b9e77")

### Functions ###
pretty_p_values <- function(p, breaks = c(0.001, 0.01, 0.05), markdown_exp=FALSE, prefix_p=FALSE){
  breaks <- sort(breaks, decreasing = TRUE)
  break_str <- as.character(breaks)
  gt <- '>'
  lt <- '<'
  
  if (markdown_exp){
    break_str <- str_replace(break_str, "^1?e(-?[0-9]*\\.?[0-9]*)", "10<sup>\\1</sup>")
    gt <- '&gt;'
    lt <- '&lt;'
  }
  
  p_out <- rep(str_c(gt, ' ', breaks[length(breaks)]), length(p))
  for (i in 1:length(breaks)){
    p_out[p < breaks[i]] <- str_c(lt, ' ', break_str[i])
  }
  break_levels <- c(str_c(gt, ' ', break_str[1]), unlist(map(break_str, ~str_c(lt, ' ', .))))
  
  if (prefix_p){
    p_out <- str_c('p ', p_out)
    break_levels <- str_c('p ', break_levels)
  }
  
  return(factor(p_out, levels = break_levels))
}

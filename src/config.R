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

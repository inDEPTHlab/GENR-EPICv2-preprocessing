#!/usr/bin/Rscript

###################################
## Script: DNAm outlier handling
## 
## Purpose: Trim/winsorize DNAm data to manage outliers
##
## Author Name: Anna Grossbach (based on Isabel's script)
##
## Date: 24-03-2026
##
## Last modified: 13-05-2026
###################################

# ---- Dependencies ----
library(tidyverse)

input = "/PATH/TO/INPUT/FOLDER/"
output = "/PATH/TO/OUTPUT/FOLDER/"

# ---- Functions ----
trim <- function(x, probs = NULL, cutpoints = NULL , replace = c(NA, NA), verbose = TRUE){
  dummy = is.integer(x)
  if (!is.null(probs)){
    stopifnot(is.null(cutpoints))
    stopifnot(length(probs)==2)
    cutpoints <- quantile(x, probs, type = 1, na.rm = TRUE)
  } else if (is.null(cutpoints)){
    l <- quantile(x, c(0.25, 0.50, 0.75), type = 1, na.rm = TRUE)
    cutpoints <- c(l[1]-1.5*(l[3]-l[1]), l[3]+1.5*(l[3]-l[1]))
  } else{
    stopifnot(length(cutpoints)==2)
  }
  if (is.integer(x)) cutpoints <- round(cutpoints)
  bottom <- x < cutpoints[1]
  top <- x > cutpoints[2]
  if (verbose){
    message(paste(sum(bottom, na.rm = TRUE),sum(top, na.rm = TRUE), sep='\t'))
  }
  x[bottom] <- replace[1]
  x[top] <- replace[2]
  if (dummy){
    x <- as.integer(x)
  }
  x
}

winsorize <- function(x, probs = NULL, cutpoints = NULL, replace = c(cutpoints[1], cutpoints[2]), verbose = TRUE){
  dummy = is.integer(x)
  if (!is.null(probs)){
    stopifnot(is.null(cutpoints))
    stopifnot(length(probs) == 2)
    cutpoints <- quantile(x, probs, type = 1, na.rm = TRUE)
  } else if (is.null(cutpoints)){
    cutpoints <- quantile(x, c(0.005, 0.995), type = 1, na.rm = TRUE)
  } else{
    stopifnot(length(cutpoints) == 2)
  }
  if (is.integer(x)) cutpoints <- round(cutpoints)
  bottom <- x < cutpoints[1]
  top <- x > cutpoints[2]
  if (verbose){
    message(paste(sum(bottom, na.rm = TRUE), sum(top, na.rm = TRUE), sep = '\t'))
  }
  x[bottom] <- replace[1]
  x[top] <- replace[2]
  if (dummy){
    x <- as.integer(x)
  }
  x
}


# process each data batch and timepoint separately
for(batch in c("450k_age0", "450k_age5", "450k_age9", "epicv1_age0", "epicv2_age0", "epicv2_age5", "epicv2_age9", "epicv2_age13", "epicv2_age17")){

  Age = str_split(batch, "_", simplify=T)[2]
  Array = str_split(batch, "_", simplify=T)[1]
  
  # ---- Input ----
  print(paste("Loading Data"))
  dnam_functional <- readRDS(paste0(input, Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctionalALL_", Age, "_20260406.rds"))
  dnam_quantile <- readRDS(paste0(input, Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantileALL_", Age, "_20260406.rds"))
  
  # ---- Trim ----
  print(paste("Trimming Functionally Normalised Data"))
  dnam_functional_trim <- t(apply(dnam_functional, 1, trim))
  saveRDS(dnam_functional_trim, 
          file = paste0(output, Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctional3IQRNA_", Age, "_20260512.rds"))
  rm(dnam_functional_trim)
  
  print(paste("Trimming Quantile Normalised Data"))
  dnam_quantile_trim <- t(apply(dnam_quantile, 1, trim))
  saveRDS(dnam_quantile_trim, 
          file = paste0(output, Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantile3IQRNA_", Age, "_20260512.rds"))
  rm(dnam_quantile_trim)
  
  # ---- Winsorize ----
  print(paste("Winsorizing Functionally Normalised Data"))
  dnam_functional_win <- t(apply(dnam_functional, 1, winsorize))
  saveRDS(dnam_functional_win, 
          file = paste0(output, Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctional1pctWins_", Age, "_20260423.rds"))
  rm(dnam_functional_win)
  
  print(paste("Winsorizing Quantile Normalised Data"))
  dnam_quantile_win <- t(apply(dnam_quantile, 1, winsorize))
  saveRDS(dnam_quantile_win, 
          file = paste0(output, Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantile1pctWins_", Age, "_20260406.rds"))
  rm(dnam_quantile_win)
}


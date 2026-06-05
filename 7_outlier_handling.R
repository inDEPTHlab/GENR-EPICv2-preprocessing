#!/usr/bin/Rscript

###################################
## Script: DNAm outlier handling
## 
## Purpose: Trim/winsorize DNAm data to manage outliers
##
## Author Name: Anna Grossbach (based on Isabels script)
##
## Date: 24-03-2026
##
## Last modified: 13-05-2026
###################################

# ---- Dependencies ----
library(tidyverse)


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


for(batch in c("genr_450k_age0", "genr_450k_age5", "genr_450k_age9", "genr_epicv1_age0", "genrnext_epicv1_age0")){

  if(str_detect(batch, "age0")){
    Age = "Birth"
  } else if(str_detect(batch, "age5")){
    Age = "5y"
  } else if(str_detect(batch, "age9")){
    Age = "9y"
  }
  
  if(str_detect(batch, "450k")){
    Array = "450k"
  } else if(str_detect(batch, "genr_epicv1")){
    Array = "EPICv1"
  } else if(str_detect(batch, "genrnext_epicv1")){
    Array = "GENRNXT_EPICv1"
  } 
  
  
  # ---- Input ----
  print(paste("Loading Data"))
  dnam_functional <- readRDS(paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctionalALL_", Age, "_20260406.rds"))
  dnam_quantile <- readRDS(paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantileALL_", Age, "_20260406.rds"))
  
  # ---- Trim ----
  print(paste("Trimming Functionally Normalised Data"))
  dnam_functional_trim <- t(apply(dnam_functional, 1, trim))
  saveRDS(dnam_functional_trim, 
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctional3IQRNA_", Age, "_20260512.rds"))
  rm(dnam_functional_trim)
  
  print(paste("Trimming Quantile Normalised Data"))
  dnam_quantile_trim <- t(apply(dnam_quantile, 1, trim))
  saveRDS(dnam_quantile_trim, 
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantile3IQRNA_", Age, "_20260512.rds"))
  rm(dnam_quantile_trim)
  
  # ---- Winsorize ----
  print(paste("Winsorizing Functionally Normalised Data"))
  dnam_functional_win <- t(apply(dnam_functional, 1, winsorize))
  saveRDS(dnam_functional_win, 
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctional1pctWins_", Age, "_20260423.rds"))
  rm(dnam_functional_win)
  
  print(paste("Winsorizing Quantile Normalised Data"))
  dnam_quantile_win <- t(apply(dnam_quantile, 1, winsorize))
  saveRDS(dnam_quantile_win, 
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantile1pctWins_", Age, "_20260406.rds"))
  rm(dnam_quantile_win)
}


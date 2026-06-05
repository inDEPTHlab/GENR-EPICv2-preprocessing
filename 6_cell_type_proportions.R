#!/usr/bin/Rscript

###################################
## Script: Cell type proportions
## 
## Purpose: Calculate cell type proportions
##          from DNAm data.
##
## Author Name: Anna Grossbach
##
## Date: 18-03-2026
##
## Last modified: 23-04-2026
###################################

# ---- Dependencies ----
library(meffil)
library(EpiDISH)
library(argparse)
library(tidyverse)
options(mc.cores=15)


# ---- Functions ----
# Package each cell type proportion calculation in a function
# an generate a plot that can be used as quality/sanity check later  

## ---- UniLIFE ----
unilife <- function(data, type){

  out <- EpiDISH::epidish(beta.m = data, ref.m = centUniLIFE.m, method = "RPC")$estF
  out <- as.data.frame(out)
  out$SampleID <- rownames(out)
  
  saveRDS(out,
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLUniLIFE_", Age, "_20260423.rds"))
  
  p <- out |> 
    pivot_longer(cols = !SampleID,
                 names_to = "celltypes",
                 values_to = "proportion") |> 
    ggplot(aes(y = proportion,
               x = celltypes)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(p, 
         file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLUniLIFE_", Age, "_20260423.png"))

}

## ---- EpiDISH ----
epidish <- function(data, type){

  out <- EpiDISH::epidish(beta.m = data, ref.m = cent12CT450k.m, method = "RPC")$estF
  out <- as.data.frame(out)
  out$SampleID <- rownames(out)
  
  saveRDS(out,
          file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLEpiDISH_", Age, "_20260423.rds"))
  
  p <- out |> 
    pivot_longer(cols = !SampleID,
                 names_to = "celltypes",
                 values_to = "proportion") |> 
    ggplot(aes(y = proportion,
               x = celltypes)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(p, 
         file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLEpiDISH_", Age, "_20260423.png"))
  
}

## ---- Meffil ----
meffil_cells <- function(idats){
  
  # specify which references to use for batch
  if(str_detect(batch, "age0")){
    
    meffil_cell_type_refs <- c("combined cord blood")
  
  } else {
    
    meffil_cell_type_refs <- c("blood gse167998", "blood gse35069 complete", "blood idoloptimized")
    
  }
  
  # create samplesheet from idats
  samplesheet = meffil.create.samplesheet(idats, recursive=TRUE)
  
  # calculate cell type proportions 
  for(ct in meffil_cell_type_refs){
    if(ct == "combined cord blood"){
      ct_label = "CombinedCordBlood" 
    } else if(ct == "blood gse167998") {
      ct_label = "Salas"
    } else if(ct == "blood gse35069 complete") {
      ct_label = "ReiniusComplete"
    } else if(ct == "blood idoloptimized"){
      ct_label = "BloodIdolOptimized"
    }
    
    print(paste(ct))
    
    qc.objects = meffil.qc(samplesheet, verbose=T, cell.type.reference=ct)
    
    out <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
    out <- data.frame(SampleID=row.names(out), out)
    
    saveRDS(out,
            file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4Idats", ct_label, "_", Age, "_20260406.rds"))
    
    p <- out |> 
      pivot_longer(cols = !SampleID,
                   names_to = "celltypes",
                   values_to = "proportion") |> 
      ggplot(aes(y = proportion,
                 x = celltypes)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90))
    
    ggsave(p, 
           file = paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4Idats", ct_label, "_", Age, "_20260406.png"))
    
    
  }
}


# ---- Cell Type Proportion Calculations ----
for(batch in c("genr_450k_age0", "genr_450k_age5", "genr_450k_age9", "genr_epicv1_age0", "genrnext_epicv1_age0")){
  
  print(paste("Calculating cell type proportions for batch", batch))
  
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
  
  
  dnam_functional <- readRDS(paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctionalALL_", Age, "_20260423.rds"))
  dnam_quantile <- readRDS(paste0("/home/a.grossbach/data/genr_dnam/release4_final/", Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantileALL_", Age, "_20260406.rds"))
  
  if(str_detect(batch, "epic")){
    if(batch == "genr_epicv1_age0"){folder = "genr_epicv1"}
    if(batch == "genrnext_epicv1_age0"){folder = "genrnext_epicv1"}
    idats <- paste0("/home/a.grossbach/data/genr_dnam/", folder, "/idats_selected/")
  } else if(str_detect(batch, "450k")){
    timepoint <- strsplit(batch, "_")[[1]][3]
    idats <- paste0("/home/a.grossbach/data/genr_dnam/genr_450k/idats_selected_", timepoint, "/")
  }

  if(str_detect(batch, "age0")){
    
    print(paste("UniLife"))
    unilife(data = dnam_functional, type = "Functional")
    unilife(data = dnam_quantile, type = "Quantile")
    
    rm(dnam_functional)
    rm(dnam_quantile)
    
    print(paste("Meffil"))
    meffil_cells(idats = idats)
    
  } else {
    
    print(paste("UniLife"))
    unilife(data = dnam_functional, type = "Functional")
    unilife(data = dnam_quantile, type = "Quantile")
    
    print(paste("EpiDISH"))
    epidish(data = dnam_functional, type = "Functional")
    epidish(data = dnam_quantile, type = "Quantile")
    
    rm(dnam_functional)
    rm(dnam_quantile)
    
    print(paste("Meffil"))
    meffil_cells(idats = idats)
    
  } 
}

print("Completed")
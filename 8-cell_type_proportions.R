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


input = "/PATH/TO/INPUT/FOLDER/"
output = "/PATH/TO/OUTPUT/FOLDER/"


# ---- Functions ----
## ---- UniLIFE ----
unilife <- function(data, type){

  out <- EpiDISH::epidish(beta.m = data, ref.m = centUniLIFE.m, method = "RPC")$estF
  out <- as.data.frame(out)
  out$SampleID <- rownames(out)
  
  saveRDS(out,
          file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLUniLIFE_", Age, "_20260423.rds"))
  
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
         file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLUniLIFE_", Age, "_20260423.png"))

}

## ---- EpiDISH ----
epidish <- function(data, type){

  out <- EpiDISH::epidish(beta.m = data, ref.m = cent12CT450k.m, method = "RPC")$estF
  out <- as.data.frame(out)
  out$SampleID <- rownames(out)
  
  saveRDS(out,
          file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLEpiDISH_", Age, "_20260423.rds"))
  
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
         file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4", type, "ALLEpiDISH_", Age, "_20260423.png"))
  
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
            file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4Idats", ct_label, "_", Age, "_20260406.rds"))
    
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
           file = paste0(output, Age , "/", Array, "/CellTypes/Child_", Array, "DNAmRelease4Idats", ct_label, "_", Age, "_20260406.png"))
    
    
  }
}


# ---- Cell Type Proportion Calculations ----
# process each data batch and timepoint separately                    
for(batch in c("450k_age0", "450k_age5", "450k_age9", "epicv1_age0", "epicv2_age0", "epicv2_age5", "epicv2_age9", "epicv2_age13", "epicv2_age17")){
  
  print(paste("Calculating cell type proportions for batch", batch))
  
  Age = str_split(batch, "_", simplify=T)[2]
  Array = str_split(batch, "_", simplify=T)[1]  
  idats <- paste0(input, "idats_", Array, "_", Age, "/")
  
  dnam_functional <- readRDS(paste0(input, Age , "/", Array, "/Functional/Child_", Array, "DNAmRelease4BetasFunctionalALL_", Age, "_20260423.rds"))
  dnam_quantile <- readRDS(paste0(input, Age , "/", Array, "/Quantile/Child_", Array, "DNAmRelease4BetasQuantileALL_", Age, "_20260406.rds"))

  if(Age == "age0"){
    
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

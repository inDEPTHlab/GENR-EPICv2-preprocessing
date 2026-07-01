#!/usr/bin/Rscript

###################################
## Script: DNAm Functional Normalisation
## 
## Purpose: Reduce technical variation in DNAm levels
##
## Author Name: Anna Grossbach
##
## Date: 04-03-2026
##
## Last modified: 07-04-2026
###################################

# ---- Dependencies ----
library(argparse)
library(tidyverse)
library(haven)
library(data.table)
library(meffil)
options(mc.cores=15)


# ---- Input ----
parser <- ArgumentParser()

# DNAm batch, e.g. 450k, EPICv1, EPICv2
parser$add_argument(
  "-b",
  "--batch",
  type = "character"
)

# path to IDAT files of data batch
parser$add_argument(
  "-i",
  "--idats",
  type = "character"
)

# path to output files, i.e. normalized dataset
parser$add_argument(
  "-o",
  "--output",
  type = "character"
)

# number of PCs to use for functional normalization
parser$add_argument(
  "-n_pcs",
  "--n_pcs",
  type = "numeric"
)

args <- parser$parse_args()

batch <- args$batch
idats <- args$idats
output <- args$output 
n_pcs <- args$n_pcs

# generate meffil sample sheet
paste("Processing", batch, "samples in directory:", idats)
samplesheet <- meffil.create.samplesheet(idats)

samplesheet$Slide <- as.factor(samplesheet$Slide)

paste("Number of samples:", dim(samplesheet)[1])
paste("Number of PCs to include in normalisation:", n_pcs)


# ---- Quality filtering ----

# default parameters used are
# detection.threshold=0.01
# bead.threshold=3
# sex.cutoff=-2

qc.parameters <- meffil.qc.parameters(
  detectionp.samples.threshold = 0.1,
  detectionp.cpgs.threshold = 0.1,
  beadnum.samples.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 3
)

qc.objects <- meffil.qc(samplesheet, verbose=TRUE)
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)

# Remove outlier samples
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)


# ---- Normalisation ----
# performs two different sets of normlization,
# one without inclusion of random effects
# and one with a random effect for slide

# Perform quantile normalization
paste("Processing quantile normalisation")

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=n_pcs)
norm.objects.rslide <- meffil.normalize.quantiles(qc.objects, random.effects="Slide", number.pcs=n_pcs)

# Generate normalized probe values
paste("Processing functional normalisation")

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
norm.beta.rslide <- meffil.normalize.samples(norm.objects.rslide, cpglist.remove=qc.summary$bad.cpgs$name)

saveRDS(norm.beta, file=file.path(output, paste0(batch, "_norm.rds")))
saveRDS(norm.beta.rslide, file=file.path(output, paste0(batch, "_norm_rslide.rds")))

# Generate normalization report
pcs <- meffil.methylation.pcs(norm.beta)
pcs.rslide <- meffil.methylation.pcs(norm.beta.rslide)

norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs)
meffil.normalization.report(norm.summary, output.file=paste0(output, "reports/", batch, "_normalisation_report.html"))

norm.summary.rslide <- meffil.normalization.summary(norm.objects.rslide, pcs=pcs.rslide)
meffil.normalization.report(norm.summary.rslide, output.file=paste0(output, "reports/", batch, "_normalisation_rslide_report.html"))

paste("Normalisation completed")


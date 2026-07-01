#!/usr/bin/Rscript

###################################
## Script: DNAm QC and PCA
## 
## Purpose: Prepare data for functional normalization:
##          Identify probes to exclude and number of principal components to normalize on.
##
## Author Name: Anna Grossbach
##
## Date: 04-03-2026
##
## Last modified: 06-03-2026
###################################

# ---- Dependencies ----
library(argparse)
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

# path to output files, i.e. QC report
parser$add_argument(
  "-o",
  "--output",
  type = "character"
)

args <- parser$parse_args()

batch <- args$batch
idats <- args$idats
output <- args$output 

# generate meffil sample sheet
paste("Processing", batch, "samples in directory:", idats)
samplesheet <- meffil.create.samplesheet(idats)

paste("Number of samples:", dim(samplesheet)[1])


# ---- Quality control ----
# Background and dye bias correction, sexprediction, cell counts estimates
paste("Generating QC report")

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
meffil.qc.report(qc.summary, output.file=paste0(output, "reports/", batch, "_qc_report.html"))

paste(length(qc.summary$bad.cpgs), "outlier probes identified")
writeLines(qc.summary$bad.cpgs$name, paste0(output, "reports/", batch, "_outlier_probes.txt"))

paste(length(qc.summary$bad.samples), "outlier samples identified")
writeLines(qc.summary$bad.samples$sample.name, paste0(output, "reports/", batch, "_outlier_samples.txt"))

# Remove outlier samples
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)


# ---- PCA ----
# Plot residuals remaining after fitting control matrix to decide on the number PCs
# to include in the normalization below.
paste("Calculating principal components")

pdf(paste0(output, "reports/", batch, "_pc_plot.pdf"))
print(meffil.plot.pc.fit(qc.objects)$plot)
dev.off()

paste("Check plot to decide how many PCs to include in normalization")

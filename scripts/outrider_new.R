#!/usr/bin/env Rscript

library(OUTRIDER)
library(dplyr)
library(SummarizedExperiment)

# ---------------- Command-line arguments ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript outrider.R <matrix_file> <output_file>")
}

mat_file <- args[1]
output_file <- args[2]

# ---------------- Load matrix ----------------
mat <- read.csv(mat_file, sep = "\t", header = TRUE, row.names = 1)

# ---------------- Build SummarizedExperiment ----------------
se <- SummarizedExperiment(assays = list(counts = as.matrix(mat)))
ods <- OutriderDataSet(se)

# ---------------- Filter and run OUTRIDER ----------------
ods <- filterExpression(ods, minCounts = TRUE, filterGenes = TRUE)
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff = NA, zScoreCutoff = 2, all = TRUE)
res$sampleID <- sub("^X", "", res$sampleID)
# ---------------- Save results ----------------
write.table(res, file = output_file, sep = "\t", row.names = TRUE, col.names = NA)

cat("✅ OUTRIDER analysis completed. Results saved to:", output_file, "\n")


#!/usr/bin/env Rscript
library(OUTRIDER)
library(dplyr)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(yaml)

# ---------------- Command-line arguments ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript outrider_new.R <matrix_file> <output_file> [config.yaml]")
}

mat_file    <- args[1]
output_file <- args[2]
config_file <- if (length(args) >= 3) args[3] else NULL

# ---------------- Load config (YAML or defaults) ----------------
defaults <- list(
    q              = "auto",
    max_iterations = 15L,
    pvalue_cutoff  = 0.05,
    min_l2fc       = 1.0
)

if (!is.null(config_file) && file.exists(config_file)) {
    cfg_raw    <- yaml::read_yaml(config_file)
    cfg        <- modifyList(defaults, cfg_raw[["outrider"]])
    cat("Config loaded from:", config_file, "\n")
} else {
    cfg <- defaults
    cat("No config file provided — using defaults.\n")
}

outrider_q     <- as.character(cfg$q)
max_iterations <- as.integer(cfg$max_iterations)
pval_cutoff    <- as.numeric(cfg$pvalue_cutoff)
min_l2fc       <- as.numeric(cfg$min_l2fc)

cat("=== OUTRIDER standalone run ===\n")
cat("Matrix        :", mat_file,       "\n")
cat("Output        :", output_file,    "\n")
cat("q             :", outrider_q,     "\n")
cat("Max iterations:", max_iterations, "\n")
cat("pval cutoff   :", pval_cutoff,    "\n")
cat("Min |l2fc|    :", min_l2fc,       "\n")

# ---------------- Load matrix ----------------
cat("Loading count matrix...\n")
mat <- read.csv(mat_file, sep = "\t", header = TRUE, row.names = 1,
                check.names = FALSE)   # preserve exact sample IDs

cat("Matrix dimensions:", nrow(mat), "genes x", ncol(mat), "samples\n")

# Log which samples come from the current run vs historical runs
# Current run samples are passed via OUTRIDER_CURRENT_SAMPLES env var (comma-separated)
current_samples_env <- Sys.getenv("OUTRIDER_CURRENT_SAMPLES", unset = "")
if (nchar(current_samples_env) > 0) {
    current_samples  <- strsplit(current_samples_env, ",")[[1]]
    historical_samples <- setdiff(colnames(mat), current_samples)
    cat("Current run samples  :", length(current_samples),  "\n")
    cat("Historical samples   :", length(historical_samples), "\n")
    cat("Total samples in model:", ncol(mat), "\n")
} else {
    cat("Total samples in model:", ncol(mat), "\n")
}

# Defensive: require at least 5 samples for a meaningful model
if (ncol(mat) < 5) warning("Fewer than 5 samples — OUTRIDER results may be unreliable.")

# ---------------- Build OutriderDataSet ----------------
se  <- SummarizedExperiment(assays = list(counts = as.matrix(mat)))
ods <- OutriderDataSet(se)

# ---------------- Filter low-expressed genes ----------------
ods <- filterExpression(ods, minCounts = TRUE, filterGenes = TRUE)
cat("Genes after expression filter:", nrow(ods), "\n")
if (nrow(ods) == 0) stop("No genes remaining after expression filter.")

# ---------------- Parallel backend ----------------
n_threads <- max(1L, as.integer(Sys.getenv("OUTRIDER_THREADS", unset = "1")))

select_bpparam <- function(n) {
    tryCatch({
        bp <- SnowParam(n)
        bpstart(bp)
        cat("Parallel backend: SnowParam with", n, "workers\n")
        return(bp)
    }, error = function(e) {
        warning("SnowParam failed; trying MulticoreParam...")
    })
    if (.Platform$OS.type == "unix" && n > 1) {
        tryCatch({
            bp <- MulticoreParam(n)
            cat("Parallel backend: MulticoreParam with", n, "workers\n")
            return(bp)
        }, error = function(e) {
            warning("MulticoreParam failed; falling back to SerialParam.")
        })
    }
    cat("Parallel backend: SerialParam (single-threaded)\n")
    SerialParam()
}

bp <- select_bpparam(n_threads)
register(bp)

# ---------------- Estimate or set encoding dimension (q) ----------------
# Config comment mentions findEncodingDim() but the correct exported API
# is estimateBestQ() — they are equivalent; findEncodingDim() is internal.
if (outrider_q == "auto") {
    cat("Estimating optimal encoding dimension via estimateBestQ() (OHT)...\n")
    ods   <- estimateBestQ(ods, BPPARAM = bp)
    q_val <- getBestQ(ods)
    cat("Estimated optimal q:", q_val, "\n")
} else {
    q_val <- as.integer(outrider_q)
    cat("Using fixed q:", q_val, "\n")
}

# ---------------- Fit OUTRIDER model ----------------
set.seed(42)
cat("Fitting OUTRIDER (q =", q_val, ", iterations =", max_iterations, ")...\n")
ods <- OUTRIDER(
    ods,
    q          = q_val,
    iterations = max_iterations,
    BPPARAM    = bp
)

tryCatch(bpstop(bp), error = function(e) NULL)
register(SerialParam())

# ---------------- Extract results (unfiltered) ----------------
cat("Extracting results...\n")
res <- as.data.table(
    results(ods,
            padjCutoff = 1,
            l2fcCutoff = 0,
            all        = TRUE)
)

# Strip leading "X" added by R to numeric sample names
res[, sampleID := sub("^X", "", sampleID)]

# Normalise column names across OUTRIDER versions
# older: rawcounts / padjust   newer: rawCounts / padjValue
setnames(res,
         old = c("rawcounts", "padjust"),
         new = c("rawCounts", "padjValue"),
         skip_absent = TRUE)

# Recompute aberrant flag using config thresholds
# (results() with padjCutoff=1 sets aberrant=FALSE everywhere)
res[, aberrant := padjValue < pval_cutoff & abs(l2fc) >= min_l2fc]

# ---------------- Attach gene symbol if available ----------------
if (!"hgnc_symbol" %in% names(res)) {
    rd      <- as.data.table(as.data.frame(rowData(ods)), keep.rownames = "rd_geneID")
    sym_col <- intersect(c("hgnc_symbol", "gene_name", "gene_symbol"), names(rd))
    if (length(sym_col) > 0) {
        rd_sub <- rd[, c("rd_geneID", sym_col[1]), with = FALSE]
        setnames(rd_sub, c("rd_geneID", sym_col[1]), c("geneID", "hgnc_symbol"))
        res <- merge(res, rd_sub, by = "geneID", all.x = TRUE)
    } else {
        res[, hgnc_symbol := NA_character_]
    }
}

# ---------------- Canonical column order ----------------
required_cols    <- c("sampleID", "geneID", "hgnc_symbol",
                      "rawCounts", "normcounts", "l2fc",
                      "zScore", "pValue", "padjValue", "aberrant")
present_required <- intersect(required_cols, names(res))
extra_cols       <- setdiff(names(res), required_cols)
setcolorder(res, c(present_required, extra_cols))

# ---------------- Write unfiltered results (full table) ----------------
# Mirrors {case_id}.outrider.tab in the Snakemake pipeline
output_tab <- sub("(\\.[^.]+)?$", "_all.tsv", output_file, perl = TRUE)
fwrite(res, output_tab, sep = "\t")
cat("Unfiltered results (all genes):", nrow(res), "→", output_tab, "\n")

# ---------------- Write filtered results (aberrant only) ----------------
# Mirrors {case_id}_outrider_results.csv in the Snakemake pipeline
res_filtered <- res[padjValue < pval_cutoff & abs(l2fc) >= min_l2fc]
cat("Aberrant genes (padjValue <", pval_cutoff, "& |l2fc| >=", min_l2fc, "):",
    nrow(res_filtered), "\n")
fwrite(res_filtered, output_file, sep = "\t")

cat("✅ OUTRIDER analysis completed. Results saved to:", output_file, "\n")

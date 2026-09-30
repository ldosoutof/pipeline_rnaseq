#!/usr/bin/env Rscript
# fraser.R — FRASER2 analysis for the normal (HTSeq) pipeline
#
# Usage:
#   Rscript fraser.R <count_dir> <config_file> <output_file> [hyperparams.yaml]
#
# Env vars:
#   FRASER_THREADS  — number of parallel workers (default: 4)
#
# Output:
#   <output_file>   — TSV of all junctions (unfiltered, padjCutoff=NA)
#                     Compatible with annotation_fraser_test4.py (fwrite, no row index)

library(FRASER)
library(data.table)
library(BiocParallel)

# ── CLI arguments ─────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript fraser.R <count_dir> <config_file> <output_file> [hyperparams.yaml]")
}

count_dir   <- args[1]
config_file <- args[2]
output_file <- args[3]
yaml_file   <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL

# ── YAML config (optional, same format as fraser_newVer.R) ───────────────────
defaults <- list(
    strand_specific              = "reverse",
    min_expression_in_one_sample = 10L,
    min_delta_psi                = 0.0,
    do_filter                    = TRUE,
    q                            = 2L,
    implementation               = "PCA",
    grch_version                 = 38L,
    padj_cutoff                  = NA_real_,     # unfiltered output
    delta_psi_cutoff             = NA_real_
)

if (!is.null(yaml_file) && file.exists(yaml_file)) {
    cfg_raw <- yaml::read_yaml(yaml_file)
    cfg     <- modifyList(defaults, cfg_raw[["fraser"]])
    cat("Config loaded from:", yaml_file, "\n")
} else {
    cfg <- defaults
    cat("No YAML config — using defaults.\n")
}

strand_specific              <- as.character(cfg$strand_specific)
min_expression_in_one_sample <- as.integer(cfg$min_expression_in_one_sample)
min_delta_psi                <- as.numeric(cfg$min_delta_psi)
do_filter                    <- as.logical(cfg$do_filter)
fraser_q                     <- as.integer(cfg$q)
implementation               <- as.character(cfg$implementation)
grch_version                 <- as.integer(cfg$grch_version)

# ── Parallel backend ──────────────────────────────────────────────────────────
n_workers <- max(1L, as.integer(Sys.getenv("FRASER_THREADS", unset = "4")))

select_bpparam <- function(n) {
    if (n > 1L) {
        tryCatch({
            bp <- SnowParam(n)
            cat("Parallel backend: SnowParam with", n, "workers\n")
            return(bp)
        }, error = function(e) {
            warning("SnowParam failed; trying MulticoreParam...")
        })
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

bp <- select_bpparam(n_workers)
register(bp)

cat("=== FRASER normal run ===\n")
cat("Count dir   :", count_dir,   "\n")
cat("Config      :", config_file, "\n")
cat("Output      :", output_file, "\n")
cat("Strand      :", strand_specific, "\n")
cat("Min expr    :", min_expression_in_one_sample, "\n")
cat("Min dPSI    :", min_delta_psi, "\n")
cat("q           :", fraser_q, "\n")
cat("Workers     :", n_workers, "\n")

# ── FRASER analysis ───────────────────────────────────────────────────────────
sampleTable <- fread(config_file)
if (nrow(sampleTable) < 5) warning("Fewer than 5 samples — FRASER results may be unreliable.")

settings <- FraserDataSet(colData = sampleTable, workingDir = count_dir)
strandSpecific(settings) <- strand_specific

fds <- countRNAData(settings)
fds <- calculatePSIValues(fds)
fds <- filterExpressionAndVariability(
    fds,
    minExpressionInOneSample = min_expression_in_one_sample,
    minDeltaPsi              = min_delta_psi,
    filter                   = do_filter
)

fds <- FRASER(fds, q = fraser_q, implementation = implementation, BPPARAM = bp)

tryCatch(bpstop(bp), error = function(e) NULL)
register(SerialParam())

fds <- annotateRanges(fds, GRCh = grch_version)

# ── Extract and write results (unfiltered) ───────────────────────────────────
res <- as.data.table(results(fds, padjCutoff = NA, deltaPsiCutoff = NA))
res[, sampleID := sub("^X", "", sampleID)]

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(res, output_file, sep = "\t")
cat("✅ FRASER done —", nrow(res), "junctions →", output_file, "\n")

#!/usr/bin/env Rscript
library(FRASER)
library(data.table)
library(BiocParallel)
library(yaml)

# ─────────────────────────────────────────────
# CLI ARGUMENTS
# ─────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript fraser_newVersion.R <count_dir> <config_file> [out_dir] [hyperparams.yaml]")
}

count_dir   <- args[1]
config_file <- args[2]

# args[3] is out_dir if it looks like a path (contains / or starts with .)
# args[3] is yaml_file if it ends with .yaml/.yml
# Supports both old 3-arg form (count_dir config yaml) and new 4-arg form
# (count_dir config out_dir yaml)
if (length(args) >= 3 && grepl("\\.(yaml|yml)$", args[3], ignore.case = TRUE)) {
    yaml_file <- args[3]
    out_dir   <- count_dir
} else if (length(args) >= 3) {
    out_dir   <- args[3]
    yaml_file <- if (length(args) >= 4) args[4] else NULL
} else {
    out_dir   <- count_dir
    yaml_file <- NULL
}

# ─────────────────────────────────────────────
# LOAD CONFIG (YAML or defaults)
# ─────────────────────────────────────────────

defaults <- list(
    strand_specific              = "reverse",
    min_expression_in_one_sample = 10L,
    min_delta_psi                = 0.0,
    do_filter                    = TRUE,
    q                            = 2L,
    implementation               = "PCA",
    grch_version                 = 38L,
    padj_cutoff                  = 0.05,
    delta_psi_cutoff             = 0.1,
    strip_leading_x              = TRUE
)

if (!is.null(yaml_file) && file.exists(yaml_file)) {
    cfg_raw <- yaml::read_yaml(yaml_file)
    cfg     <- modifyList(defaults, cfg_raw[["fraser"]])
    cat("Config loaded from:", yaml_file, "\n")
} else {
    cfg <- defaults
    cat("No YAML config provided — using defaults.\n")
}

strand_specific              <- as.character(cfg$strand_specific)
min_expression_in_one_sample <- as.integer(cfg$min_expression_in_one_sample)
min_delta_psi                <- as.numeric(cfg$min_delta_psi)
do_filter                    <- as.logical(cfg$do_filter)
fraser_q                     <- as.integer(cfg$q)
implementation               <- as.character(cfg$implementation)
grch_version                 <- as.integer(cfg$grch_version)
padj_cutoff                  <- as.numeric(cfg$padj_cutoff)
delta_psi_cutoff             <- as.numeric(cfg$delta_psi_cutoff)
strip_leading_x              <- as.logical(cfg$strip_leading_x)

# Output paths — written to out_dir (Snakemake-declared paths when called from pipeline,
# or count_dir for backward-compatible standalone use)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_file     <- file.path(out_dir, "fraser_results_aberrant.tsv")
output_file_all <- file.path(out_dir, "fraser_results_all.tsv")

cat("=== FRASER standalone run ===\n")
cat("Count dir              :", count_dir,                      "\n")
cat("Output dir             :", out_dir,                        "\n")
cat("Sample config          :", config_file,                    "\n")
cat("Strand specificity     :", strand_specific,                "\n")
cat("Min expression         :", min_expression_in_one_sample,   "\n")
cat("Min delta PSI          :", min_delta_psi,                  "\n")
cat("Filter                 :", do_filter,                      "\n")
cat("q (latent factors)     :", fraser_q,                       "\n")
cat("Implementation         :", implementation,                  "\n")
cat("GRCh version           :", grch_version,                   "\n")
cat("padjCutoff             :", padj_cutoff,                    "\n")
cat("deltaPsiCutoff         :", delta_psi_cutoff,               "\n")

# ─────────────────────────────────────────────
# PARALLEL BACKEND
# ─────────────────────────────────────────────

n_workers <- max(1L, as.integer(Sys.getenv("FRASER_THREADS", unset = "10")))

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

bp <- select_bpparam(n_workers)
register(bp)

# ─────────────────────────────────────────────
# LOAD DATA & COUNT
# ─────────────────────────────────────────────

cat("Loading sample table...\n")
sampleTable <- fread(config_file)
if (nrow(sampleTable) < 5) warning("Fewer than 5 samples — FRASER results may be unreliable.")

settings <- FraserDataSet(colData = sampleTable, workingDir = count_dir)
strandSpecific(settings) <- strand_specific

cat("Counting RNA data...\n")
fds <- countRNAData(settings)

# ─────────────────────────────────────────────
# PSI VALUES & FILTERING
# ─────────────────────────────────────────────

cat("Calculating PSI values...\n")
fds <- calculatePSIValues(fds)

cat("Filtering by expression and variability...\n")
fds <- filterExpressionAndVariability(
    fds,
    minExpressionInOneSample = min_expression_in_one_sample,
    minDeltaPsi              = min_delta_psi,
    filter                   = do_filter
)
cat("Junctions after filter:", nrow(fds), "\n")
if (nrow(fds) == 0) stop("No junctions remaining after expression filter.")

# ─────────────────────────────────────────────
# FIT FRASER MODEL
# ─────────────────────────────────────────────

# Re-register in case the filter step changed state
bp <- select_bpparam(n_workers)
register(bp)

set.seed(42)
cat("Fitting FRASER model (q =", fraser_q, ", implementation =", implementation, ")...\n")
fds <- FRASER(fds, q = fraser_q, implementation = implementation)

tryCatch(bpstop(bp), error = function(e) NULL)
register(SerialParam())

# ─────────────────────────────────────────────
# ANNOTATE & EXTRACT RESULTS
# ─────────────────────────────────────────────

cat("Annotating ranges (GRCh", grch_version, ")...\n")
fds <- annotateRanges(fds, GRCh = grch_version)

cat("Extracting results (unfiltered)...\n")
res <- as.data.table(results(fds, padjCutoff = NA, deltaPsiCutoff = NA))

if (strip_leading_x) {
    res[, sampleID := sub("^X", "", sampleID)]
}

# Normalise column names across FRASER versions
setnames(res,
         old = c("padjust", "psi"),
         new = c("padjValue", "psiValue"),
         skip_absent = TRUE)

# Recompute aberrant flag using config thresholds
res[, aberrant := padjValue < padj_cutoff & abs(deltaPsi) >= delta_psi_cutoff]

# ─────────────────────────────────────────────
# CANONICAL COLUMN ORDER
# ─────────────────────────────────────────────

required_cols    <- c("sampleID", "seqnames", "start", "end", "strand",
                      "type", "psiValue", "deltaPsi",
                      "pValue", "padjValue", "zScore", "aberrant")
present_required <- intersect(required_cols, names(res))
extra_cols       <- setdiff(names(res), required_cols)
setcolorder(res, c(present_required, extra_cols))

# ─────────────────────────────────────────────
# WRITE OUTPUTS
# ─────────────────────────────────────────────

fwrite(res, output_file_all, sep = "\t")
cat("Unfiltered results (all junctions):", nrow(res), "→", output_file_all, "\n")

res_filtered <- res[padjValue < padj_cutoff & abs(deltaPsi) >= delta_psi_cutoff]
cat("Aberrant junctions (padjValue <", padj_cutoff, "& |deltaPsi| >=", delta_psi_cutoff, "):",
    nrow(res_filtered), "\n")
fwrite(res_filtered, output_file, sep = "\t")

cat("✅ FRASER analysis completed. Results saved to:", output_file, "\n")

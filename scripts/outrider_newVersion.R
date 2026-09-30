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
    # Les colonnes de la matrice ont un format LONG (26D0643-ROU-Jul-...-PUROMOINS)
    # alors que current_samples est au format COURT (26D0643). On normalise les
    # deux au préfixe court (tout avant le 1er tiret/underscore) avant de comparer.
    short_id        <- function(x) sub("[-_].*$", "", x)
    col_short       <- short_id(colnames(mat))
    cur_short       <- short_id(current_samples)
    matched_current <- intersect(cur_short, col_short)
    historical_cols <- colnames(mat)[!(col_short %in% cur_short)]
    cat("Current run samples  :", length(current_samples),  "\n")
    cat("  ...dont présents dans la matrice :", length(matched_current),
        "/", length(current_samples), "\n")
    cat("Historical samples   :", length(historical_cols), "\n")
    cat("Total samples in model:", ncol(mat), "\n")
    # Garde-fou : si AUCUN sample courant ne matche (même après normalisation),
    # l'appariement d'ID est réellement rompu -> on avertit.
    if (length(matched_current) == 0L) {
        warning("AUCUN échantillon du run courant ne correspond aux colonnes de la matrice ",
                "(même après normalisation au préfixe court). ",
                "Exemples colonnes : ", paste(head(colnames(mat), 3), collapse = ", "),
                " | attendus : ", paste(head(current_samples, 3), collapse = ", "))
    }
} else {
    cat("Total samples in model:", ncol(mat), "\n")
}

# Defensive: require at least 5 samples for a meaningful model
if (ncol(mat) < 5) warning("Fewer than 5 samples — OUTRIDER results may be unreliable.")

# ---------------- Build OutriderDataSet ----------------
# CORRECTIF : fournir explicitement colData$sampleID (= colonnes de la matrice)
# pour supprimer le warning "No sampleID was specified" et garantir que results()
# porte les vrais identifiants quelle que soit la version d'OUTRIDER.
sample_ids <- colnames(mat)
se  <- SummarizedExperiment(
    assays  = list(counts = as.matrix(mat)),
    colData = S4Vectors::DataFrame(sampleID  = sample_ids,
                                   row.names = sample_ids)
)
ods <- OutriderDataSet(se)
# Garde-fou : échouer franchement si l'appariement ID est rompu plutôt que de
# produire des sorties faussement nommées (exigence diagnostic).
stopifnot(identical(as.character(colData(ods)$sampleID), sample_ids))

# ---------------- Filter low-expressed genes ----------------
# Remove genes with zeros in more than 30% of samples before OUTRIDER
# This prevents "Every gene contains at least one zero" error with mixed cohorts
count_mat <- assay(ods, "counts")

# ── Liste COMPLÈTE des gènes de l'annotation (AVANT tout filtre) ─────────────
# Sert à réindexer la sortie finale : tous les patients / tous les runs auront
# EXACTEMENT les mêmes lignes de gènes (mêmes gènes, même ordre), pour permettre
# la comparaison ligne à ligne. Les gènes non testés (retirés par les filtres
# ci-dessous, nécessaires au bon fonctionnement du modèle OUTRIDER) seront
# présents dans la sortie avec des valeurs NA (= "non évalué", honnête).
ALL_GENES <- rownames(count_mat)
cat("Gènes totaux dans l'annotation (avant filtre) :", length(ALL_GENES), "\n")

zero_frac  <- rowMeans(count_mat == 0)
ods        <- ods[zero_frac <= 0.3, ]
cat("Genes after zero-fraction filter (<=30% zeros):", nrow(ods), "\n")

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

# ---------------- Size factors (poscounts) ----------------
# Avec une grande cohorte (centaines d'échantillons), CHAQUE gène a au moins un
# zéro. La méthode estimateSizeFactors d'OUTRIDER calcule
# loggeomeans = rowMeans(log(counts)) AVANT de filtrer les zéros : un seul zéro
# dans une ligne -> log(0) = -Inf -> loggeomeans = -Inf pour ce gène. Si tous
# les gènes ont >=1 zéro, all(is.infinite(loggeomeans)) -> stop :
#   "Every gene contains at least one zero, cannot compute log geometric means"
#
# CORRECTIF : calculer les size factors avec la fonction DESeq2
# estimateSizeFactorsForMatrix(type = "poscounts"), qui gère nativement les
# zéros (moyenne géométrique sur les comptes positifs uniquement), puis les
# injecter dans l'ODS. La décomposition controlForConfounders()/fit() qui suit
# empêche OUTRIDER de ré-estimer (et donc de replanter).
# NB : on n'utilise PAS library(DESeq2) (qui masquerait results()/counts() etc.
# d'OUTRIDER) ; on appelle la fonction qualifiée DESeq2::estimateSizeFactorsForMatrix.
cat("Estimating size factors (DESeq2 poscounts, gère les zéros)...\n")
sf <- DESeq2::estimateSizeFactorsForMatrix(counts(ods), type = "poscounts")
# Remplacer d'éventuels NA / 0 par la médiane des facteurs valides
bad <- is.na(sf) | sf <= 0
if (any(bad)) {
    sf[bad] <- median(sf[!bad])
    cat("  ", sum(bad), "size factor(s) NA/0 remplacé(s) par la médiane\n")
}
sizeFactors(ods) <- sf
cat("  size factors (3 premiers) :", paste(round(head(sf, 3), 3), collapse = ", "), "\n")

# ---------------- Estimate or set encoding dimension (q) ----------------
# Config comment mentions findEncodingDim() but the correct exported API
# is estimateBestQ() — they are equivalent; findEncodingDim() is internal.
if (outrider_q == "auto") {
    cat("Estimating optimal encoding dimension via estimateBestQ() (OHT)...\n")
    # Compatibilité multi-versions d'OUTRIDER : les versions récentes de
    # estimateBestQ() n'acceptent plus l'argument BPPARAM (parallélisme géré
    # en interne), les anciennes oui. On teste la signature et on appelle en
    # conséquence, plutôt que d'échouer sur "unused argument (BPPARAM = bp)".
    if ("BPPARAM" %in% names(formals(estimateBestQ))) {
        ods <- estimateBestQ(ods, BPPARAM = bp)
    } else {
        ods <- estimateBestQ(ods)
    }
    q_val <- getBestQ(ods)
    cat("Estimated optimal q:", q_val, "\n")
} else {
    q_val <- as.integer(outrider_q)
    cat("Using fixed q:", q_val, "\n")
}

# ---------------- Fit OUTRIDER model ----------------
# On décompose le wrapper OUTRIDER() en ses étapes explicites pour PRÉSERVER
# les facteurs de taille calculés en mode "poscounts" ci-dessus. Le wrapper
# OUTRIDER() rappelle estimateSizeFactors() avec la méthode "ratio" par défaut,
# ce qui écraserait nos facteurs et reproduirait l'erreur de moyenne géométrique.
set.seed(42)
cat("Controlling for confounders (q =", q_val, ")...\n")
ods <- controlForConfounders(
    ods,
    q          = q_val,
    iterations = max_iterations,
    BPPARAM    = bp
)

cat("Fitting OUTRIDER model and computing p-values...\n")
ods <- fit(ods, BPPARAM = bp)
ods <- computePvalues(ods, BPPARAM = bp)
ods <- computeZscores(ods)

tryCatch(bpstop(bp), error = function(e) NULL)
register(SerialParam())

# ---------------- Extract results (unfiltered) ----------------
# CORRECTIF : results() (OutriderDataSet) est QUALIFIÉ OUTRIDER:: car DESeq2
# (utilisé via DESeq2::estimateSizeFactorsForMatrix) exporte aussi un results().
# all = TRUE => assemble la table COMPLÈTE samples x genes (cf. doc OUTRIDER :
# "If TRUE all results are assembled resulting in a data.table of length
# samples x genes"). INDISPENSABLE pour obtenir TOUS les échantillons dans la
# table _all.tsv (et donc un fichier par échantillon dans filesbysample/).
# Sans all=TRUE, results() ne renvoie que les événements significatifs -> seuls
# les échantillons porteurs d'un hit apparaissent (~16 au lieu de la cohorte).
cat("Extracting results...\n")
res <- as.data.table(
    OUTRIDER::results(ods,
            padjCutoff   = 1,
            zScoreCutoff = 0,
            all          = TRUE)
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

# ── Réindexation sur TOUS les gènes de l'annotation ─────────────────────────
# Garantit un ensemble de gènes IDENTIQUE (même nombre, même ordre) pour tous
# les patients et tous les runs → comparaison ligne à ligne possible. Chaque
# échantillon reçoit la grille complète ALL_GENES ; les gènes non testés par
# OUTRIDER (filtrés pour la validité du modèle) sont présents avec padjValue,
# l2fc, zScore, pValue = NA et aberrant = FALSE (non évalué, pas "non aberrant").
all_samples <- unique(res$sampleID)
grid <- CJ(sampleID = all_samples, geneID = ALL_GENES, sorted = FALSE)
res_full <- merge(grid, res, by = c("sampleID", "geneID"), all.x = TRUE)
# Les gènes non testés : marquer aberrant = FALSE (ils n'ont pas été évalués,
# donc ne sont pas déclarés aberrants), le reste reste NA.
res_full[is.na(aberrant), aberrant := FALSE]
# Ordre des colonnes identique à res
setcolorder(res_full, intersect(names(res), names(res_full)))
# Ordre des lignes : par échantillon puis ordre canonique des gènes
res_full[, geneID := factor(geneID, levels = ALL_GENES)]
setorder(res_full, sampleID, geneID)
res_full[, geneID := as.character(geneID)]

fwrite(res_full, output_tab, sep = "\t")
cat("Résultats réindexés (tous les gènes de l'annotation) :",
    length(ALL_GENES), "gènes ×", length(all_samples), "échantillons =",
    nrow(res_full), "lignes →", output_tab, "\n")

# ---------------- Write filtered results (aberrant only) ----------------
# Mirrors {case_id}_outrider_results.csv in the Snakemake pipeline
res_filtered <- res[padjValue < pval_cutoff & abs(l2fc) >= min_l2fc]
cat("Aberrant genes (padjValue <", pval_cutoff, "& |l2fc| >=", min_l2fc, "):",
    nrow(res_filtered), "\n")
fwrite(res_filtered, output_file, sep = "\t")

cat("✅ OUTRIDER analysis completed. Results saved to:", output_file, "\n")

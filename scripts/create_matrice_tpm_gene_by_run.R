library(dplyr)
library(purrr)
library(rtracklayer)

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript create_matrice_tpm_gene_by_run.R <fastq_parent> <matrix_dir> <gtf> <current_run_tag>")
}

fastq_parent    <- args[1]
matrix_dir      <- args[2]
gtf             <- args[3]
current_run_tag <- args[4]   # ex: 20260610_RUN50_NextSeq_High_16RNASEQ

# ── GTF mapping ───────────────────────────────────────────────────────────────
gtf_data  <- import(gtf)
gene_info <- unique(as.data.frame(gtf_data)[, c("gene_id", "gene_name")])

# ── Collect abundance files (run courant uniquement) ──────────────────────────
root_dir    <- file.path(fastq_parent, "pipeline_v0", "kallisto_bed")
directories <- list.dirs(root_dir, full.names = TRUE)

sampleTable <- list()
for (dir in directories) {
  sample <- basename(dir)
  files  <- list.files(dir, full.names = TRUE,
                       pattern = "abundance_gene_level_counts.tsv")
  for (file in files) {
    sampleTable[[length(sampleTable) + 1]] <- list(sample = sample, path = file)
  }
}

if (length(sampleTable) == 0) {
  stop(paste0("[ERROR] Aucun fichier abundance_gene_level_counts.tsv trouvé sous ", root_dir))
}

sampleTable <- as.data.frame(do.call(rbind, sampleTable))

# ── Validation : tous les échantillons du run courant doivent être présents ───
# Les sous-dossiers kallisto_bed/{sample}/ correspondent aux wildcards {sample}
# du run courant (les FASTQ n'existent que pour ce run).
expected_samples <- basename(directories[directories != root_dir])
expected_samples <- expected_samples[nchar(expected_samples) > 0]
found_samples    <- as.character(sampleTable$sample)
missing_samples  <- setdiff(expected_samples, found_samples)

if (length(missing_samples) > 0) {
  msg <- paste0(
    "[ERROR] ", length(missing_samples), " échantillon(s) du run courant (",
    current_run_tag, ") absents de la matrice TPM Kallisto :\n",
    "  ", paste(missing_samples, collapse = ", "), "\n",
    "  → Vérifiez que kallisto2gene a tourné pour ces échantillons ",
    "(logs sous log/kallisto/<sample>/)"
  )
  message(msg)
  quit(status = 2)
} else {
  message(paste0("[OK] Tous les échantillons du run courant sont présents dans la matrice TPM (",
                 length(found_samples), " échantillons)"))
}

# ── Build TPM matrix ──────────────────────────────────────────────────────────
lst        <- setNames(as.character(sampleTable$path),
                       as.character(sampleTable$sample))

get_sample_id <- function(x) basename(dirname(x))

dfList <- lapply(lst, function(x) {
  message(x)
  read.csv(x, sep = "\t", header = TRUE, row.names = 1,
           col.names = c("genes", get_sample_id(x)))
})

mat <- bind_cols(dfList)

dir.create(matrix_dir, showWarnings = FALSE, recursive = TRUE)
write.table(mat,
            file  = file.path(matrix_dir, "matrice_gene_tpm.tsv"),
            sep   = "\t", quote = FALSE, col.names = NA)

# ── Gene-name mapping ─────────────────────────────────────────────────────────
rownames_mat <- row.names(mat)
mat_bind     <- cbind(mat)
rownames(mat_bind) <- rownames_mat
colnames(mat_bind) <- gsub("^X", "", colnames(mat_bind))

row_names     <- data.frame(Ensembl_ID = rownames(mat_bind))
matched_genes <- merge(gene_info, row_names, by.x = "gene_id", by.y = "Ensembl_ID")
matched_genes$gene_name[is.na(matched_genes$gene_name)] <-
  as.character(matched_genes$gene_id[is.na(matched_genes$gene_name)])

rownames(mat_bind) <- make.unique(matched_genes$gene_name)
colnames(mat_bind) <- gsub("^X", "", colnames(mat_bind))

write.table(mat_bind,
            file  = file.path(matrix_dir, "matrice_gene_tpm_gene.tsv"),
            sep   = "\t", quote = FALSE, col.names = NA)

message(paste0("✅ Matrices TPM sauvegardées dans ", matrix_dir))

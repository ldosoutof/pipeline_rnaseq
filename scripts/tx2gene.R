# tx2gene.R
# Usage: Rscript tx2gene.R <abundance.h5> <gtf_file> <output_counts.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript tx2gene.R <abundance_h5> <gtf_file> <output_counts.tsv>")
}

abundance_file <- args[1]
gtf_file       <- args[2]
output_counts  <- args[3]

if (!file.exists(abundance_file)) {
  stop(paste("Abundance file not found:", abundance_file))
}
if (!file.exists(gtf_file)) {
  stop(paste("GTF file not found:", gtf_file))
}

library(tximport)
library(GenomicFeatures)

sample_path <- dirname(abundance_file)

# Build tx2gene mapping from GTF
txdb    <- makeTxDbFromGFF(file = gtf_file)
k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Import kallisto data
files        <- file.path(sample_path, "abundance.h5")
names(files) <- basename(sample_path)

txi_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                     ignoreTxVersion = TRUE)

# Save gene-level TPM (path owned by Snakemake output:)
# CORRECTIF : on écrit txi_gene$abundance (= TPM), PAS txi_gene$counts.
# L'ancien pipeline écrivait le TPM ici ; le calcul des metrics attend du TPM
# (seuil "TPM > 10", colonne HBA_total en TPM). Écrire $counts gonflait les
# valeurs (~7x, variable selon la profondeur) et faussait HBA_total et %DI_green.
dir.create(dirname(output_counts), recursive = TRUE, showWarnings = FALSE)
write.table(txi_gene$abundance,
            file  = output_counts,
            sep   = "\t", quote = FALSE)

cat("✅ tx2gene done:", output_counts, "\n")

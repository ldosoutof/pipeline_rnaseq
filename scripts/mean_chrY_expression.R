# mean_chrY_expression.R
# Usage: Rscript mean_chrY_expression.R <matrice_tpm_gene.tsv> <chrY_genes.txt> <output.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript mean_chrY_expression.R <matrice_tpm_gene.tsv> <chrY_genes.txt> <output.tsv>")
}

matrix_file <- args[1]
genes_file  <- args[2]
out_file    <- args[3]

if (!file.exists(matrix_file)) stop(paste("Matrix file not found:", matrix_file))
if (!file.exists(genes_file))  stop(paste("Genes file not found:",  genes_file))

mat <- read.table(matrix_file, header = TRUE, row.names = 1,
                  sep = "\t", check.names = FALSE)

genes     <- readLines(genes_file)
genes_ok  <- intersect(genes, rownames(mat))
genes_miss <- setdiff(genes, rownames(mat))

if (length(genes_miss) > 0) {
  message("[WARN] ", length(genes_miss), " chrY gene(s) absent from matrix: ",
          paste(genes_miss, collapse = ", "))
}
if (length(genes_ok) == 0) {
  stop("No chrY genes found in matrix — cannot infer sex. ",
       "Check that matrice_gene_tpm_gene.tsv uses the same gene IDs as ",
       genes_file)
}

chrY_mat <- mat[genes_ok, , drop = FALSE]
mean_tpm <- colMeans(chrY_mat)

out <- data.frame(
  sample        = names(mean_tpm),
  mean_chrY_TPM = mean_tpm,
  row.names     = NULL
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ mean_chrY_expression done:", nrow(out), "samples,",
    length(genes_ok), "/", length(genes), "chrY genes used\n")

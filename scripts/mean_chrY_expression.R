args <- commandArgs(trailingOnly=TRUE)

matrix_file <- args[1]
genes_file  <- args[2]
out_file    <- args[3]

mat <- read.table(matrix_file, header=TRUE, row.names=1,
                  sep="\t", check.names=FALSE)

genes <- readLines(genes_file)
genes <- intersect(genes, rownames(mat))

chrY_mat <- mat[genes, , drop=FALSE]

mean_tpm <- colMeans(chrY_mat)

out <- data.frame(
  sample = names(mean_tpm),
  mean_chrY_TPM = mean_tpm
)

write.table(out, out_file, sep="\t",
            quote=FALSE, row.names=FALSE)


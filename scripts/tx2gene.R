### example : for i in /data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/*RNASEQ/pipeline_v0/kallisto_strand;do echo $i; for j in $i/*; do echo $j;sudo /home/ldosouto/miniconda3/envs/fraser/bin/Rscript tx2gene.R $j/abundance.h5; done; done  ###

# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 1) {
	  stop("Please provide the abundance file path as a command line argument.")
}

# Retrieve the abundance file path from command line argument
abundance_file <- commandArgs(trailingOnly = TRUE)[1]

# Check if the file exists
if (!file.exists(abundance_file)) {
	  stop("The specified abundance file does not exist.")
}

# Extracting directory path
sample_path <- dirname(abundance_file)
print(sample_path)  # This will print the directory path

# Libraries loading
library(GenomicFeatures)
library(tximport)
library(rhdf5)

# Creating a transcript database
txdb <- makeTxDbFromGFF(file = "/dataref/bank/human/annotation/GRCh38/ensembl/current/Homo_sapiens.GRCh38.106.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Run tximport for gene-level counts
txi_gene <- tximport(abundance_file, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Return the gene-level counts or save them to a file
#abundance = tpm
write.table(txi_gene$abundance, file = paste0(sample_path, "/abundance_gene_level_counts.tsv"), sep = "\t", quote = FALSE)

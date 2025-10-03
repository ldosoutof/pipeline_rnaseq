### example : for i in /data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/*RNASEQ/pipeline_v0/kallisto_strand;do echo $i; for j in $i/*; do echo $j;sudo /home/ldosouto/miniconda3/envs/fraser/bin/Rscript tx2gene.R $j/abundance.h5; done; done  ###

library(dplyr)
library(purrr)
library(rtracklayer)


# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 3) {
	  stop("Please provide the output directory of run analysis.")
}

# Retrieve the abundance file path from command line argument
dir <- commandArgs(trailingOnly = TRUE)[1]
folder <- commandArgs(trailingOnly = TRUE)[2]
gtf <- commandArgs(trailingOnly = TRUE)[3]


# Read the GTF file
gtf_data <- import(gtf)
# Extract gene name from the GTF data
gene_info <- as.data.frame(gtf_data)[, c("gene_id", "gene_name")]
# Filter unique gene entries
gene_info <- unique(gene_info)


# Check if the file exists
#if (!file.exists(abundance_file)) {
#	  stop("The specified abundance file does not exist.")
#}

# Initialize an empty list to store sample names and file paths
sampleTable <- list()

# Specify the root directory where gene count files are located
root_dir <-  paste0(dir, "/pipeline_v0/kallisto_bed/")

# List all directories in the root directory
directories <- list.dirs(root_dir, full.names = TRUE)

# Loop through each directory
for (dir in directories) {
  # Get the sample name from the directory name
  sample <- basename(dir)
  # List all files in the directory
  files <- list.files(dir, full.names = TRUE, pattern = "abundance_gene_level_counts.tsv")
  #print(files)
  # Loop through each file
  for (file in files) {
     # Add the sample name and file path to the sampleTable list
     sampleTable[[length(sampleTable) + 1 ]] <- list(sample = sample, path = file)
  }
}
#print(sampleTable)
# Convert the list to a data frame
matrice <- do.call(rbind, sampleTable)
#print(matrice)
sampleTable=as.data.frame(matrice)
#print(sampleTable$sample)

lst <- sampleTable$path
names(lst) <- sampleTable$sample

# Function to extract sample ID from the file path
get_sample_id <- function(x) {
	  dirname <- dirname(x)
  basename(dirname)
}


## making a list of data frames with row and column names while removing last five lines from HTSeq that contains the features and ambiguous read infromation
dfList <- lapply(lst, function(x) {
	   print(x)
           read.csv(x, sep = "\t", header = TRUE, row.names = 1, col.names = c("genes", get_sample_id(x)))})


## combining all the dataframes into single dataframe
mat <- bind_cols(dfList)
mat<-mat[,-1]
print(head(mat))
# Return the gene-level counts or save them to a file
write.table(mat, file = paste0(root_dir, "/marice_gene_tpm.tsv"), sep = "\t", quote = FALSE,col.names=NA)



rownames_mat <-row.names(mat)

   mat_bind=cbind(mat)
   rownames(mat_bind)<-rownames_mat
   colnames(mat_bind) <- gsub("^X", "", colnames(mat_bind))
   row_names <- as.data.frame(rownames(mat_bind))
   names(row_names) <- "Ensembl_ID"

   matched_genes <- merge(gene_info, row_names, by.x = "gene_id", by.y = "Ensembl_ID")
   matched_genes$gene_name[is.na(matched_genes$gene_name)] <- as.character(matched_genes$gene_id[is.na(matched_genes$gene_name)])
   rownames(mat_bind) <- make.unique(matched_genes$gene_name)
   matriceTPM2 <- mat_bind[,-1]
   colnames(matriceTPM2) <- gsub("^X", "", colnames(matriceTPM2))
   #print(head(mat_bind))
   write.table(matriceTPM2, file = paste0(root_dir, "/marice_gene_tpm_gene.tsv"), sep = "\t", quote = FALSE, col.names=NA)


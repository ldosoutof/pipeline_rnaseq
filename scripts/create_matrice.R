### example : for i in /data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/*RNASEQ/pipeline_v0/kallisto_strand;do echo $i; for j in $i/*; do echo $j;sudo /home/ldosouto/miniconda3/envs/fraser/bin/Rscript tx2gene.R $j/abundance.h5; done; done  ###

library(dplyr)


# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
	  stop("Please provide the output directory of run analysis.")
}

# Retrieve the abundance file path from command line argument
dir <- commandArgs(trailingOnly = TRUE)[1]
folder <- commandArgs(trailingOnly = TRUE)[2]
# Check if the file exists
#if (!file.exists(abundance_file)) {
#	  stop("The specified abundance file does not exist.")
#}

# Initialize an empty list to store sample names and file paths
sampleTable <- list()

# Specify the root directory where gene count files are located
root_dir <-  paste0(dir, "/pipeline_v0/htseq/")

# List all directories in the root directory
directories <- list.dirs(root_dir, full.names = TRUE)

# Loop through each directory
for (dir in directories) {
  # Get the sample name from the directory name
  sample <- basename(dir)
  # List all files in the directory
  files <- list.files(dir, full.names = TRUE, pattern = "*gene_counts.txt")
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
print(matrice)
sampleTable=as.data.frame(matrice)
#print(sampleTable$sample)

lst <- sampleTable$path
names(lst) <- sampleTable$sample

## making a list of data frames with row and column names while removing last five lines from HTSeq that contains the features and ambiguous read infromation
dfList <- lapply(lst, function(x) {
	   print(x)
           read.csv(x, nrows = length(count.fields(x)) - 5, sep = "\t", header = FALSE, row.names = 1, col.names = c("genes", tools::file_path_sans_ext(basename(x))))})

## combining all the dataframes into single dataframe
mat <- bind_cols(dfList)


# Return the gene-level counts or save them to a file
write.table(mat, file = paste0(root_dir, "/matrice_gene_counts.tsv"), sep = "\t", quote = FALSE)


sampleTable <- list()

# Loop through each directory
for (dir in directories) {
  # Get the sample name from the directory name
  sample <- basename(dir)
  print(sample)
  # List all files in the directory
  files <- list.files(dir, full.names = TRUE, pattern = "*trans_counts.txt")
  print(files)
      # Loop through each file
      for (file in files) {
	  # Add the sample name and file path to the sampleTable list
	  sampleTable[[length(sampleTable) + 1 ]] <- list(sample = sample, path = file)
      }
}
#head(sampleTable)
# Convert the list to a data frame
matrice <- do.call(rbind, sampleTable)

sampleTable=as.data.frame(matrice)

lst <- sampleTable$path
names(lst) <- sampleTable$sample

## making a list of data frames with row and column names while removing last five lines from HTSeq that contains the features and ambiguous read infromation
dfList <- lapply(lst, function(x) {
                    print(x)
                    read.csv(x, nrows = length(count.fields(x)) - 5, sep = "\t", header = FALSE, row.names = 1, col.names = c("transcripts", tools::file_path_sans_ext(basename(x))))})

## combining all the dataframes into single dataframe
mat_trans <- bind_cols(dfList)


# Return the gene-level counts or save them to a file
write.table(mat_trans, file = paste0(root_dir, "/matrice_trans_counts.tsv"), sep = "\t", quote = FALSE)



# List all files in the directory
files <- list.files(folder, pattern = "_matrice.txt", full.names = TRUE)


# Get the current date
current_date <- Sys.Date()

# Format the date as "YYYYMMDD"
formatted_date <- format(current_date, "%Y%m%d")

# Concatenate the formatted date with the rest of the file name
file_name <- paste0(formatted_date, "_matrice.txt")

print(files)
print(file_name)


# Check if there are any files in the folder
if (length(files) == 0) {
   write.table(mat, paste0(folder, file_name), sep = "\t", quote = FALSE)
} else {
    # Extract dates from filenames
    dates <- as.Date(gsub(".*?([0-9]{8}).*", "\\1", files), format = "%Y%m%d")

    # Find the index of the file with the most recent date
    latest_index <- which.max(dates)

   # Get the filename of the file with the most recent date
   latest_file <- files[latest_index]

   # Print the filename of the latest file
   print(latest_file)
   mat_file=read.csv(latest_file, sep = "\t", header = TRUE, row.names = 1)
   print(head(mat_file))
   mat_bind=cbind(mat_file,mat)
   colnames(mat_bind) <- gsub("^X", "", colnames(mat_bind))
   write.table(mat_bind, paste0(folder, file_name), sep = "\t", quote = FALSE)

}

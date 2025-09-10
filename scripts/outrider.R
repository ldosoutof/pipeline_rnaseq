library(OUTRIDER)
library(dplyr)


# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 3) {
	          stop("Please provide the output directory of run analysis.")
}


mat_file <- commandArgs(trailingOnly = TRUE)[1]
mat=read.csv(mat_file, sep = "\t", header = TRUE,  row.names = 1)


remove_samples_file <- commandArgs(trailingOnly = TRUE)[2]

puromoins_cols <- grep("MOINS", colnames(mat))

mat_puromoins <- mat[, puromoins_cols]

mat_filtered <- mat_puromoins[, !grepl("POLYA", colnames(test))]

# Identify column names starting with the same 6 characters
start_with_same_7 <- substr(colnames(mat_filtered), 1, 7)
duplicate_start <- start_with_same_7[duplicated(start_with_same_7)]

# Identify columns with "CMP" in the header
has_LM <- grepl("LM", colnames(mat_filtered))

# Identify columns to be removed
columns_to_remove <- which(duplicate_start %in% substr(colnames(mat_filtered)[has_LM], 1, 7))

# Remove columns
mat_final <- mat_filtered[, -columns_to_remove]

# Extract the first 7 characters of the column names
new_colnames <- substr(colnames(mat_final), 1, 7)

# Assign the new column names
colnames(mat_final) <- new_colnames

######################
#      OUTRIDER      # 
######################

se <- SummarizedExperiment(assays = list(counts = as.matrix(mat_final)))

ods <- OutriderDataSet(se)

ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

# run full OUTRIDER pipeline (control, fit model, calculate P-values)
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff=NA, zScoreCutoff=2, all=TRUE)

output_file <- commandArgs(trailingOnly = TRUE)[3]

write.table(res, file=output_file, sep="\t")

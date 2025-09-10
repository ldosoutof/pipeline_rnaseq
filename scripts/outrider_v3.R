library(OUTRIDER)
library(dplyr)


# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Please provide the output directory of run analysis.")
}


mat_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

mat=read.csv(mat_file, sep = "\t", header = TRUE, row.names=1 )

head(mat)
has_puromoins <- grepl("PUROMOINS", colnames(mat))
# Keep only columns containing puromoins
mat_bind <- mat[, has_puromoins]


to_del_run17<-c('BAB.SAN.PUROMOINS_gene_counts','X22D2375.1.BED.LIL.LYMPHO.PUROMOINS_gene_counts','X23D0686.4.BAR.EMI.LYMPHO.PUROMOINS_gene_counts','X23D0946.1.LAI.Rom.WDFY3.LYMPHO.PUROMOINS_gene_counts','X23D0981.01.LOR.ZEL.ND.LM.PUROMOINS_gene_counts', 'X23D0984.01.ELI.Fir.ND.LM.PUROMOINS_gene_counts','X23D1172.01.ZAL.Ily.MBD5.LM.PUROMOINS_gene_counts', 'X23D1223.05.LEC.Jul.INTS6.CMP.PUROMOINS_gene_counts', 'X23D1432.01.DRO.Tip.TRIO.CMP.PUROMOINS_gene_counts')
mat_bind<- mat[, !(colnames(mat) %in% to_del_run17)]

to_del_err<-c('X23D1192.01.LUX.Lou.GLUT1.LM.PUROMOINS_gene_counts', 'X23D1375.01.CLE.Luc.ND.LM.PUROMOINS_gene_counts', 'X23D2564.CHE.Lin.PUROMOINS_gene_counts', 'X24D1179.LOU.Lou.PUROMOINS_gene_counts')
mat_bind<- mat_bind[, !(colnames(mat_bind) %in% to_del_err)]

to_del_gtf2i<-c('X23D2544.MOH.JAY.PUROMOINS_gene_counts',
'X23D2564.CHE.Lin.PUROMOINS_gene_counts',
'X23D2578.REY.Jas.PUROMOINS_gene_counts',
'X23D2617.COR.Emm.PUROMOINS_gene_counts',
'X23D2636.SCH.May.PUROMOINS_gene_counts',
'X23D2668.MER.Elo.PUROMOINS_gene_counts',
'X24D0044.CAL.Flo.PUROMOINS_gene_counts',
'X24D0086.DU.Del.PUROMOINS_gene_counts',
'X24D0414.BOU.Abd.PUROMOINS_gene_counts',
'X24D0474.POL.Lyl.PUROMOINS_gene_counts',
'X24D0493.CHA.Lol.PUROMOINS_gene_counts',
'X24D0513.DES.Zel.PUROMOINS_gene_counts',
'X24D0529.GAR.Ber.PUROMOINS_gene_counts',
'X24D0558.RAM.Emm.PUROMOINS_gene_counts',
'X24D0589.BER.Aro.PUROMOINS_gene_counts',
'X24D0660.FOU.Orp.PUROMOINS_gene_counts',
'X24D0698.GAB.Eph.PUROMOINS_gene_counts')

mat_bind<- mat_bind[, !(colnames(mat_bind) %in% to_del_gtf2i)]


colnames(mat_bind) <- gsub("^X", "", colnames(mat_bind))

colnames(mat_bind) <- gsub('_gene_counts', '', colnames(mat_bind))
sort(colnames(mat_bind))
# Identify column names starting with the same 6 characters
start_with_same_6 <- substr(colnames(mat_bind), 1, 7)
duplicate_start <- start_with_same_6[duplicated(start_with_same_6)]

## Identify columns with "CMP" in the header
has_CMP <- grepl("LM|POLYA|BAB|EDTA", colnames(mat_bind))

common_cols <- intersect(duplicate_start, substr(colnames(mat_bind), 1, 7))

columns_to_remove <- which(substr(colnames(mat_bind), 1, 7) %in% common_cols)

columns_to_remove <- columns_to_remove[has_CMP[columns_to_remove]]

# Extract the subset of columns that contain "CMP" and are marked for removal
cmp_cols_to_remove <- colnames(mat_bind)[columns_to_remove]

mat_filtered <- mat_bind[, !colnames(mat_bind) %in% cmp_cols_to_remove]

colnames(mat_filtered) <- substr(colnames(mat_filtered), 1, 7)

#cols_to_remove <- which(colnames(mat_filtered) == "23D2700")

cols_to_remove <- which(colnames(mat_filtered) == "23D2700")
if(length(cols_to_remove) > 0) {
  mat_filteredata2 <- mat_filtered[, -cols_to_remove]
} else {
  mat_filteredata2 <- mat_filtered
}

# Remove columns
#mat_filteredata2 <- mat_filtered[, -cols_to_remove]

summary(rowSums(mat_filteredata2))
summary(colSums(mat_filteredata2))

mat_numeric <- matrix(as.numeric(unlist(mat_filteredata2)),
                      nrow = nrow(mat_filteredata2),
                      dimnames = dimnames(mat_filteredata2))

se <- SummarizedExperiment(assays = list(counts = as.matrix(mat_numeric)))
ods <- OutriderDataSet(se)

ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

# run full OUTRIDER pipeline (control, fit model, calculate P-values)
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff=NA, zScoreCutoff=2, all=TRUE)

# Get the current date
current_date <- Sys.Date()
formatted_date <- format(current_date, "%y%m%d")

# Construct the output file path with the date
#output_file <- paste0("/data2/Exome_analysis_BC/Server_S_exome/RNA-Seq_bioinformatique/5-outrider/outrider_", formatted_date, "_htseq.tab")

# Write the results to a file
write.table(res, file=output_file, sep="\t", row.names=TRUE, col.names=NA)


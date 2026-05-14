library(FRASER)
library(data.table)

if (length(commandArgs(trailingOnly = TRUE)) != 3) {
    stop("Please provide the output directory of run analysis.")
}

count_dir <- commandArgs(trailingOnly = TRUE)[1]
config_file <- commandArgs(trailingOnly = TRUE)[2]
output_file <- commandArgs(trailingOnly = TRUE)[3]


#register(MulticoreParam(workers = 12))
#print(config_file)
BiocParallel::register(SnowParam(workers = 10))

sampleTable <- fread(config_file)
settings <- FraserDataSet(colData=sampleTable, workingDir=count_dir)
strandSpecific(settings) <- "reverse"
fds <- countRNAData(settings)
fds <- calculatePSIValues(fds)
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=10,minDeltaPsi=0.0, filter=TRUE)
BiocParallel::register(SnowParam(workers = 10))
#strandSpecific(fds) <- "reverse"
fds <- FRASER(fds, q=2, implementation="PCA")
fds <- annotateRanges(fds,GRCh=38)
res <- results(fds, padjCutoff=NA, deltaPsiCutoff=NA)
res$sampleID <- sub("^X", "", res$sampleID)
write.table(res, file=output_file, sep="\t")


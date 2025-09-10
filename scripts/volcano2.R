library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(stringr)

# Check if the correct number of command line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Please provide the output directory of run analysis.")
}

outrider_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

# Read the input file
outrider <- read.csv(outrider_file, sep = "\t", header = TRUE)

# Data cleaning and transformation
outrider$pValue <- as.numeric(gsub(",", ".", outrider$pValue))
outrider$zScore <- as.numeric(gsub(",", ".", outrider$zScore))
outrider$l2fc <- as.numeric(gsub(",", ".", outrider$l2fc))
outrider$fc <- 2^outrider$l2fc  # Calculate fold change

genesSpe <- outrider$gene_name
top <- outrider[outrider$gene_name != "", ]

# Plot 1: zScore Volcano Plot
sub_top <- ifelse(
    top$pValue < 0.01 & top$zScore > 3, "up",
    ifelse(top$pValue < 0.01 & top$zScore < -3, "not sig", "down")
)

data <- top %>%
    mutate(
        Expression = case_when(
            zScore > 3 & pValue <= 0.01 ~ "Up-regulated",
            zScore < -3 & pValue <= 0.01 ~ "Down-regulated",
            TRUE ~ "Unchanged"
        )
    )

label_genes <- data %>%
    filter(Expression %in% c("Up-regulated", "Down-regulated")) %>%
    filter(!is.na(Model_Of_Inheritance) & Model_Of_Inheritance != "")

g <- ggplot(data, aes(x = zScore, y = -log10(pValue))) +
    geom_point(aes(color = sub_top)) +
    scale_color_manual("Legend", values = c("grey", "blue", "red"), labels = c("Not sig", "Down", "Up")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top-right", legend.text = element_text(size = 20)) +
    xlim(-10, 10) +
    geom_hline(yintercept = 2, linetype = "dotted") +
    geom_vline(xintercept = c(-3, 3), linetype = "dashed") +
    geom_label_repel(
        data = label_genes,
        aes(label = gene_name),
        size = 5,
        max.overlaps = 12,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.5, "lines")
    )

# Save first plot
png(filename = output_file, width = 2800, height = 1800, res = 300)
print(g)
dev.off()

# Plot 2: l2fc Volcano Plot
data_fc <- top %>%
    mutate(
        Expression = case_when(
            l2fc > 1 & pValue <= 0.01 ~ "Up-regulated",
            l2fc < -1 & pValue <= 0.01 ~ "Down-regulated",
            TRUE ~ "Unchanged"
        )
    )
#head(data_fc)
label_genes_fc <- data_fc %>%
    filter(Expression %in% c("Up-regulated", "Down-regulated")) %>%
    filter(!is.na(Model_Of_Inheritance) & Model_Of_Inheritance != "")
#print(label_genes_fc)

sub_top <- ifelse(
    top$pValue < 0.01 & top$l2fc > 1, "up",
    ifelse(top$pValue < 0.01 & top$l2fc < -1, "not sig", "down")
)

g_fc <- ggplot(data_fc, aes(x = l2fc, y = -log10(pValue))) +
    geom_point(aes(color = sub_top)) +
    scale_color_manual("Legend", values = c("grey", "blue", "red"), labels = c("Not sig", "Down", "Up")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top-right", legend.text = element_text(size = 20)) +
    xlim(-5, 5) +
    geom_hline(yintercept = 2, linetype = "dotted") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_label_repel(
        data = label_genes_fc,
        aes(label = gene_name),
        size = 5,
        max.overlaps = 12,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.5, "lines")
    )

# Save second plot
output_file_fc <- gsub("\\.png$", "_fc.png", output_file)
png(filename = output_file_fc, width = 2800, height = 1800, res = 300)
print(g_fc)
dev.off()


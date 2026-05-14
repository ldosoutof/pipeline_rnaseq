# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]  # Directory containing all files
output_file <- args[2]  # Output file for the first boxplot
output_file_filter <- args[3]  # Output file for the filtered data boxplot
actual_samples <- unlist(strsplit(args[4], ","))  # Comma-separated actual sample IDs

# List all files in the directory (assuming .txt format)
file_list <- list.files(path = directory, pattern = "\\.outrider.tab$", full.names = TRUE)

# Debugging: Print detected files
print("Detected files in directory:")
print(file_list)

# Read files and process data
outrider_list <- lapply(file_list, function(file) {
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }
  print(paste("Reading file:", file))
  tryCatch({
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    colnames(data) <- make.names(colnames(data), unique = TRUE)
    # Convert padjust to numeric, replacing non-numeric values with NA
    if ("padjust" %in% colnames(data)) {
      data$padjust <- as.numeric(gsub(",", ".", data$padjust))
    }
    # Force chromosome to be character
    if ("chromosome" %in% colnames(data)) {
      data$chromosome <- as.character(data$chromosome)
    }

    if ("l2fc" %in% colnames(data)) {
      data$l2fc <- as.numeric(gsub(",", ".", data$l2fc))
    }
    if ("ID_SPICE" %in% colnames(data)) {
      data$ID_SPICE <- as.character(data$ID_SPICE)
    }

    return(data)
  }, error = function(e) {
    message("Error reading file ", file, ": ", e$message)
    return(NULL)
  })
})

# Remove NULL entries from failed file reads
outrider_list <- Filter(Negate(is.null), outrider_list)

# Combine all data into one dataframe
if (length(outrider_list) > 0) {
  outrider <- bind_rows(outrider_list)
} else {
  stop("No valid data files found, exiting.")
}

# Ensure 'sampleID' column exists
if (!"sampleID" %in% colnames(outrider)) {
  stop("Error: 'sampleID' column missing in input data.")
}

# Categorize samples as "Actual Run" or "Older Runs"
outrider <- outrider %>%
  mutate(RunCategory = ifelse(sampleID %in% actual_samples, "Actual Run", "Older Runs"))

# Summarize counts for each sample
count_data <- outrider %>%
  group_by(sampleID, RunCategory) %>%
  summarise(count = n(), .groups = 'drop')

# Identify outliers
outliers <- count_data %>%
  group_by(RunCategory) %>%
  filter(count %in% boxplot.stats(count)$out) %>%
  ungroup()

# Debugging: Print table of sample categories
print("Sample count by category:")
print(table(count_data$RunCategory))

# Create boxplot
boxplot <- ggplot(count_data, aes(x = RunCategory, y = count, fill = RunCategory)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16) +
  geom_text_repel(data = outliers, aes(x = RunCategory, y = count, label = sampleID),
                  size = 3, color = "red", max.overlaps = 10) +
  theme_minimal() +
  labs(title = "Boxplots of Events Counts in Outrider",
       x = "Sample Group",
       y = "Count") +
  scale_fill_manual(values = c("Actual Run" = "lightblue", "Older Runs" = "lightgreen"))

# Save plot
png(filename = output_file, width = 2800, height = 1800, res = 300)
print(boxplot)
dev.off()

# Filter data on padj < 0.05
filtered_outrider <- outrider %>% filter(padjust < 0.05)

# Summarize counts for each sample (Filtered Data)
filtered_count_data <- filtered_outrider %>%
  group_by(sampleID, RunCategory) %>%
  summarise(count = n(), .groups = 'drop')

# Identify outliers for both categories (Filtered Data)
filtered_outliers <- filtered_count_data %>%
  group_by(RunCategory) %>%
  filter(count %in% boxplot.stats(count)$out) %>% 
  ungroup()

# Create boxplot with outlier annotations (Filtered Data)
boxplot_filt <- ggplot(filtered_count_data, aes(x = RunCategory, y = count, fill = RunCategory)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16) +
  geom_text_repel(data = filtered_outliers, aes(x = RunCategory, y = count, label = sampleID),
                  size = 3, color = "red", max.overlaps = 10) +  # Annotate outliers
  theme_minimal() +
  labs(title = "Boxplots of Events Counts in Outrider (Filtered Data: padj < 0.05)",
       x = "Sample Group",
       y = "Count") +
  scale_fill_manual(values = c("Actual Run" = "lightblue", "Older Runs" = "lightgreen"))

# Save the filtered data plot as PNG
png(filename = output_file_filter, width = 2800, height = 1800, res = 300)
print(boxplot_filt)  # Print the plot to the device
dev.off()  # Close the device and save the plot

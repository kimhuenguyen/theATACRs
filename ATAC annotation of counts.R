# Load necessary packages, install if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}

# List of required packages
required_packages <- c("GenomicRanges", "rtracklayer", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "tidyverse")

# Install missing packages
install_if_missing(required_packages)

# Load libraries
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(dplyr)

# Load ATAC-seq log counts data
ATAC_log_counts <- read.table(gzfile("GSE162515_ATACseq_logCounts.txt.gz"), header = TRUE)

# Explore dimensions of the dataframe
#dim(ATAC_log_counts)

# Remove useless column 'rowNames'
ATAC_log_counts <- ATAC_log_counts %>% select(-rowNames)

# Display the first few rows of the dataframe
head(ATAC_log_counts)

# Create GRanges object
gr <- GRanges(seqnames = ATAC_log_counts$Chr,
              ranges = IRanges(start = ATAC_log_counts$Start, end = ATAC_log_counts$End))

# Load the transcript database
txs <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Display sequence lengths from the transcript database
head(seqlengths(txs))

# Annotate peaks with nearest gene information
annotate_result <- annotatePeak(gr, TxDb = txs)

# Convert annotated results to a dataframe
annotate_df <- as.data.frame(annotate_result)

# Merge annotated results with ATAC_log_counts dataframe
merged_df <- cbind(ATAC_log_counts, annotate_df)

# Reshape the merged dataframe to tidy format
tidy_df <- merged_df %>%
  pivot_longer(cols = starts_with("F"), names_to = "Patient", values_to = "logCounts") %>%
  mutate(counts = 2^(logCounts)) %>% # Convert log counts to normal counts
  filter(!grepl("\\.", Patient)) %>% # Remove rows with patient names containing a dot
  select(Patient, Chr, Start, End, counts, everything()) # Reorder columns

# Display the first few rows of the tidy dataframe
head(tidy_df)


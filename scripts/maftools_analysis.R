# Load necessary libraries
library(maftools)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read MAF files
maf1 <- read.maf(maf="/path/to/file1.filtered.annotated.maf.gz")
maf2 <- read.maf(maf="/path/to/file2.filtered.annotated.maf.gz")
maf3 <- read.maf(maf="/path/to/file3.filtered.annotated.maf.gz")
maf4 <- read.maf(maf="/path/to/file4.filtered.annotated.maf.gz")
maf5 <- read.maf(maf="/path/to/file5.filtered.annotated.maf.gz")
maf6 <- read.maf(maf="/path/to/file6.filtered.annotated.maf.gz")

# Manually assign Tumor_Sample_Barcode to each MAF object
maf1@data$Tumor_Sample_Barcode <- "Sample1"
maf2@data$Tumor_Sample_Barcode <- "Sample2"
maf3@data$Tumor_Sample_Barcode <- "Sample3"
maf4@data$Tumor_Sample_Barcode <- "Sample4"
maf5@data$Tumor_Sample_Barcode <- "Sample5"
maf6@data$Tumor_Sample_Barcode <- "Sample6"

# Create a list of MAF data frames
maf_list <- list(maf1, maf2, maf3, maf4, maf5, maf6)

# Merge the MAF objects
merged_maf <- merge_mafs(mafs = maf_list, verbose = TRUE)

# Filter out any rows where Tumor_Sample_Barcode is "__UNKNOWN__"
merged_maf@data <- merged_maf@data[merged_maf@data$Tumor_Sample_Barcode != "__UNKNOWN__", ]

# Plot MAF summary
plotmafSummary(maf = merged_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes = TRUE)

# Plot barplot of top 20 mutations
mafbarplot(maf = merged_maf, n = 20)

# Oncoplot for specific genes
oncoplot(maf = merged_maf, 
         genes = c("VHL", "PBMR1", "SETD2", "BAP1", "SMARC1", "KDM5C", "KDM6A", "MTOR", 
                   "CDKN2A", "CDKN2B", "MET", "NRF2", "NF2", "FH", "FLCN", "TP53", "TERT", 
                   "PTEN", "PI3KCA", "MITF", "TFE3", "TFEB", "SDHA", "SDHB", "SDHC", "SDHD", 
                   "ARID1A", "TSC1", "TSC2", "RB1", "ERBB4", "PDH"), 
         showTumorSampleBarcodes = TRUE, 
         sampleOrder = c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"), 
         SampleNamefontSize = 1.5, 
         titleText = "RCC-Associated Genes", 
         showPct = FALSE)

# Compare mutation load with TCGA cohorts
mutload_comparison <- tcgaCompare(maf = merged_maf, cohortName = 'Cohort', logscale = TRUE, 
                                  tcga_cohorts = c('KIRC', 'KICH', 'KIRP', 'PRAD', 'SKCM'))

# List of Hugo Symbols to filter by
hugo_symbols <- c("VHL", "PBMR1", "SETD2", "BAP1", "SMARC1", "KDM5C", "KDM6A", "MTOR", 
                  "CDKN2A", "CDKN2B", "MET", "NRF2", "NF2", "FH", "FLCN", "TP53", "TERT", 
                  "PTEN", "PI3KCA", "MITF", "TFE3", "TFEB", "SDHA", "SDHB", "SDHC", "SDHD", 
                  "ARID1A", "TSC1", "TSC2", "RB1", "ERBB4", "PDH")

# Filter merged_maf@data for rows where Hugo_Symbol is in the list of symbols
filtered_genes <- merged_maf@data %>% dplyr::filter(Hugo_Symbol %in% hugo_symbols)

# View the filtered data
filtered_genes <- filtered_genes %>% dplyr::select(Hugo_Symbol, Variant_Classification, Variant_Type, 
                                                   Tumor_Sample_Barcode, Genome_Change, Protein_Change)
print(filtered_genes)

# Create a binary mutation matrix for heatmap
mutation_matrix <- filtered_genes %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  distinct() %>%  # Remove duplicate mutations of the same gene in the same sample
  mutate(Present = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = Present, values_fill = list(Present = 0))

# Convert to long format for ggplot
mutation_long <- mutation_matrix %>%
  pivot_longer(cols = -Tumor_Sample_Barcode, names_to = "Gene", values_to = "Mutation")

# Plot the heatmap
ggplot(mutation_long, aes(x = Tumor_Sample_Barcode, y = Gene, fill = factor(Mutation))) +
  geom_tile() +
  scale_fill_manual(values = c("lightblue", "red"), name = "Mutation", labels = c("Absent", "Present")) +
  theme_minimal() +
  labs(title = "RCC-Related Mutations", x = "Gene", y = "Tumor Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# List of Hugo Symbols to filter by (for EBV analysis)
EBV_hugo_symbols <- c("SOCS1", "CD58", "NOTCH2", "NOTCH1", "B2M", "FOX01", "MTOR", "KMT2C", "ARID1A", "TP53")

# Filter merged_maf@data for rows where Hugo_Symbol is in the list of EBV-related symbols
EBV_filtered_genes <- merged_maf@data %>% dplyr::filter(Hugo_Symbol %in% EBV_hugo_symbols)

# Lollipop plot for VHL gene
VHL_loll <- lollipopPlot(maf = merged_maf, gene = 'VHL', AACol = 'Protein_Change', 
                         showMutationRate = FALSE, labelPos = "all")

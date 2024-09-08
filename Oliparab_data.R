#=================================================================
#Oliparab dataset processing 
#=================================================================
# install the required packages 
# install.packages("nanostringr")
# devtools::install_github("Nanostring-Biostats/NanoStringNCTools", 
#                      force = TRUE, ref = "master")
# load packages 
library(NanoStringNCTools)
library(nanostringr)
library(affy)
library(GEOquery)
library(tidyverse)
library(limma)
library(AgiMicroRna)
library(ggthemes)
library(ggiraph)
library(pheatmap)

# get supplementary files
getGEOSuppFiles("GSE121682")

# untar files
untar("GSE121682/GSE121682_RAW.tar", exdir = 'GSE121682/')

list.files(system.file(package = "NanoStringNCTools"))

# set your file location
datadir <- system.file("exdata" , "GSE121682", 
                       package = "NanoStringNCTools")
# Specify the directory where your .RCC.gz files are located
datadir <- "GSE121682/"

# List all .RCC.gz files in the directory
rcc_files <- list.files(datadir, pattern = "\\.RCC\\.gz$", full.names = TRUE)

# Process each file using parse_counts and parse_attributes
results <- lapply(rcc_files, function(rcc_file) {
  counts <- parse_counts(rcc_file)
  attributes <- parse_attributes(rcc_file)
  
  list(counts = counts, attributes = attributes)
})



# Specify the directory where your .RCC.gz files are located
datadir <- "GSE121682/"

# List all .RCC.gz files in the directory
rcc_files <- list.files(datadir, pattern = "\\.RCC\\.gz$", full.names = TRUE)

# Initialize lists to store counts and metadata
all_counts <- list()
all_metadata <- list()

# Loop through each file and extract the counts and metadata
for (rcc_file in rcc_files) {
  # Extract expression counts
  counts <- parse_counts(rcc_file)
  counts1 <-  DataFrame(counts$Name , counts[[4]])
  # Extract metadata
  metadata <- parse_attributes(rcc_file)
  
  # Store them in the lists
  all_counts[[rcc_file]] <- counts1
  all_metadata[[rcc_file]] <- metadata
}

# Combine all counts into a single data frame
combined_counts <- do.call(cbind, all_counts)

# Combine all metadata into a single data frame
combined_metadata <- do.call(rbind, all_metadata)

# View the combined data
head(combined_counts)
head(combined_metadata)

combined_counts <- as.data.frame(combined_counts)
# Identify columns that do not end with '.Name'
cols_to_keep <- grep("\\.Name$", colnames(combined_counts), invert = TRUE)

# Subset the data frame to keep only the columns that do not end with '.Name'
combined_counts <- combined_counts[, cols_to_keep]

rownames(combined_counts) <- counts$Name

# Use gsub to extract only the sample ID (e.g., GSM3443299)
colnames(combined_counts) <- gsub(".*\\.(GSM\\d+)_.*", "\\1", colnames(combined_counts))
write.csv(combined_counts, "olaparib_full_matrix.csv")

# Load necessary library
library(readr)

# Read the CSV file
df <- combined_counts

# Determine the number of genes
num_genes <- nrow(df)  # Subtract 1 for the header row
nrow(df)
# Create a vector of metabolin numbers
metabolon_numbers <- paste0("M", rep(1:ceiling(num_genes / 10), each = 10))

# Trim the metabolin numbers to match the number of genes
metabolon_numbers <- metabolon_numbers[1:num_genes]

# Add a new column with metabolin numbers
df$Metabolon <- c(metabolon_numbers)  # Add NA for the header row

# Optionally, write the modified data frame to a new CSV file
write_csv(df, 'modified_olaparib_full_matrix.csv')



# ==========================================================================
# Meta data manipulations ==================================================
#===========================================================================
# Load necessary library
library(dplyr)
pheno.olaparib <- read_csv("olaparib_metadata - Sheet1.csv")
rownames(pheno.olaparib) <- pheno.olaparib$Accession

# ==========================================================================
# ===================== Normalization unsing DESeq2 ======================== 
#===========================================================================
library(DESeq2)

# making the rownames and column names identical
all(rownames(pheno.olaparib) %in% colnames(combined_counts))
all(rownames(pheno.olaparib) == colnames(combined_counts))

count_matrix <- round(combined_counts)
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = pheno.olaparib,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

# dds <- dds[rowSums(counts(dds) >= 10) >= 12,]
nrow(dds) #  453 genes


# perform variance stabilization
dds_norm <- varianceStabilizingTransformation(dds)

# get normalized counts
normalized.count.matrix <- assay(dds_norm) %>% 
  t()

# norm.counts <- norm.counts[tumor$ids,]
normalized.count.matrix <- t(normalized.count.matrix) %>% as.data.frame()

df <- normalized.count.matrix
# Determine the number of genes
num_genes <- nrow(df)  # Subtract 1 for the header row
nrow(df)
# Create a vector of metabolin numbers
metabolon_numbers <- paste0("M", rep(1:ceiling(num_genes / 10), each = 10))

# Trim the metabolin numbers to match the number of genes
metabolon_numbers <- metabolon_numbers[1:num_genes]

# Add a new column with metabolin numbers
df$Metabolon <- c(metabolon_numbers)  # Add NA for the header row

# Optionally, write the modified data frame to a new CSV file
write_csv(df2, 'modified_olaparib_full_matrix.csv')

df2 <- df
df2$Gene_names <- rownames(df2)

# write.csv(normalized.count.matrix, "normalized.mRNAcount.Olaparib.csv")

##================================================================
# Add  the metabolon number column for normalized.count.matrix 
##================================================================

list_D <- read.csv("/home/admin/Sunchul/List_D.csv")
# Get the intersected gene symbols
intersected.genes <- intersect(rownames(normalized.count.matrix), list_D$Gene_name)
length(intersected.genes)
# Get the indices of intersected genes in metabolin data frame
intersected.indices <- which(list_D$Gene_name %in% intersected.genes)
list_D <- list_D[intersected.indices,]
length(list_D$Gene_name)
library(dplyr)

#normalized.count.matrix<- normalized.count.matrix[which(rownames(normalized.count.matrix) %in% intersected.genes),]
normalized.count.matrix$metabolon_number <- df$Metabolon
normalized.count.matrix <- normalized.count.matrix %>% relocate(metabolon_number, .before = 1)

normalized.count.matrix$metabolon_number %>% table()
normalized.count.matrix <- normalized.count.matrix[-66,]
# M12 M15 M19 M35 M37 M46 M48 M49 M50 M57 M58  M6 M60 M80 
# 6   8   4   6   8  12   3   2   2   6   4   2   2   1 
#  

#-------------------------------------------------------------------------------
# 01 subset the first androgen treated samples 
#-------------------------------------------------------------------------------

treated  <- pheno.olaparib[which(pheno.olaparib$Title == "olaparib treatment"),]

treated.counts <-df[,treated$Accession]
treated.counts$metabolon_number <- df$Metabolon
treated.counts <- treated.counts  %>% relocate(metabolon_number, .before = 1)

untreated <- pheno.olaparib[which(pheno.olaparib$Title == "no treatment"),]

untreated.counts <- df[,untreated$Accession]
untreated.counts$metabolon_number <- df$Metabolon
untreated.counts <- untreated.counts  %>% relocate(metabolon_number, .before = 1)

##==============================================================================
# --------------------- one function foe correlations --------------------------
#-------------------------------------------------------------------------------
##==============================================================================
library(ggplot2)
library(reshape2)

generate_heatmap_pdf <- function(count_matrix, pdf_filename) {
  # Load necessary libraries
  library(ggplot2)
  library(reshape2)
  
  # Split the data by metabolin
  metabolin_data <- split(count_matrix, count_matrix$metabolon_number)
  
  # Function to calculate R^2 matrix for a single metabolin
  calculate_r_squared <- function(metabolin_df) {
    metabolin_df <- metabolin_df[, -1]
    cor_matrix <- cor(metabolin_df, use = "pairwise.complete.obs")
    r_squared_matrix <- cor_matrix^2
    return(r_squared_matrix)
  }
  
  # Calculate R^2 matrices for each metabolin
  r_squared_list <- lapply(metabolin_data, calculate_r_squared)
  
  # Function to plot R^2 heatmap
  plot_r_squared_heatmap <- function(r_squared_matrix, metabolin) {
    melted_r_squared <- melt(r_squared_matrix)
    melted_r_squared$is_diagonal <- melted_r_squared$Var1 == melted_r_squared$Var2
    
    ggplot(melted_r_squared, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = ifelse(is_diagonal, "Diagonal", 
                                  ifelse(value > 0.55, "Above 0.55", "Below 0.55"))), 
                color = "black", size = 0.5) +
      geom_text(aes(label = ifelse(value > 0.55, sprintf("%.2f", value), "")), 
                size = 4, color = "white") +
      scale_fill_manual(values = c("Above 0.55" = "#FF0000", 
                                   "Below 0.55" = "#FFFFFF",
                                   "Diagonal" = "#808080"), 
                        name = "R^2 Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            legend.position = "bottom",
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1)) +
      labs(title = paste("Metabolon", metabolin),
           x = NULL, y = NULL) +
      coord_fixed(ratio = 1)
  }
  
  # Generate plots for each metabolin and save to PDF
  pdf(pdf_filename, width = 15, height = 15)
  
  for (metabolin in names(r_squared_list)) {
    print(plot_r_squared_heatmap(r_squared_list[[metabolin]], metabolin))
  }
  
  dev.off()
}


generate_heatmap_pdf(untreated.counts, "Heatmap_Prostate_untreated_olaparip_6samples.pdf")
generate_heatmap_pdf(treated.counts, "Heatmap_Prostate_treated_olaparip_6samples.pdf")

##==============================================================================
# ----------------- correlation 6h untreated and 6h treated ------------------
# ------------------------------------------------------------------------------
##==============================================================================

treated6h<- pheno.olaparib[which(pheno.olaparib$Title == "olaparib treatment" & pheno.olaparib$`Time point`=="6h"),]

treated6h.counts <- normalized.count.matrix[,treated6h$Accession]
treated6h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
treated6h.counts <- treated6h.counts  %>% relocate(metabolon_number, .before = 1)

untreated6h <- pheno.olaparib[which(pheno.olaparib$Title == "no treatment" & pheno.olaparib$`Time point`=="6h"),]

untreated6h.counts <- normalized.count.matrix[,untreated6h$Accession]
untreated6h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
untreated6h.counts <- untreated6h.counts  %>% relocate(metabolon_number, .before = 1)


generate_heatmap_pdf(untreated6h.counts, "Heatmap_Prostate_untreated_olaparip_6hr.pdf")
generate_heatmap_pdf(treated6h.counts, "Heatmap_Prostate_treated_olaparip_6hr.pdf")


##==============================================================================
# ----------------- correlation 12h untreated and 12h treated ------------------
# ------------------------------------------------------------------------------
##==============================================================================

treated12h<- pheno.olaparib[which(pheno.olaparib$Title == "olaparib treatment" & pheno.olaparib$`Time point`=="12h"),]

treated12h.counts <- normalized.count.matrix[,treated12h$Accession]
treated12h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
treated12h.counts <- treated12h.counts  %>% relocate(metabolon_number, .before = 1)

untreated12h <- pheno.olaparib[which(pheno.olaparib$Title == "no treatment" & pheno.olaparib$`Time point`=="12h"),]

untreated12h.counts <- normalized.count.matrix[,untreated12h$Accession]
untreated12h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
untreated12h.counts <- untreated12h.counts  %>% relocate(metabolon_number, .before = 1)


generate_heatmap_pdf(untreated12h.counts, "Heatmap_Prostate_untreated_olaparip_12hr.pdf")
generate_heatmap_pdf(treated12h.counts, "Heatmap_Prostate_treated_olaparip_12hr.pdf")



##==============================================================================
# ----------------- correlation 18h untreated and 18h treated ------------------
# ------------------------------------------------------------------------------
##==============================================================================

treated12h<- pheno.olaparib[which(pheno.olaparib$Title == "olaparib treatment" & pheno.olaparib$`Time point`=="18h"),]

treated12h.counts <- normalized.count.matrix[,treated12h$Accession]
treated12h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
treated12h.counts <- treated12h.counts  %>% relocate(metabolon_number, .before = 1)

untreated12h <- pheno.olaparib[which(pheno.olaparib$Title == "no treatment" & pheno.olaparib$`Time point`=="18h"),]

untreated12h.counts <- normalized.count.matrix[,untreated12h$Accession]
untreated12h.counts$metabolon_number <- normalized.count.matrix$metabolon_number
untreated12h.counts <- untreated12h.counts  %>% relocate(metabolon_number, .before = 1)


generate_heatmap_pdf(untreated12h.counts, "Heatmap_Prostate_untreated_olaparip_18hr.pdf")
generate_heatmap_pdf(treated12h.counts, "Heatmap_Prostate_treated_olaparip_18hr.pdf")

















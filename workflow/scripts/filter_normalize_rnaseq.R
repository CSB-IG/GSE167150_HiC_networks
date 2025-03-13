library(SummarizedExperiment)
library(DESeq2)
library(GenomicRanges)
library(tidyverse)

# se list has diffferent summarized experiments. We want the one with all
se.list <- readRDS("results/merged_se/se_list.Rds")
t <- se.list$all

# filter and normalize 
# this is how you get the counts
# figure out which samples we have
# filter and normalize (use the function from lymphoma MRA)
filter_low_count_genes_SE <- function(se_object, threshold = 10, 
                                      proportion = 0.8, assay_name = "unstranded") {
  # Extract count matrix from SE object
  count_matrix <- assay(se_object, assay_name)
  
  # Calculate the minimum number of samples that must exceed the threshold
  min_samples <- floor(ncol(count_matrix) * proportion)
  
  # Filter the genes based on the given threshold and proportion
  keep <- rowSums(count_matrix <= threshold) < min_samples
  
  # Subset the SE object using the logical vector
  filtered_se <- se_object[keep, ]
  
  # Return the filtered SE object
  return(filtered_se)
}

# filter lowly expressed genes
t <- filter_low_count_genes_SE(t)

# variance stabilizing transformation to the filtered matrix for normalization.
perform_vst_SE <- function(se_object, assay_name = "unstranded") {
  # Extract count matrix
  counts <- assay(se_object, assay_name)
  
  # create DESeqDataSet object
  # we need a minimal design formula, even if we don't use it
  dds <- DESeqDataSet(se_object, design = ~ 1)
  
  vst_data <- vst(dds)
  
  # Create new SE object with transformed counts
  # This preserves the original structure but updates the counts
  assay(se_object, "vst") <- assay(vst_data)
  
  return(se_object)
}

# perform vst
t <- perform_vst_SE(t)

# load gene annotation to annotate row data 
tcga_annotation <- readRDS("~/GSE167150_hiedge.git/results/igraph_intra/tcga_annotation.Rds")

# First, clean up ENSEMBL IDs by removing version numbers
# Function to remove version numbers from ENSEMBL IDs
clean_ensembl <- function(ids) {
  gsub("\\.[0-9]+$", "", ids)
}

# Clean rownames and gene_id column in rowData
rownames(t) <- clean_ensembl(rownames(t))
rowData(t)$gene_id <- clean_ensembl(rowData(t)$gene_id)
rownames(t@assays@data$unstranded) <- clean_ensembl(rownames(t@assays@data$unstranded))
rownames(t@assays@data$vst) <- rownames(t@assays@data$unstranded)

# Check for duplicates
rownames(t)[duplicated(rownames(t))]

# because we removed everything but coding genes and miRNA/lncRNA from the
# tcga annotation for the hiedge networks, some RNA-seq assay features are not there
table(rowData(t)[!(rowData(t)$gene_id %in% tcga_annotation$gene_id), "gene_type"])

# Subset SummarizedExperiment object
t <- t[rowData(t)$gene_id %in% tcga_annotation$gene_id, ]

# ANNOTATE
# Convert tcga_annotation to a data frame for easier manipulation
annot_df <- as.data.frame(tcga_annotation)
annot_df$strand <- NULL
colnames(annot_df)[1] <- "chromosome"

# Create a new DataFrame with all annotation information
new_rowData <- DataFrame(rowData(t))
new_rowData$type <- NULL
new_rowData$score <- NULL
new_rowData$phase <- NULL

# Remove duplicated columns before merging
annot_df <- annot_df[, !colnames(annot_df) %in% colnames(new_rowData) |
                             colnames(annot_df) == "gene_id"]

# Merge annotation information
merged_rowdata <- merge(
  as.data.frame(new_rowData),
  annot_df,
  by.x = "gene_id",
  by.y = "gene_id",
  all.x = TRUE
)

rm(annot_df, new_rowData)
# Convert back to DataFrame and set rownames
merged_rowdata <- DataFrame(merged_rowdata)
rownames(merged_rowdata) <- merged_rowdata$gene_id

# Update the rowData of the SummarizedExperiment
rowData(t) <- merged_rowdata[rownames(t), ]

rm(merged_rowdata)

colnames(rowData(t))

# Check for any genes that didn't get annotations
# missing_annot <- is.na(rowData(t)$seqnames)
# if(any(missing_annot)) {
#   message("Number of genes without annotations: ", sum(missing_annot))
#   message("First few genes without annotations:")
#   print(head(rownames(t)[missing_annot]))
# }


######## PCA COLORS tcga brest cancer, tcga normal, SYMBOLS gse tnbc, gse normal
library(SummarizedExperiment)
library(ggplot2)
library(ggfortify) 
library(DESeq2)    

# Extract assay data
assay_unstranded <- assay(t, "unstranded")
assay_vst <- assay(t, "vst")

# Perform PCA using prcomp (transpose is needed to have samples as rows)
pca_unstranded <- prcomp(t(assay_unstranded), center = TRUE, scale. = TRUE)
pca_vst <- prcomp(t(assay_vst), center = TRUE, scale. = TRUE)

# Extract sample metadata
metadata <- as.data.frame(colData(t))

# Create PCA data frames for each assay
pca_unstranded_df <- as.data.frame(pca_unstranded$x)
pca_vst_df <- as.data.frame(pca_vst$x)

# Add sample metadata for coloring
pca_unstranded_df$sample <- rownames(pca_unstranded_df)
pca_vst_df$sample <- rownames(pca_vst_df)

# Merge with metadata
pca_unstranded_df <- merge(pca_unstranded_df, metadata, by.x = "sample", by.y = "row.names")
pca_vst_df <- merge(pca_vst_df, metadata, by.x = "sample", by.y = "row.names")


# Function to plot PCA
plot_pca <- function(pca_obj, metadata, title, color_var, shape_var = NULL) {
  # Get explained variance
  explained_var <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)
  
  # Convert PCA results into a data frame
  pca_df <- as.data.frame(pca_obj$x)
  pca_df$sample <- rownames(pca_df)
  
  # Merge with metadata
  pca_df <- merge(pca_df, metadata, by.x = "sample", by.y = "row.names")
  
  # Identify GSE samples for labeling
  pca_df$label <- ifelse(pca_df$dataset == "GSE", pca_df$sample, NA)
  
  # Create base plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = !!sym(color_var)))
  
  # Add shape aesthetic if provided
  if (!is.null(shape_var)) {
    p <- p + aes(shape = !!sym(shape_var))
  }
  
  # Add points
  p <- p + geom_point(size = 3, alpha = 0.7) +
    
    # Add labels only for GSE samples
    geom_text(aes(label = label), hjust = 1.2, vjust = 1.2, size = 3, na.rm = TRUE) +
    
    # Axis labels with variance explained
    labs(
      title = title,
      x = paste0("PC1 (", explained_var[1], "% variance)"),
      y = paste0("PC2 (", explained_var[2], "% variance)")
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}


plot_pca(pca_unstranded, metadata, "PCA of Unstranded Assay", "BRCA_Subtype_PAM50")
plot_pca(pca_vst, metadata, "PCA of VST Assay", "BRCA_Subtype_PAM50")


# filter to keep just tcga normal and tcga 
# export per chromosome matrices (normal and tnbc)
# Subset for Basal samples
t_basal_tcga <- t[, colData(t)$sample == "Basal"]

# Subset for Normal_TCGA samples
t_normal_tcga <- t[, colData(t)$sample == "Normal_TCGA"]


write_assay_by_chromosome <- function(se_obj, assay_name, output_dir, phenotype) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract assay matrix
  assay_data <- assay(se_obj, assay_name)
  
  # Extract chromosome information
  row_data <- as.data.frame(rowData(se_obj))
  
  # Check if chromosome exists
  if (!"chromosome" %in% colnames(row_data)) {
    stop("Column 'chromosome' not found in rowData.")
  }
  
  # Loop through unique chromosomes and write separate CSVs
  for (chr in unique(row_data$chromosome)) {
    # Subset rows for this chromosome
    chr_rows <- row_data$chromosome == chr
    
    # Extract relevant assay data
    chr_matrix <- assay_data[chr_rows, , drop = FALSE]
    
    # Convert to data frame and add row names as a column
    chr_df <- as.data.frame(chr_matrix)
    chr_df <- cbind(geneID = rownames(chr_df), chr_df)  # Rename rownames as "Gene_ID"
    
    # Define output file path
    file_path <- file.path(output_dir, paste0(phenotype, "_assay_", assay_name, "_", chr, ".csv"))
    
    # Write to CSV with row names handled
    write.csv(chr_df, file = file_path, row.names = FALSE, quote = FALSE)  # Row names already in "Gene_ID" column
    
    # Print message
    message("Saved: ", file_path)
  }
}


# Write out the assay matrices separately for each chromosome
write_assay_by_chromosome(t_basal_tcga, "vst", "results/aracne_input_se/cancer_tcga", "BASAL-TCGA")
write_assay_by_chromosome(t_normal_tcga, "vst", "results/aracne_input_se/normal_tcga", "NORMAL-TCGA")

saveRDS(t_basal_tcga, file = "results/aracne_input_se/BASAL-TCGA_SE.Rds")
saveRDS(t_normal_tcga, file = "results/aracne_input_se/NORMAL-TCGA_SE.Rds")

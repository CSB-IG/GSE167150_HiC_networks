###########################################################################
#
# Summary Statistics of TNBC and Normal Networks
# Resolution is 40kb
#
###########################################################################
# setwd("/efs/users/hreyes/GSE167150_hiedge.git")
#
#
# 1. Read in HiEdge output (csv of interactions)
# 2. clean up data table
# 3. create networks
# 4. read in gencode annotation
# 5. read in lambert TF annotation
# 5. add gencode node attributes
# 6. read in tad annotation
# 7. add tad node attributes
# 8. read in compartment annotation
# 9. add compartment node attribute
# 10. export networks to R object and to graphml format
#
#
# libraries
library(data.table)
library(dplyr)
library(tidyverse)
library(igraph)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(ggraph)     # For network visualization
#library(visNetwork) # For interactive visualizations

######### set up variables
phenotype="Normal"
#phenotype="TNBC"

# beggers cant be choosers lol
if(phenotype == "Normal") {
  file_path <- "results/hiedge_intra/normal/40000/normal_40000_verbose.csv"
} else if(phenotype == "TNBC") {
  file_path <- "results/hiedge_intra/tnbc/40000/tnbc_40000_verbose.csv"
} else {
  stop("No more phenotypes")
}

my_resolution = 40000

#my_significance = 0.01
my_significance = 0.001

# Ensure chromosomes are in the correct order
chromosomes <- factor(paste0("chr", c(seq(1:22), "X")))  

######### FUNCTIONS
# function to create intrachromosomal hic networks obtained with HiEdge 
#
create_hic_networks_function <- function(hic_data, chromosomes, resolution, qvalue_threshold = 0.01) {
  # Validate resolution parameter
  if(!is.numeric(resolution) || resolution <= 0) {
    stop("resolution must be a positive number representing bin size in base pairs")
  }
  
  # Initialize list to store networks
  networks <- list()
  
  # Convert tibble to dataframe if needed
  hic_df <- as.data.frame(hic_data)
  
  # Create network for each chromosome
  for(chr in chromosomes) {
    # Filter data for current chromosome and q-value
    chr_data <- hic_df %>%
      filter(chr_1 == chr & 
               chr_2 == chr &  # Only intrachromosomal
               q_value <= qvalue_threshold)  # Apply q-value threshold
    
    # Create edges dataframe
    edges <- data.frame(
      source = chr_data$idx_1,
      target = chr_data$idx_2,
      qvalue = chr_data$q_value,
      zscore = chr_data$zscore_chr,
      distance = chr_data$genomic_distance,
      count = chr_data$interaction_count
    )
    
    # Create nodes dataframe with unique bins
    nodes <- data.frame(
      id = sort(unique(c(chr_data$idx_1, chr_data$idx_2))),
      chr = chr,
      midpoint = NA,
      start = NA,
      end = NA
    )
    
    # Fill in positions
    nodes$midpoint <- sapply(nodes$id, function(idx) {
      if(idx %in% chr_data$idx_1) {
        return(chr_data$midpoint_1[chr_data$idx_1 == idx][1])
      } else {
        return(chr_data$midpoint_2[chr_data$idx_2 == idx][1])
      }
    })
    
    # Calculate start and end using provided resolution
    half_resolution <- resolution / 2
    nodes$start <- nodes$midpoint - half_resolution
    nodes$end <- nodes$midpoint + half_resolution
    
    # Create igraph object (unweighted)
    g <- graph_from_data_frame(
      d = edges,
      vertices = nodes,
      directed = FALSE
    )
    
    # Store additional attributes in the graph
    E(g)$zscore <- edges$zscore  # Retain z-score as an edge attribute
    E(g)$qvalue <- edges$qvalue
    E(g)$distance <- edges$distance
    E(g)$count <- edges$count
    
    # Store network
    networks[[chr]] <- g
  }
  
  return(networks)
}

# function to clean up the gencode annotation 
# Step 1: Process gene annotations
process_gene_annotations_function <- function(annotation_gr, chromosomes, gene_types = c("protein_coding", "lncRNA", "miRNA")) {
  # Filter for desired chromosomes and gene types
  filtered_gr <- annotation_gr[
    seqnames(annotation_gr) %in% chromosomes & 
      annotation_gr$type == "gene" &
      annotation_gr$gene_type %in% gene_types
  ]
  
  # Remove version numbers from gene_id
  filtered_gr$gene_id <- gsub("\\.[0-9]+$", "", filtered_gr$gene_id)
  
  # Create a data frame for easier processing
  genes_df <- data.frame(
    chr = seqnames(filtered_gr),
    start = start(filtered_gr),
    end = end(filtered_gr),
    gene_name = filtered_gr$gene_name,
    gene_id = filtered_gr$gene_id,
    gene_type = filtered_gr$gene_type,
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates based on gene_name and coordinates
  genes_df <- genes_df %>%
    distinct(chr, start, end, gene_name, gene_id, gene_type)
  
  # Convert back to GRanges
  clean_gr <- GRanges(
    seqnames = genes_df$chr,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id,
    gene_type = genes_df$gene_type
  )
  
  return(clean_gr)
}

# function to add node annotations to networks (based on coordinates)
add_gene_annotations_function <- function(networks, gene_annotations) {
  # Validate inputs
  if (!is.list(networks) || !all(sapply(networks, inherits, "igraph"))) {
    stop("networks must be a named list of igraph objects.")
  }
  
  if (!inherits(gene_annotations, "GRanges")) {
    stop("gene_annotations must be a GRanges object.")
  }
  
  # Annotate each network
  annotated_networks <- lapply(names(networks), function(chr) {
    g <- networks[[chr]]
    
    # Check if network has necessary attributes
    if (!all(c("chr", "start", "end") %in% vertex_attr_names(g))) {
      stop(sprintf("Network for chromosome %s is missing 'chr', 'start', or 'end' node attributes.", chr))
    }
    
    # Get chromosome-specific annotations
    chr_genes <- gene_annotations[seqnames(gene_annotations) == chr]
    
    if (length(chr_genes) == 0) {
      warning(sprintf("No gene annotations found for chromosome %s", chr))
      return(g)
    }
    
    # Use only the start positions of the genes for overlap
    gene_starts <- GRanges(
      seqnames = seqnames(chr_genes),
      ranges = IRanges(start = start(chr_genes), end = start(chr_genes)),
      gene_name = chr_genes$gene_name,
      gene_id = chr_genes$gene_id,
      gene_type = chr_genes$gene_type
    )
    
    # Create GRanges for network nodes
    nodes_gr <- GRanges(
      seqnames = V(g)$chr,
      ranges = IRanges(start = V(g)$start, end = V(g)$end)
    )
    
    # Find overlaps (only the gene start positions are considered)
    overlaps <- findOverlaps(nodes_gr, gene_starts)
    
    # Initialize gene annotation attributes
    V(g)$genes <- NA_character_
    V(g)$gene_ids <- NA_character_
    V(g)$gene_types <- NA_character_
    
    # Process overlaps
    for (i in seq_along(nodes_gr)) {
      # Find all genes whose start position overlaps with this node
      gene_idx <- subjectHits(overlaps[queryHits(overlaps) == i])
      
      if (length(gene_idx) > 0) {
        V(g)$genes[i] <- paste(unique(gene_starts$gene_name[gene_idx]), collapse = ";")
        V(g)$gene_ids[i] <- paste(unique(gene_starts$gene_id[gene_idx]), collapse = ";")
        V(g)$gene_types[i] <- paste(unique(gene_starts$gene_type[gene_idx]), collapse = ";")
      }
    }
    
    return(g)
  })
  
  names(annotated_networks) <- names(networks)
  return(annotated_networks)
}

# add node type attribute to nodes
really_create_node_type <- function(graph) {
  # Get the gene types vector
  gene_types <- V(graph)$gene_types
  
  # Create node_type attribute
  V(graph)$node_type <- sapply(gene_types, function(gt) {
    # Check if gene_type is NA
    if (is.na(gt)) {
      return("N")
    }
    
    # Check if gene_type contains "protein_coding"
    if (grepl("protein_coding", gt)) {
      return("C")
    } else {
      return("R")
    }
  })
  
  return(graph)
}

# Function to process all networks in a list
add_ntype_function <- function(network_list) {
  # Apply the create_node_type function to each network in the list
  processed_networks <- lapply(network_list, really_create_node_type)
  return(processed_networks)
}

# function to add edge type
add_edge_types <- function(graph) {
  edge_list <- get.edgelist(graph, names=FALSE)  # Get numeric indices
  type1 <- V(graph)$node_type[edge_list[,1]]
  type2 <- V(graph)$node_type[edge_list[,2]]
  
  edge_types <- paste(pmin(type1, type2), pmax(type1, type2), sep="-")
  E(graph)$edge_type <- edge_types
  return(graph)
}


######### Load the data
hic_data <- fread(file_path)
# fread is for regular delimited files
# When you use fread from data.table package to read your CSV file, you get a data.table object
# A data.table is similar to a data.frame but optimized for faster data manipulation.

# Check the first few rows
head(hic_data)

# removeduplicate chr1 column
hic_data <- hic_data[, !duplicated(colnames(hic_data)), with = FALSE]
head(hic_data)

################ Threshold ################
# my only reason to keep 0.01 would be if I lose long range interactions
hic_data_01 <- hic_data[hic_data$q_value < 0.01, ]
hic_data_001 <- hic_data[hic_data$q_value < 0.001, ]

# Ensure the chromosome order matches the factor "chromosomes"
hic_data_01 <- hic_data_01 %>%
  mutate(chr_1 = factor(chr_1, levels = levels(chromosomes)))

hic_data_001 <- hic_data_001 %>%
  mutate(chr_1 = factor(chr_1, levels = levels(chromosomes)))

# Summarize data for counts and genomic_distance
summary_data_01 <- hic_data_01 %>%
  group_by(chr_1) %>%
  summarise(
    interaction_count = n(),
    median_distance = median(genomic_distance, na.rm = TRUE),
    .groups = "drop"
  )

summary_data_001 <- hic_data_001 %>%
  group_by(chr_1) %>%
  summarise(
    interaction_count = n(),
    median_distance = median(genomic_distance, na.rm = TRUE),
    .groups = "drop"
  )

# Combine genomic_distance for boxplots
combined_data <- bind_rows(
  hic_data_01 %>% mutate(dataset = "hic_data_01"),
  hic_data_001 %>% mutate(dataset = "hic_data_001")
)

# Boxplot of genomic_distance by chromosome
distances.boxp <- ggplot(combined_data, aes(x = factor(chr_1, levels = chromosomes),
                                            y = genomic_distance/1000, fill = dataset)) +
  geom_boxplot(alpha = 0.7, outliers = FALSE) +
  #scale_y_log10() +  # Log scale if needed for better visualization
  labs(
    title = paste("Comparison of Genomic Distance by Chromosome\n", phenotype),
    x = "Chromosome",
    y = "Genomic Distance (Kb)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Display counts
summary_counts <- left_join(
  summary_data_01,
  summary_data_001,
  by = "chr_1",
  suffix = c("_01", "_001")
)

summary_counts <- summary_counts[unlist(lapply(chromosomes, function(x) { which(summary_counts$chr_1 == x) })), ]

# plots
if(phenotype == "Normal") {
  pdf("results/hiedge_intra/normal/40000/significance_threshold_comparison.pdf", width = 13, height = 8)
} else if(phenotype == "TNBC") {
  pdf("results/hiedge_intra/tnbc/40000/significance_threshold_comparison.pdf", width = 13, height = 8)
} else {
  stop("No more phenotypes")
}
print(distances.boxp)

par(oma=c(1,2,1,1))
barplot(t(as.matrix(summary_counts[,c(2,4)]) ), las=1, names.arg = summary_counts$chr_1, beside=T, 
        main=paste("Comparison of Interactions number by Chromosome\n", phenotype))
dev.off()

#print(summary_counts)

rm(combined_data, distances.boxp, hic_data_001, hic_data_01, summary_counts, summary_data_001, summary_data_01)
############# Threshold ends ##############
###########################################

############## remove nonsignificant edges
#hic_data <- hic_data[hic_data$q_value < qval_thres, ]

# check interaction count
#quantile(hic_data$interaction_count, probs = seq(0,1,0.01))
#quantile(hic_data$genomic_distance, probs = seq(0,1,0.01))/1000000

#data <- hic_data
#data$interaction_count <- log2(data$interaction_count)  

############## calculate zscore
# interaction_count values have a heavy-tailed distribution
# Most of the values are very low but there are extreme outliers in the higher percentiles 
# This kind of distribution is common in Hi-C data due to 
## 1 Sparse contacts between many genomic regions
## 2 A few highly interacting regions
#
# The z-score is a statistical measure that indicates how many standard deviations 
# a data point is from the mean of the dataset
#
# zscores are sensitive to the underlying data distribution.
# Including non-significant interactions can skew the mean and standard deviation,
# reducing the ability of the z-scores to highlight meaningful deviations in significant interactions.
#
# we calculated the zscore on a chromosome by chromosome basis because each chromosome has
# its own interaction pattern due to differences in size, gene density, overall chromatin
# structure. 
# 
# calculate z scores within each chromosome
stopifnot(all(hic_data$chr_1 == hic_data$chr_2))

# group by chromosome, compute mean and sd of interaction count within each chr
# apply z score formula within groups
hic_data <- hic_data %>%
  group_by(chr_1) %>% # Group by the first chromosome column
  mutate(
    zscore_chr = (interaction_count - mean(interaction_count)) / sd(interaction_count) 
  ) %>%
  ungroup() # Remove grouping after calculation

############## create networks
intra_graphs <- create_hic_networks_function(hic_data,
                                             chromosomes,
                                             qvalue_threshold = my_significance,
                                             resolution = my_resolution)

############## get and process annotation 
# get the genes if they dont exist 
tcga_annotation <- readGFFAsGRanges("resources/tcga_gencode/gencode.v36.annotation.gtf")

tcga_annotation <- process_gene_annotations_function(tcga_annotation,
                                                     chromosomes = chromosomes)

############## annotate nodes with the gene data
# # Add annotations to networks
annotated_networks <- add_gene_annotations_function(intra_graphs, 
                                                    gene_annotations = tcga_annotation)


# save annotation
saveRDS(tcga_annotation, file = "results/igraph_intra/tcga_annotation.Rds")

# # Check results for one chromosome
# head(V(annotated_networks$chr1)$genes)
# V(annotated_networks$chr1)$gene_types
#
# empty strings are when there arent genes in the node

############## annotate nodes with node type
annotated_networks <- add_ntype_function(annotated_networks)

############## annotate edges with interaction type
annotated_networks <- lapply(annotated_networks, add_edge_types)


############## save object 
if(phenotype == "Normal") {
  saveRDS(hic_data, file = "results/igraph_intra/normal/40000/normal_intra_data_q001.Rds")
  saveRDS(annotated_networks, file = "results/igraph_intra/normal/40000/normal_intra_graphs_q001_genes.Rds")
} else if(phenotype == "TNBC") {
  saveRDS(hic_data, file = "results/igraph_intra/tnbc/40000/tnbc_intra_data_q001.Rds")
  saveRDS(annotated_networks, file = "results/igraph_intra/tnbc/40000/tnbc_intra_graphs_q001_genes.Rds")
} else {
  stop("No more phenotypes")
}



# #lyt <- layout_with_lgl(chromosome_graphs$chr22)
# plot(chromosome_graphs$chr22, layout = layout.fruchterman.reingold, vertex.label=NA,
#      main = "Chromosome 22 Normal")
# 
# 
# plot(chromosome_graphs$chr1, layout = layout.fruchterman.reingold)


# export for cytoscape
# Export an igraph object to a GraphML file
# write_graph(chromosome_graphs[["chr19"]], file = "chr19_network.graphml", format = "graphml")

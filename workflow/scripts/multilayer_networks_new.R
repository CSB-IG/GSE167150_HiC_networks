library(igraph)
library(muxViz) 
library(ggplot2)
library(tidyverse)
library(GenomicRanges)

normal_path ="results/igraph_intra/normal/40000/normal_intra_graphs_q001_genes.Rds"
tnbc_path="results/igraph_intra/tnbc/40000/tnbc_intra_graphs_q001_genes.Rds"
normal_graphs <- readRDS(normal_path)
tnbc_graphs <- readRDS(tnbc_path)

# read in matrices
MI_normal_path <- "results/aracne_networks/networks_normal_tcga/"

MI_tnbc_path <- "results/aracne_networks/networks_cancer_tcga/"

chromosomes <- paste0("chr", c(seq(1:22), "X"))  

##### Read normal mutual information/
lapply(c(seq(1, 22), "X"), function(c) { 
  
  file.path(MI_normal_path, 
            paste0("chromosome", c, "_network.txt")
  ) %>%
    lapply(read.table, header=TRUE, check.names=FALSE) -> df
  
  data.frame(df, Chromosome=paste0("chr", c)) -> df
  
  df$Chromosome <- factor(df$Chromosome, levels=chromosomes)
  colnames(df)[1:2] <- c("gene1", "gene2")
  
  return(df)
  
}) -> mi_normal
names(mi_normal) <- chromosomes

##### Read tnbc mutual information
lapply(c(seq(1, 22), "X"), function(c) { 
  
  file.path(MI_tnbc_path, 
            paste0("cancer_chromosome", c, "_network.txt")
  ) %>%
    lapply(read.table, header=TRUE, check.names=FALSE) -> df
  
  data.frame(df, Chromosome=paste0("chr", c)) -> df
  df$Chromosome <- factor(df$Chromosome, levels=chromosomes)
  colnames(df)[1:2] <- c("gene1", "gene2")
  
  return(df)
  
}) -> mi_tnbc
names(mi_tnbc) <- chromosomes


# load annotation 
tcga_annotation <- readRDS("results/igraph_intra/tcga_annotation.Rds")

################# Annotate MI networks
# FUNCTION

lookup_gene_positions <- function(mi_table, annotation) {
  # Convert annotation to a data frame for easier manipulation
  annotation_df <- as.data.frame(mcols(annotation))
  annotation_df$start <- start(annotation)
  annotation_df$end <- end(annotation)
  annotation_df$chromosome <- as.character(seqnames(annotation))
  
  # Rename columns for clarity
  colnames(annotation_df)[colnames(annotation_df) == "gene_id"] <- "gene"
  colnames(annotation_df)[colnames(annotation_df) == "gene_name"] <- "gene_name"
  
  # Merge mi_table with annotation to get genomic positions and names of gene1
  mi_table <- merge(mi_table, annotation_df[, c("gene", "start", "end", "chromosome", "gene_name")], 
                    by.x = "gene1", by.y = "gene", all.x = TRUE)
  colnames(mi_table)[(ncol(mi_table)-3):(ncol(mi_table))] <- c("start_gene1", "end_gene1", "chromosome_gene1", "gene_name1")
  
  # Merge mi_table with annotation to get genomic positions and names of gene2
  mi_table <- merge(mi_table, annotation_df[, c("gene", "start", "end", "chromosome", "gene_name")], 
                    by.x = "gene2", by.y = "gene", all.x = TRUE)
  colnames(mi_table)[(ncol(mi_table)-3):(ncol(mi_table))] <- c("start_gene2", "end_gene2", "chromosome_gene2", "gene_name2")
  
  # Set chromosome order
  mi_table$chromosome_gene1 <- factor(mi_table$chromosome_gene1, levels = paste0("chr", c(1:22, "X")))
  mi_table$chromosome_gene2 <- factor(mi_table$chromosome_gene2, levels = paste0("chr", c(1:22, "X")))
  
  return(mi_table)
}

lapply(chromosomes, function(c) {
  lookup_gene_positions(mi_normal[[c]], tcga_annotation)
}) -> mi_normal

names(mi_normal) <- chromosomes

lapply(chromosomes, function(c) {
  lookup_gene_positions(mi_tnbc[[c]], tcga_annotation)
}) -> mi_tnbc

names(mi_tnbc) <- chromosomes
################# 
#' Create a multilayer network from HiC and gene expression (MI) data using MuxViz
#'
#' @param hic_graph An igraph object representing the HiC network
#' @param mi_data A data frame with MI data (gene1, gene2, MI, etc.)
#' @return A list containing MuxViz-compatible multilayer network components
#' @importFrom igraph edge_attr as_edgelist vertex_attr V E
#' @importFrom GenomicRanges GRanges findOverlaps
create_muxviz_multilayer_network <- function(hic_graph, mi_data) {
  require(GenomicRanges)
  require(igraph)
  
  # Define layer names
  layers <- c("HiC", "MI")
  
  # ---- STEP 1: Map genes to HiC regions using GenomicRanges ----
  
  # Extract HiC node attributes for mapping
  hic_nodes <- data.frame(
    node_id = V(hic_graph)$name,
    chr = V(hic_graph)$chr,
    start = V(hic_graph)$start,
    end = V(hic_graph)$end,
    genes = V(hic_graph)$genes,
    gene_types = V(hic_graph)$gene_types,
    gene_ids = V(hic_graph)$gene_ids,
    node_type = V(hic_graph)$node_type,
    stringsAsFactors = FALSE
  )
  
  # Create GRanges for HiC nodes
  hic_gr <- GRanges(
    seqnames = hic_nodes$chr,
    ranges = IRanges(start = hic_nodes$start, end = hic_nodes$end),
    node_id = hic_nodes$node_id
  )
  
  # Extract gene coordinates from MI data
  # Create a data frame of all unique genes with their positions
  all_genes_info <- rbind(
    data.frame(gene = mi_data$gene_name1, chr = mi_data$Chromosome, start = mi_data$start_gene1, end = mi_data$end_gene1),
    data.frame(gene = mi_data$gene_name2, chr = mi_data$Chromosome, start = mi_data$start_gene2, end = mi_data$end_gene2)
  )
  
  # Remove duplicates to get unique gene entries
  all_genes_info <- all_genes_info[!duplicated(all_genes_info$gene), ]
  
  # Process genes from MI data with proper genomic coordinates
  gene1_gr <- GRanges(
    seqnames = mi_data$Chromosome,
    ranges = IRanges(start = mi_data$start_gene1, width = 1),
    gene = mi_data$gene_name1
  )
  
  gene2_gr <- GRanges(
    seqnames = mi_data$Chromosome,
    ranges = IRanges(start = mi_data$start_gene2, width = 1),
    gene = mi_data$gene_name2
  )
  
  # Combine all genes and find unique entries
  all_genes_gr <- c(gene1_gr, gene2_gr)
  all_genes_gr <- unique(all_genes_gr)
  
  # Find overlaps between genes and HiC regions
  overlaps <- findOverlaps(all_genes_gr, hic_gr)
  
  # Create mapping dictionary: gene to HiC node
  gene_to_hic_mapping <- data.frame(
    gene = all_genes_gr$gene[queryHits(overlaps)],
    hic_node = hic_gr$node_id[subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )
  
  # Handle duplicates (if a gene maps to multiple HiC nodes)
  gene_to_hic_mapping <- gene_to_hic_mapping[!duplicated(gene_to_hic_mapping$gene), ]
  
  # ---- STEP 2: Prepare MuxViz-compatible network data ----
  
  # 1. Create layer-specific edge lists
  
  # HiC layer edges - capture all edge attributes
  hic_edge_attrs <- edge_attr(hic_graph)
  hic_edges <- as.data.frame(as_edgelist(hic_graph))
  names(hic_edges) <- c("from", "to")
  
  # Add all edge attributes to the data frame
  for (attr_name in names(hic_edge_attrs)) {
    hic_edges[[attr_name]] <- hic_edge_attrs[[attr_name]]
  }
  
  # Ensure these specific attributes exist
  if (!"zscore" %in% names(hic_edges)) hic_edges$zscore <- NA
  if (!"distance" %in% names(hic_edges)) hic_edges$distance <- NA
  if (!"count" %in% names(hic_edges)) hic_edges$count <- NA
  if (!"edge_type" %in% names(hic_edges)) hic_edges$edge_type <- NA
  
  # MI layer edges
  mi_edges <- data.frame(
    from = mi_data$gene_name1,
    to = mi_data$gene_name2,
    weight = mi_data$MI,
    stringsAsFactors = FALSE
  )
  
  # 2. Create interlayer edges (for MuxViz's coupling)
  interlayer_edges <- data.frame(
    node = gene_to_hic_mapping$gene,
    nodeLayer = "MI",
    connectedNode = gene_to_hic_mapping$hic_node,
    connectedLayer = "HiC",
    weight = 1,  # Default weight for interlayer connections
    stringsAsFactors = FALSE
  )
  
  # 3. Prepare node metadata
  
  # All unique nodes in HiC layer with all attributes
  hic_layer_nodes <- data.frame(
    node = V(hic_graph)$name,
    layer = "HiC",
    chr = V(hic_graph)$chr,
    start = V(hic_graph)$start,
    end = V(hic_graph)$end,
    genes = V(hic_graph)$genes,
    gene_types = V(hic_graph)$gene_types,
    gene_ids = V(hic_graph)$gene_ids,
    node_type = V(hic_graph)$node_type,
    stringsAsFactors = FALSE
  )
  
  # All unique nodes in MI layer with proper genomic coordinates
  all_genes <- unique(c(mi_data$gene_name1, mi_data$gene_name2))
  
  # Create MI layer nodes with proper start and end positions
  mi_layer_nodes <- merge(
    data.frame(node = all_genes, layer = "MI"),
    all_genes_info,
    by.x = "node", 
    by.y = "gene", 
    all.x = TRUE
  )
  
  # Ensure the structure matches the HiC nodes
  # Add empty columns for HiC-specific attributes
  mi_layer_nodes$genes <- NA
  mi_layer_nodes$gene_types <- NA
  mi_layer_nodes$gene_ids <- NA
  mi_layer_nodes$node_type <- "gene"
  
  # Reorder columns to match HiC layer nodes
  mi_layer_nodes <- mi_layer_nodes[, c("node", "layer", "chr", "start", "end", 
                                       "genes", "gene_types", "gene_ids", "node_type")]
  
  # Combine node metadata
  all_nodes_metadata <- rbind(
    cbind(hic_layer_nodes, type="HiC"),
    cbind(mi_layer_nodes, type="gene")
  )
  
  # 4. Prepare layer information
  layer_info <- data.frame(
    layerId = 1:2,
    layerName = layers,
    stringsAsFactors = FALSE
  )
  
  # ---- STEP 3: Format for MuxViz ----
  
  # Create edge lists in MuxViz format (layerID, node1, node2, weight)
  hic_edges_muxviz <- data.frame(
    layerID = 1,  # HiC layer
    node1 = hic_edges$from,
    node2 = hic_edges$to,
    weight = hic_edges$count,  # Use count as the primary weight
    zscore = hic_edges$zscore,
    distance = hic_edges$distance,
    edge_type = hic_edges$edge_type,
    stringsAsFactors = FALSE
  )
  
  mi_edges_muxviz <- data.frame(
    layerID = 2,  # MI layer
    node1 = mi_edges$from,
    node2 = mi_edges$to,
    weight = mi_edges$weight,
    zscore = NA,
    distance = NA,
    edge_type = "MI",
    stringsAsFactors = FALSE
  )
  
  # Combine edges
  all_edges_muxviz <- rbind(hic_edges_muxviz, mi_edges_muxviz)
  
  # Return MuxViz components
  result <- list(
    edges = all_edges_muxviz,
    nodes = all_nodes_metadata,
    interlayer_edges = interlayer_edges,
    layer_info = layer_info,
    gene_to_hic_mapping = gene_to_hic_mapping
  )
  
  return(result)
}

#' Apply MuxViz multilayer network creation to all chromosomes
#'
#' @param hic_graphs List of igraph objects for HiC data by chromosome
#' @param mi_data_list List of MI data frames by chromosome
#' @return List of MuxViz-compatible multilayer networks by chromosome
create_all_muxviz_networks <- function(hic_graphs, mi_data_list) {
  # Check that both lists have the same chromosomes
  common_chrs <- intersect(names(hic_graphs), names(mi_data_list))
  
  # Create a list to store results
  muxviz_networks <- list()
  
  # Process each chromosome
  for (chr in common_chrs) {
    cat("Processing", chr, "...\n")
    tryCatch({
      muxviz_networks[[chr]] <- create_muxviz_multilayer_network(
        hic_graphs[[chr]], 
        mi_data_list[[chr]]
      )
      cat("Successfully processed", chr, "\n")
    }, error = function(e) {
      cat("Error processing", chr, ":", conditionMessage(e), "\n")
    })
  }
  
  return(muxviz_networks)
}

#' Write MuxViz output files for a multilayer network
#'
#' @param muxviz_net A MuxViz-compatible multilayer network list
#' @param output_dir Directory to save the files
#' @param prefix Prefix for output filenames
#' @param write_extended_attributes Logical, whether to write extended edge attributes
write_muxviz_files <- function(muxviz_net, output_dir, prefix = "multilayer", write_extended_attributes = TRUE) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define file paths
  edge_file <- file.path(output_dir, paste0(prefix, "_edges.txt"))
  extended_edge_file <- file.path(output_dir, paste0(prefix, "_edges_extended.txt"))
  layer_file <- file.path(output_dir, paste0(prefix, "_layers.txt"))
  coupling_file <- file.path(output_dir, paste0(prefix, "_couplings.txt"))
  node_metadata_file <- file.path(output_dir, paste0(prefix, "_nodes_metadata.txt"))
  
  # Write basic edge list for MuxViz standard format if it doesn't exist
  if (!file.exists(edge_file)) {
    write.table(
      muxviz_net$edges[, c("layerID", "node1", "node2", "weight")],
      file = edge_file,
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " "
    )
  }
  
  # Write extended edge attributes if requested and if the file doesn't exist
  if (write_extended_attributes && !file.exists(extended_edge_file)) {
    write.table(
      muxviz_net$edges,
      file = extended_edge_file,
      row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
    )
  }
  
  # Write layer information if it doesn't exist
  if (!file.exists(layer_file)) {
    write.table(
      muxviz_net$layer_info,
      file = layer_file,
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " "
    )
  }
  
  # Write interlayer edges (coupling) if it doesn't exist
  if (!file.exists(coupling_file)) {
    write.table(
      muxviz_net$interlayer_edges[, c("node", "nodeLayer", "connectedNode", "connectedLayer", "weight")],
      file = coupling_file,
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " "
    )
  }
  
  # Write node metadata if it doesn't exist
  if (!file.exists(node_metadata_file)) {
    write.table(
      muxviz_net$nodes,
      file = node_metadata_file,
      row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
    )
  }
  
  return(invisible(output_dir))
}


# Create MuxViz-compatible networks
normal_muxviz_networks <- create_all_muxviz_networks(normal_graphs, mi_normal)

tnbc_muxviz_networks <- create_all_muxviz_networks(tnbc_graphs, mi_tnbc)

# For a single chromosome
#chr1_muxviz <- create_muxviz_multilayer_network(normal_graphs$chr1, mi_normal$chr1)

# Export files for MuxViz GUI
lapply(chromosomes, function (c) {
  write_muxviz_files(normal_muxviz_networks[[c]], "results/multilayer/normal", 
                     paste0(c, "_multilayer_normal_40kb"))
})


lapply(chromosomes, function (c) {
  write_muxviz_files(tnbc_muxviz_networks[[c]], "results/multilayer/tnbc", 
                     paste0(c, "_multilayer_tnbc_40kb"))
})

###################################################################################
###################################################################################
# ---------------- STRUCTURE-EXPRESSION CORRELATION ANALYSIS ----------------

analyze_structure_expression_correlation <- function(muxviz_network) {
  # Convert to data.table for faster processing
  require(data.table)
  
  # Early validation
  if (is.null(muxviz_network) || !all(c("edges", "nodes", "gene_to_hic_mapping") %in% names(muxviz_network))) {
    return(list(paired_data = data.table(), correlation = NA))
  }
  
  # Extract HiC nodes and gene nodes from the nodes dataframe
  hic_nodes <- subset(muxviz_network$nodes, layer == "HiC")
  gene_nodes <- subset(muxviz_network$nodes, layer == "MI")
  
  # Create gene to HiC mapping from the provided mapping
  gene_to_hic <- data.table(muxviz_network$gene_to_hic_mapping)
  setkey(gene_to_hic, gene)  # Index for fast lookups
  
  # Convert edge data to data.tables for faster operations
  edges_dt <- data.table(muxviz_network$edges)
  hic_edges <- edges_dt[layerID == 1]  # Hi-C layer
  mi_edges <- edges_dt[layerID == 2]   # MI layer
  
  # If no MI edges, return NA
  if (nrow(mi_edges) == 0) {
    return(list(paired_data = data.table(), correlation = NA))
  }
  
  # Set keys for faster joins
  setkey(hic_edges, node1, node2)
  
  # Initialize empty result data.table
  paired_interactions <- data.table()
  
  # Process batches of edges for better memory management
  batch_size <- 1000
  n_batches <- ceiling(nrow(mi_edges) / batch_size)
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(mi_edges))
    
    batch_edges <- mi_edges[start_idx:end_idx]
    
    # For each MI edge, get the corresponding HiC regions
    # This avoids the problematic j expression
    temp_genes1 <- batch_edges$node1
    temp_genes2 <- batch_edges$node2
    temp_weights <- batch_edges$weight
    
    # Process all edges in the batch
    batch_results <- data.table()
    
    for (i in 1:length(temp_genes1)) {
      gene1 <- temp_genes1[i]
      gene2 <- temp_genes2[i]
      mi_weight <- temp_weights[i]
      
      # Get HiC regions for these genes
      hic_regions1 <- gene_to_hic[.(gene1), nomatch=0]
      hic_regions2 <- gene_to_hic[.(gene2), nomatch=0]
      
      if (nrow(hic_regions1) > 0 && nrow(hic_regions2) > 0) {
        # Create all combinations of HiC regions for this gene pair
        region_pairs <- CJ(hic_region1 = hic_regions1$hic_node,
                           hic_region2 = hic_regions2$hic_node)
        
        # Add gene pair info and MI strength
        region_pairs[, `:=`(
          gene1 = gene1,
          gene2 = gene2,
          mi_strength = mi_weight
        )]
        
        # Add to batch results
        batch_results <- rbindlist(list(batch_results, region_pairs), 
                                   use.names = TRUE, fill = TRUE)
      }
    }
    
    # Append batch results to overall results
    if (nrow(batch_results) > 0) {
      paired_interactions <- rbindlist(list(paired_interactions, batch_results), 
                                       use.names = TRUE, fill = TRUE)
    }
  }
  
  # If no paired interactions, return NA
  if (nrow(paired_interactions) == 0) {
    return(list(paired_data = data.table(), correlation = NA))
  }
  
  # Standardize region order to ensure proper matching with HiC edges
  paired_interactions[, `:=`(
    region1 = pmin(hic_region1, hic_region2),
    region2 = pmax(hic_region1, hic_region2)
  )]
  
  # Create lookup table for HiC edges with standardized node order
  hic_lookup <- copy(hic_edges)
  hic_lookup[, `:=`(
    region1 = pmin(node1, node2),
    region2 = pmax(node1, node2)
  )]
  setkey(hic_lookup, region1, region2)
  
  # Match with HiC weights
  setkey(paired_interactions, region1, region2)
  paired_interactions <- hic_lookup[paired_interactions, nomatch=0]
  
  # Calculate correlation
  if (nrow(paired_interactions) > 0) {
    correlation <- cor(paired_interactions$weight,
                       paired_interactions$mi_strength,
                       method = "spearman",
                       use = "pairwise.complete.obs")
  } else {
    correlation <- NA
  }
  
  return(list(paired_data = paired_interactions, correlation = correlation))
}

# ---------------- COMPARE CONDITIONS FUNCTION ----------------

compare_conditions <- function(normal_network, tnbc_network, chromosome) {
  require(data.table)
  
  # Analyze each condition
  normal_analysis <- analyze_structure_expression_correlation(normal_network)
  tnbc_analysis <- analyze_structure_expression_correlation(tnbc_network)
  
  # Early return if either analysis is empty
  if (is.na(normal_analysis$correlation) || is.na(tnbc_analysis$correlation)) {
    return(list(
      comparison_data = data.table(),
      normal_correlation = normal_analysis$correlation,
      tnbc_correlation = tnbc_analysis$correlation,
      change_correlation = NA,
      chromosome = chromosome
    ))
  }
  
  # Extract paired data
  normal_paired <- normal_analysis$paired_data
  tnbc_paired <- tnbc_analysis$paired_data
  
  # Convert to data.table if not already
  if (!is.data.table(normal_paired)) normal_paired <- as.data.table(normal_paired)
  if (!is.data.table(tnbc_paired)) tnbc_paired <- as.data.table(tnbc_paired)
  
  # Early return if either paired data is empty
  if (nrow(normal_paired) == 0 || nrow(tnbc_paired) == 0) {
    return(list(
      comparison_data = data.table(),
      normal_correlation = normal_analysis$correlation,
      tnbc_correlation = tnbc_analysis$correlation,
      change_correlation = NA,
      chromosome = chromosome
    ))
  }
  
  # Create pair identifiers for matching
  normal_paired[, pair_id := ifelse(gene1 < gene2, paste(gene1, gene2, sep = "_"), paste(gene2, gene1, sep = "_"))]
  tnbc_paired[, pair_id := ifelse(gene1 < gene2, paste(gene1, gene2, sep = "_"), paste(gene2, gene1, sep = "_"))]
  
  # Remove duplicate pairs
  normal_paired <- unique(normal_paired, by = "pair_id")
  tnbc_paired <- unique(tnbc_paired, by = "pair_id")
  
  # Set keys for joining
  setkey(normal_paired, pair_id)
  setkey(tnbc_paired, pair_id)
  
  # Find common gene pairs
  common_pairs <- normal_paired[tnbc_paired, nomatch=0, .(
    pair_id,
    hic_strength_normal = weight,
    mi_strength_normal = mi_strength,
    hic_strength_tnbc = i.weight,
    mi_strength_tnbc = i.mi_strength
  )]
  
  # Calculate changes
  common_pairs[, `:=`(
    hic_change = hic_strength_tnbc - hic_strength_normal,
    mi_change = mi_strength_tnbc - mi_strength_normal,
    hic_fold_change = log2((hic_strength_tnbc + 0.001) / (hic_strength_normal + 0.001)),
    mi_fold_change = log2((mi_strength_tnbc + 0.001) / (mi_strength_normal + 0.001))
  )]
  
  # Calculate correlation between changes
  change_correlation <- NA
  if (nrow(common_pairs) > 0) {
    change_correlation <- cor(common_pairs$hic_change, 
                              common_pairs$mi_change, 
                              method = "spearman",
                              use = "pairwise.complete.obs")
  }
  
  results <- list(
    comparison_data = common_pairs,
    normal_correlation = normal_analysis$correlation,
    tnbc_correlation = tnbc_analysis$correlation,
    change_correlation = change_correlation,
    chromosome = chromosome
  )
  
  return(results)
}

# ---------------- FIND SIGNIFICANT CHANGES FUNCTION ----------------

find_significant_changes <- function(comparison_results, hic_threshold = 1, mi_threshold = 0.5, output_file = NULL) {
  require(data.table)
  
  # Early validation
  if (is.null(comparison_results) || 
      !("comparison_data" %in% names(comparison_results)) ||
      nrow(comparison_results$comparison_data) == 0) {
    return(data.table())
  }
  
  # Convert to data.table if not already
  comparison_data <- comparison_results$comparison_data
  if (!is.data.table(comparison_data)) comparison_data <- as.data.table(comparison_data)
  
  # Filter for significant changes
  significant_changes <- comparison_data[
    abs(hic_fold_change) > hic_threshold & 
      abs(mi_fold_change) > mi_threshold
  ]
  
  # Add pattern classification
  if (nrow(significant_changes) > 0) {
    significant_changes[, change_pattern := 
                          ifelse(hic_fold_change > 0 & mi_fold_change > 0, "Both_Increased",
                                 ifelse(hic_fold_change < 0 & mi_fold_change < 0, "Both_Decreased",
                                        ifelse(hic_fold_change > 0 & mi_fold_change < 0, "HiC_Up_MI_Down",
                                               "HiC_Down_MI_Up")))
    ]
  }
  
  # Write results to file if output_file is provided
  if (!is.null(output_file) && nrow(significant_changes) > 0) {
    write.table(
      significant_changes,
      file = output_file,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
  }
  
  return(significant_changes)
}

# ---------------- MULTIPLEX COMMUNITY DETECTION FUNCTION ----------------

detect_multiplex_communities <- function(muxviz_network) {
  require(igraph)
  require(data.table)
  
  # Early validation
  if (is.null(muxviz_network) || 
      !all(c("edges", "nodes", "interlayer_edges") %in% names(muxviz_network))) {
    return(list(
      community_data = data.table(),
      community_stats = data.table(),
      igraph_communities = NULL
    ))
  }
  
  # Create edge tables as data.tables
  edges_dt <- as.data.table(muxviz_network$edges)
  
  # Check if edges exist for both layers
  if (nrow(edges_dt) == 0 || !all(c(1, 2) %in% edges_dt$layerID)) {
    return(list(
      community_data = data.table(),
      community_stats = data.table(),
      igraph_communities = NULL
    ))
  }
  
  # Extract edges by layer
  hic_edges <- edges_dt[layerID == 1, .(from = node1, to = node2, weight, layer = "HiC")]
  mi_edges <- edges_dt[layerID == 2, .(from = node1, to = node2, weight, layer = "MI")]
  
  # Get interlayer edges
  interlayer <- as.data.table(muxviz_network$interlayer_edges[, c("node", "connectedNode", "weight")])
  setnames(interlayer, c("from", "to", "weight"))
  interlayer[, layer := "inter"]
  
  # Combine all edges
  all_edges <- rbindlist(list(hic_edges, mi_edges, interlayer), use.names = TRUE, fill = TRUE)
  
  # Create unified graph
  tryCatch({
    unified_graph <- graph_from_data_frame(all_edges, directed = FALSE)
    
    # Apply Louvain community detection
    communities <- cluster_louvain(unified_graph)
    
    # Create node table with layer info
    nodes_dt <- data.table(
      node = V(unified_graph)$name
    )
    
    # Add layer information based on muxviz_network$nodes
    muxviz_nodes <- as.data.table(muxviz_network$nodes)
    
    # Set keys for joining
    setkey(muxviz_nodes, node)
    setkey(nodes_dt, node)
    
    # Add layer information
    nodes_dt <- muxviz_nodes[nodes_dt, .(node, layer)]
    
    # Add community assignments
    nodes_dt[, community := communities$membership]
    
    # Calculate community statistics
    community_stats <- nodes_dt[, .(
      total_nodes = .N,
      hic_nodes = sum(layer == "HiC"),
      mi_nodes = sum(layer == "MI")
    ), by = community]
    
    # Calculate layer mixing ratio
    community_stats[, mixing_ratio := pmin(
      hic_nodes / pmax(mi_nodes, 1),  # Avoid division by zero
      mi_nodes / pmax(hic_nodes, 1)
    )]
    
    return(list(
      community_data = nodes_dt,
      community_stats = community_stats,
      igraph_communities = communities
    ))
  }, error = function(e) {
    cat("Error in community detection:", conditionMessage(e), "\n")
    return(list(
      community_data = data.table(),
      community_stats = data.table(),
      igraph_communities = NULL
    ))
  })
}

# ---------------- COMPARE COMMUNITIES FUNCTION ----------------

compare_communities <- function(normal_communities, tnbc_communities, output_dir = NULL) {
  require(data.table)
  
  # Early validation
  if (is.null(normal_communities) || is.null(tnbc_communities) ||
      !("community_data" %in% names(normal_communities)) ||
      !("community_data" %in% names(tnbc_communities)) ||
      nrow(normal_communities$community_data) == 0 ||
      nrow(tnbc_communities$community_data) == 0) {
    return(list(
      node_changes = data.table(),
      summary = data.table(
        total_nodes = 0,
        changed_nodes = 0,
        percent_changed = 0,
        hic_nodes = 0,
        hic_changed = 0,
        hic_percent_changed = 0,
        mi_nodes = 0,
        mi_changed = 0,
        mi_percent_changed = 0
      )
    ))
  }
  
  # Extract community assignments
  normal_comm <- as.data.table(normal_communities$community_data)
  tnbc_comm <- as.data.table(tnbc_communities$community_data)
  
  # Set keys for joining
  setkey(normal_comm, node)
  setkey(tnbc_comm, node)
  
  # Find common nodes and join tables
  node_community_changes <- normal_comm[tnbc_comm, nomatch=0, .(
    node = node,
    normal_community = community,
    tnbc_community = i.community,
    layer = layer
  )]
  
  # Early return if no common nodes
  if (nrow(node_community_changes) == 0) {
    return(list(
      node_changes = data.table(),
      summary = data.table(
        total_nodes = 0,
        changed_nodes = 0,
        percent_changed = 0,
        hic_nodes = 0,
        hic_changed = 0,
        hic_percent_changed = 0,
        mi_nodes = 0,
        mi_changed = 0,
        mi_percent_changed = 0
      )
    ))
  }
  
  # Identify nodes that change communities
  node_community_changes[, changed := normal_community != tnbc_community]
  
  # Calculate overall statistics
  summary_stats <- node_community_changes[, .(
    total_nodes = .N,
    changed_nodes = sum(changed),
    percent_changed = 100 * sum(changed) / .N
  )]
  
  # Calculate statistics by layer
  layer_stats <- node_community_changes[, .(
    nodes = .N,
    changed = sum(changed),
    percent_changed = 100 * sum(changed) / .N
  ), by = layer]
  
  # Handle case where one layer might be missing
  hic_stats <- layer_stats[layer == "HiC"]
  mi_stats <- layer_stats[layer == "MI"]
  
  # Default values if a layer is missing
  if (nrow(hic_stats) == 0) hic_stats <- data.table(layer = "HiC", nodes = 0, changed = 0, percent_changed = 0)
  if (nrow(mi_stats) == 0) mi_stats <- data.table(layer = "MI", nodes = 0, changed = 0, percent_changed = 0)
  
  # Add layer-specific stats to summary
  summary_stats$hic_nodes <- hic_stats$nodes
  summary_stats$hic_changed <- hic_stats$changed
  summary_stats$hic_percent_changed <- hic_stats$percent_changed
  summary_stats$mi_nodes <- mi_stats$nodes
  summary_stats$mi_changed <- mi_stats$changed
  summary_stats$mi_percent_changed <- mi_stats$percent_changed
  
  # Write results to files if output_dir is provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    node_changes_file <- file.path(output_dir, "node_community_changes.csv")
    summary_file <- file.path(output_dir, "community_changes_summary.csv")
    
    fwrite(node_community_changes, node_changes_file)
    fwrite(summary_stats, summary_file)
  }
  
  return(list(
    node_changes = node_community_changes,
    summary = summary_stats
  ))
}

# ---------------- ANALYZE COMMUNITY CONTENTS FUNCTION ----------------

analyze_community_contents <- function(communities, gene_annotation) {
  require(data.table)
  
  # Early validation
  if (is.null(communities) || 
      !("community_data" %in% names(communities)) ||
      nrow(communities$community_data) == 0 ||
      is.null(gene_annotation)) {
    return(list(
      community_sizes = data.table(),
      community_analysis = data.table()
    ))
  }
  
  # Extract MI layer nodes (genes)
  gene_communities <- as.data.table(subset(communities$community_data, layer == "MI"))
  
  # Early return if no genes found
  if (nrow(gene_communities) == 0) {
    return(list(
      community_sizes = data.table(),
      community_analysis = data.table()
    ))
  }
  
  # Convert gene annotation to data.table
  gene_annotation <- as.data.table(gene_annotation)
  
  # Set keys for joining
  setkey(gene_communities, node)
  
  # Match column name for gene_id in annotation
  if ("gene_id" %in% names(gene_annotation)) {
    setkey(gene_annotation, gene_id)
    join_col <- "gene_id"
  } else if ("gene" %in% names(gene_annotation)) {
    setkey(gene_annotation, gene)
    join_col <- "gene"
  } else {
    # No matching column found
    return(list(
      community_sizes = data.table(community = unique(gene_communities$community), 
                                   gene_count = 0),
      community_analysis = data.table()
    ))
  }
  
  # Join with gene annotation
  gene_community_info <- gene_communities[gene_annotation, nomatch=0]
  
  # If no matches found, return empty results
  if (nrow(gene_community_info) == 0) {
    return(list(
      community_sizes = data.table(community = unique(gene_communities$community), 
                                   gene_count = 0),
      community_analysis = data.table()
    ))
  }
  
  # Get community sizes
  community_sizes <- gene_community_info[, .(gene_count = .N), by = community]
  
  # For each community, extract genes and analyze
  # Check if gene_biotype column exists
  has_biotype <- "gene_biotype" %in% names(gene_community_info)
  
  community_analysis <- gene_community_info[, {
    # Count by gene biotype if available
    biotype_counts <- NULL
    if (has_biotype) {
      biotype_counts <- .SD[, .N, by = gene_biotype]
    }
    
    # Get associated HiC regions for this community
    hic_regions <- communities$community_data[
      community == .BY$community & layer == "HiC", 
      node
    ]
    
    # Check if gene_name column exists
    gene_names <- NULL
    if ("gene_name" %in% names(.SD)) {
      gene_names <- .SD$gene_name
    } else {
      gene_names <- .SD$node  # Use node IDs if gene_name not available
    }
    
    list(
      community = .BY$community,
      gene_count = .N,
      genes = .SD$node,
      gene_names = gene_names,
      biotype_counts = biotype_counts,
      hic_regions = hic_regions,
      hic_region_count = length(hic_regions)
    )
  }, by = community]
  
  return(list(
    community_sizes = community_sizes,
    community_analysis = community_analysis
  ))
}

# ---------------- PLOTTING FUNCTIONS ----------------

# Plot structure-expression relationship
plot_structure_expression_relationship <- function(comparison, significant_changes = NULL) {
  require(ggplot2)
  require(ggrepel)
  
  # Early validation
  if (is.null(comparison) || 
      !("comparison_data" %in% names(comparison)) ||
      nrow(comparison$comparison_data) == 0) {
    # Return empty plot
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No data available for plotting") +
             theme_minimal() +
             theme(axis.title = element_blank()))
  }
  
  # Convert to data.table for easier handling
  plot_data <- as.data.table(comparison$comparison_data)
  
  # Check if significant changes are provided
  has_significant <- !is.null(significant_changes) && nrow(significant_changes) > 0
  
  # Prepare plot
  p <- ggplot(plot_data, aes(x = hic_fold_change, y = mi_fold_change)) +
    geom_point(alpha = 0.3, color = "darkgrey") +
    geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste("Structure-Expression Relationship -", comparison$chromosome),
      subtitle = paste("Correlation:", round(comparison$change_correlation, 3)),
      x = "Log2 Fold Change in Hi-C Interaction",
      y = "Log2 Fold Change in Gene Co-expression"
    )
  
  # Add significant points if available
  if (has_significant) {
    # Define colors for different patterns
    pattern_colors <- c(
      "Both_Increased" = "darkred",
      "Both_Decreased" = "darkblue",
      "HiC_Up_MI_Down" = "darkorange",
      "HiC_Down_MI_Up" = "darkgreen"
    )
    
    # Add significant points with color coding
    p <- p + geom_point(data = significant_changes, 
                        aes(color = change_pattern), 
                        size = 2, alpha = 0.8) +
      scale_color_manual(values = pattern_colors, name = "Change Pattern") +
      guides(color = guide_legend(title = "Change Pattern"))
    
    # Add labels for significant points using pair_id
    if ("pair_id" %in% names(significant_changes)) {
      significant_changes[, quadrant := ifelse(hic_fold_change > 0 & mi_fold_change > 0, "Q1",
                                               ifelse(hic_fold_change < 0 & mi_fold_change > 0, "Q2",
                                                      ifelse(hic_fold_change < 0 & mi_fold_change < 0, "Q3", "Q4")))]
      
      labeled_points <- significant_changes[, .SD[1:min(.N, 3)], by = quadrant]
      
      p <- p + geom_text_repel(data = labeled_points, 
                               aes(label = pair_id), 
                               size = 3, 
                               box.padding = 0.3, 
                               point.padding = 0.3, 
                               max.overlaps = Inf)
    }
  }
  
  return(p)
}

# Plot multiplex communities
plot_multiplex_communities <- function(muxviz_network, communities, title = NULL) {
  require(igraph)
  
  # Early validation
  if (is.null(communities) || 
      !("community_data" %in% names(communities)) ||
      is.null(communities$igraph_communities)) {
    # Create an empty plot
    plot.new()
    title(main = "No community data available")
    return(invisible(NULL))
  }
  
  # Extract community data
  community_data <- as.data.table(communities$community_data)
  
  # Create graph from edge data
  edges_dt <- as.data.table(muxviz_network$edges)
  
  # Combine layer edges
  hic_edges <- edges_dt[layerID == 1, .(from = node1, to = node2, weight, layer = "HiC")]
  mi_edges <- edges_dt[layerID == 2, .(from = node1, to = node2, weight, layer = "MI")]
  
  # Get interlayer edges
  interlayer <- as.data.table(muxviz_network$interlayer_edges[, c("node", "connectedNode", "weight")])
  setnames(interlayer, c("from", "to", "weight"))
  interlayer[, layer := "inter"]
  
  # Combine all edges
  all_edges <- rbindlist(list(hic_edges, mi_edges, interlayer), use.names = TRUE, fill = TRUE)
  
  # Create graph
  g <- graph_from_data_frame(all_edges, directed = FALSE)
  
  # Add community membership as vertex attribute
  community_lookup <- community_data[, .(node, community)]
  setkey(community_lookup, node)
  
  # Add vertex attributes
  V(g)$community <- communities$igraph_communities$membership[match(V(g)$name, community_data$node)]
  V(g)$layer <- community_data$layer[match(V(g)$name, community_data$node)]
  
  # Set node shapes and colors
  V(g)$shape <- ifelse(V(g)$layer == "HiC", "square", "circle")
  
  # Get a good color palette for communities
  num_communities <- length(unique(communities$igraph_communities$membership))
  community_colors <- rainbow(num_communities)
  
  # Set node colors by community
  V(g)$color <- community_colors[V(g)$community]
  
  # Set edge colors
  E(g)$color <- ifelse(all_edges$layer == "HiC", "red", 
                       ifelse(all_edges$layer == "MI", "blue", "gray"))
  
  # Set edge width based on weight
  E(g)$width <- scales::rescale(all_edges$weight, to = c(0.5, 3))
  
  # Plot the graph
  plot(g, 
       vertex.label = NA,
       vertex.size = ifelse(V(g)$layer == "HiC", 10, 6),
       edge.width = E(g)$width,
       layout = layout_with_fr(g),
       main = title
  )
  
  # Add legend for layers
  legend("topright", 
         legend = c("HiC Layer", "MI Layer", "Interlayer"),
         col = c("red", "blue", "gray"),
         lty = 1,
         lwd = 2,
         cex = 0.8)
  
  return(invisible(g))
}

# Plot aggregate results
plot_aggregate_results <- function(summary_results) {
  require(ggplot2)
  require(data.table)
  
  # Early validation
  if (is.null(summary_results) || 
      !all(c("correlations", "community_changes") %in% names(summary_results))) {
    # Return empty plots
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No summary data available") +
      theme_minimal() +
      theme(axis.title = element_blank())
    
    return(list(
      correlation_plot = empty_plot,
      community_change_plot = empty_plot
    ))
  }
  
  # Convert to data.table
  correlations <- as.data.table(summary_results$correlations)
  community_changes <- as.data.table(summary_results$community_changes)
  
  # Define genomic order
  genomic_order <- c(paste0("chr", 1:22), "chrX")
  
  # Set factor levels for chromosome column
  correlations[, chromosome := factor(chromosome, levels = genomic_order)]
  community_changes[, chromosome := factor(chromosome, levels = genomic_order)]
  
  # Plot 1: Correlations by chromosome
  correlation_plot <- ggplot(correlations) +
    geom_point(aes(x = chromosome, y = normal_correlation, color = "Normal"), size = 3) +
    geom_point(aes(x = chromosome, y = tnbc_correlation, color = "TNBC"), size = 3) +
    geom_point(aes(x = chromosome, y = change_correlation, color = "Change"), size = 3, shape = 17) +
    scale_color_manual(values = c("Normal" = "#008080", "TNBC" = "darkred", "Change" = "purple"),
                       name = "Condition") +
    labs(title = "Structure-Expression Correlations by Chromosome",
         x = "Chromosome",
         y = "Spearman Correlation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Community changes by chromosome
  community_change_plot <- ggplot(community_changes) +
    geom_bar(aes(x = chromosome, y = percent_changed), stat = "identity", 
             fill = "darkgrey", alpha = 0.7) +
    geom_point(aes(x = chromosome, y = hic_percent_changed, color = "HiC"), size = 3) +
    geom_point(aes(x = chromosome, y = mi_percent_changed, color = "MI"), size = 3) +
    scale_color_manual(values = c("HiC" = "red", "MI" = "blue"),
                       name = "Layer") +
    labs(title = "Community Changes by Chromosome",
         x = "Chromosome",
         y = "Percent Nodes Changed Community") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    correlation_plot = correlation_plot,
    community_change_plot = community_change_plot
  ))
}

# ---------------- MAIN WORKFLOW FUNCTION ----------------

run_multilayer_analysis <- function(chromosomes, normal_muxviz_networks, tnbc_muxviz_networks, gene_annotation = NULL) {
  require(data.table)
  
  # Create output directories if they don't exist
  dir.create("results/multilayer/figures", showWarnings = FALSE, recursive = TRUE)
  dir.create("results/multilayer/significant_changes", showWarnings = FALSE, recursive = TRUE)
  
  # Check which chromosomes have data in both conditions
  available_chrs <- chromosomes[chromosomes %in% names(normal_muxviz_networks) & 
                                  chromosomes %in% names(tnbc_muxviz_networks)]
  
  cat("Processing", length(available_chrs), "chromosomes:", paste(available_chrs, collapse = ", "), "\n")
  
  # Initialize results storage
  results_by_chromosome <- list()
  correlations <- data.table(chromosome = available_chrs)
  community_changes <- data.table(chromosome = available_chrs)
  
  # Process each chromosome
  for (i in 1:length(available_chrs)) {
    chr <- available_chrs[i]
    cat("Processing", chr, "...\n")
    
    # 1. Structure-expression relationship analysis
    comparison_results <- compare_conditions(
      normal_muxviz_networks[[chr]], 
      tnbc_muxviz_networks[[chr]],
      chr
    )
    
    # Store correlation results
    correlations[i, `:=`(
      normal_correlation = comparison_results$normal_correlation,
      tnbc_correlation = comparison_results$tnbc_correlation,
      change_correlation = comparison_results$change_correlation
    )]
    
    # 2. Find significant changes
    significant_changes_file <- paste0("results/multilayer/significant_changes/", chr, "_significant_changes.txt")
    significant_changes <- find_significant_changes(
      comparison_results,
      hic_threshold = 1,
      mi_threshold = 0.5,
      output_file = significant_changes_file
    )
    
    # 3. Multiplex community detection
    normal_communities <- detect_multiplex_communities(normal_muxviz_networks[[chr]])
    tnbc_communities <- detect_multiplex_communities(tnbc_muxviz_networks[[chr]])
    
    # 4. Compare communities
    community_comparison <- compare_communities(normal_communities, tnbc_communities, output_dir = "results/multilayer/community_comparison")
    
    # Store community change results
    if (nrow(community_comparison$summary) > 0) {
      community_changes[i, `:=`(
        total_nodes = community_comparison$summary$total_nodes,
        changed_nodes = community_comparison$summary$changed_nodes,
        percent_changed = community_comparison$summary$percent_changed,
        hic_nodes = community_comparison$summary$hic_nodes,
        hic_changed = community_comparison$summary$hic_changed,
        hic_percent_changed = community_comparison$summary$hic_percent_changed,
        mi_nodes = community_comparison$summary$mi_nodes,
        mi_changed = community_comparison$summary$mi_changed,
        mi_percent_changed = community_comparison$summary$mi_percent_changed
      )]
    }
    
    # 5. Analyze community contents
    normal_community_analysis <- NULL
    tnbc_community_analysis <- NULL
    
    if (!is.null(gene_annotation)) {
      normal_community_analysis <- analyze_community_contents(normal_communities, gene_annotation)
      tnbc_community_analysis <- analyze_community_contents(tnbc_communities, gene_annotation)
    }
    
    # 6. Visualizations
    plot_file_pdf <- paste0("results/multilayer/figures/", chr, "_structure_expression.pdf")
    plot_file_svg <- paste0("results/multilayer/figures/", chr, "_structure_expression.svg")
    pdf(plot_file_pdf, width = 10, height = 8)
    print(plot_structure_expression_relationship(comparison_results, significant_changes))
    dev.off()
    svg(plot_file_svg, width = 10, height = 8)
    print(plot_structure_expression_relationship(comparison_results, significant_changes))
    dev.off()
    
    plot_file_pdf <- paste0("results/multilayer/figures/", chr, "_communities_normal.pdf")
    plot_file_svg <- paste0("results/multilayer/figures/", chr, "_communities_normal.svg")
    pdf(plot_file_pdf, width = 10, height = 10)
    plot_multiplex_communities(
      normal_muxviz_networks[[chr]], 
      normal_communities, 
      title = paste("Normal -", chr)
    )
    dev.off()
    svg(plot_file_svg, width = 10, height = 10)
    plot_multiplex_communities(
      normal_muxviz_networks[[chr]], 
      normal_communities, 
      title = paste("Normal -", chr)
    )
    dev.off()
    
    plot_file_pdf <- paste0("results/multilayer/figures/", chr, "_communities_tnbc.pdf")
    plot_file_svg <- paste0("results/multilayer/figures/", chr, "_communities_tnbc.svg")
    pdf(plot_file_pdf, width = 10, height = 10)
    plot_multiplex_communities(
      tnbc_muxviz_networks[[chr]], 
      tnbc_communities, 
      title = paste("TNBC -", chr)
    )
    dev.off()
    svg(plot_file_svg, width = 10, height = 10)
    plot_multiplex_communities(
      tnbc_muxviz_networks[[chr]], 
      tnbc_communities, 
      title = paste("TNBC -", chr)
    )
    dev.off()
    
    # Store all results for this chromosome
    results_by_chromosome[[chr]] <- list(
      comparison = comparison_results,
      significant_changes = significant_changes,
      normal_communities = normal_communities,
      tnbc_communities = tnbc_communities,
      community_comparison = community_comparison,
      normal_community_analysis = normal_community_analysis,
      tnbc_community_analysis = tnbc_community_analysis
    )
    
    cat("Completed", chr, "\n")
  }
  
  # Aggregate results
  summary_results <- list(
    correlations = correlations,
    community_changes = community_changes
  )
  
  # Plot aggregate results
  aggregate_plots <- plot_aggregate_results(summary_results)
  
  plot_file_pdf <- paste0("results/multilayer/figures/aggregate_correlation_plot.pdf")
  plot_file_svg <- paste0("results/multilayer/figures/aggregate_correlation_plot.svg")
  pdf(plot_file_pdf, width = 12, height = 8)
  print(aggregate_plots$correlation_plot)
  dev.off()
  svg(plot_file_svg, width = 12, height = 8)
  print(aggregate_plots$correlation_plot)
  dev.off()
  
  plot_file_pdf <- paste0("results/multilayer/figures/aggregate_community_changes.pdf")
  plot_file_svg <- paste0("results/multilayer/figures/aggregate_community_changes.svg")
  pdf(plot_file_pdf, width = 12, height = 8)
  print(aggregate_plots$community_change_plot)
  dev.off()
  svg(plot_file_svg, width = 12, height = 8)
  print(aggregate_plots$community_change_plot)
  dev.off()
  
  # Return all results
  return(list(
    by_chromosome = results_by_chromosome,
    summary = summary_results,
    aggregate_plots = aggregate_plots
  ))
}

################
run_multilayer_analysis(chromosomes, normal_muxviz_networks, tnbc_muxviz_networks) -> multilayer_analysis_results_list

saveRDS(multilayer_analysis_results_list, file="results/multilayer/multilayer_analysis_all_chromosomes.Rds")
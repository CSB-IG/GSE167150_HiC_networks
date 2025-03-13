#!/usr/bin/env Rscript

# Complete Network Analysis Script
# This script performs comprehensive analysis of normal vs TNBC networks

library(igraph)
library(data.table)

# Set number of threads for data.table
setDTthreads(parallel::detectCores() - 3)

# 1. Basic network metrics function
get_basic_metrics <- function(normal_graph, tnbc_graph, chr) {
  dt <- data.table(
    chromosome = chr,
    normal_nodes_n = vcount(normal_graph),
    tnbc_nodes_n = vcount(tnbc_graph),
    normal_edges_n = ecount(normal_graph),
    tnbc_edges_n = ecount(tnbc_graph)
  )
  return(dt)
}

# 2. Node comparison metrics
get_node_comparison_metrics <- function(normal_graph, tnbc_graph, chr) {
  normal_nodes <- V(normal_graph)$name
  tnbc_nodes <- V(tnbc_graph)$name
  
  common_nodes <- intersect(normal_nodes, tnbc_nodes)
  union_nodes <- union(normal_nodes, tnbc_nodes)
  
  dt <- data.table(
    chromosome = chr,
    common_nodes_n = length(common_nodes),
    normal_unique_nodes_n = length(setdiff(normal_nodes, tnbc_nodes)),
    tnbc_unique_nodes_n = length(setdiff(tnbc_nodes, normal_nodes)),
    jaccard_nodes = length(common_nodes) / length(union_nodes)
  )
  return(dt)
}

# 3. Edge comparison metrics
get_edge_comparison_metrics <- function(normal_graph, tnbc_graph, chr) {
  normal_el <- as.data.table(get.edgelist(normal_graph))
  tnbc_el <- as.data.table(get.edgelist(tnbc_graph))
  
  setnames(normal_el, c("V1", "V2"))
  setnames(tnbc_el, c("V1", "V2"))
  
  normal_el[, edge_id := paste(pmin(V1, V2), pmax(V1, V2), sep="--")]
  tnbc_el[, edge_id := paste(pmin(V1, V2), pmax(V1, V2), sep="--")]
  
  common_edges <- normal_el[tnbc_el, on = .(edge_id)]
  all_edges <- unique(rbindlist(list(
    normal_el[, .(edge_id)],
    tnbc_el[, .(edge_id)]
  )))
  
  dt <- data.table(
    chromosome = chr,
    jaccard_edges = common_edges[, .N] / all_edges[, .N]
  )
  return(dt)
}

# 4. Network property metrics
get_network_properties <- function(normal_graph, tnbc_graph, chr) {
  dt <- data.table(
    chromosome = chr,
    normal_avg_degree = mean(degree(normal_graph)),
    tnbc_avg_degree = mean(degree(tnbc_graph)),
    normal_density = graph.density(normal_graph),
    tnbc_density = graph.density(tnbc_graph),
    normal_clustering = transitivity(normal_graph, type="global"),
    tnbc_clustering = transitivity(tnbc_graph, type="global"),
    normal_components = count_components(normal_graph),
    tnbc_components = count_components(tnbc_graph)
  )
  return(dt)
}

# 5. Gene metrics
get_gene_metrics <- function(normal_graph, tnbc_graph, chr) {
  normal_genes <- sapply(V(normal_graph)$genes, function(x) length(unlist(strsplit(x, ","))))
  tnbc_genes <- sapply(V(tnbc_graph)$genes, function(x) length(unlist(strsplit(x, ","))))
  
  dt <- data.table(
    chromosome = chr,
    normal_genes_total_n = sum(normal_genes),
    tnbc_genes_total_n = sum(tnbc_genes),
    normal_avg_genes_per_node = mean(normal_genes),
    tnbc_avg_genes_per_node = mean(tnbc_genes)
  )
  return(dt)
}

# 6. Complexity metrics
get_complexity_metrics <- function(normal_graph, tnbc_graph, chr) {
  # Wrapper for safe calculation
  safe_calc <- function(expr) {
    tryCatch(expr, error = function(e) NA_real_)
  }
  
  dt <- data.table(chromosome = chr)
  
  # Basic complexity
  dt[, `:=`(
    normal_max_clique_size = safe_calc(length(largest.cliques(normal_graph)[[1]])),
    tnbc_max_clique_size = safe_calc(length(largest.cliques(tnbc_graph)[[1]])),
    normal_maximal_cliques_n = safe_calc(count_max_cliques(normal_graph)),
    tnbc_maximal_cliques_n = safe_calc(count_max_cliques(tnbc_graph)),
    normal_diameter = safe_calc(diameter(normal_graph, directed=FALSE, unconnected=TRUE)),
    tnbc_diameter = safe_calc(diameter(tnbc_graph, directed=FALSE, unconnected=TRUE)),
    normal_avg_path_length = safe_calc(mean_distance(normal_graph, directed=FALSE, unconnected=TRUE)),
    tnbc_avg_path_length = safe_calc(mean_distance(tnbc_graph, directed=FALSE, unconnected=TRUE))
  )]
  
  # Centralization metrics
  dt[, `:=`(
    normal_degree_centralization = safe_calc(centralization.degree(normal_graph)$centralization),
    tnbc_degree_centralization = safe_calc(centralization.degree(tnbc_graph)$centralization),
    normal_betweenness_centralization = safe_calc(centralization.betweenness(normal_graph)$centralization),
    tnbc_betweenness_centralization = safe_calc(centralization.betweenness(tnbc_graph)$centralization)
  )]
  
  # Assortativity
  dt[, `:=`(
    normal_assortativity = safe_calc(assortativity.degree(normal_graph, directed=FALSE)),
    tnbc_assortativity = safe_calc(assortativity.degree(tnbc_graph, directed=FALSE))
  )]
  
  # Structural metrics
  dt[, `:=`(
    normal_girth = safe_calc(girth(normal_graph)$girth),
    tnbc_girth = safe_calc(girth(tnbc_graph)$girth),
    normal_edge_connectivity = safe_calc(edge.connectivity(normal_graph)),
    tnbc_edge_connectivity = safe_calc(edge.connectivity(tnbc_graph))
  )]
  
  return(dt)
}

# 7. Detailed node metrics
get_detailed_node_metrics <- function(normal_graph, tnbc_graph, chr) {
  # Get metrics for normal network nodes
  normal_nodes <- data.table(
    chromosome = chr,
    node_id = V(normal_graph)$name,
    network = "normal",
    degree = degree(normal_graph),
    betweenness = betweenness(normal_graph),
    closeness = closeness(normal_graph),
    clustering_coefficient = transitivity(normal_graph, type="local"),
    genes = V(normal_graph)$genes
  )
  
  # Get metrics for TNBC network nodes
  tnbc_nodes <- data.table(
    chromosome = chr,
    node_id = V(tnbc_graph)$name,
    network = "tnbc",
    degree = degree(tnbc_graph),
    betweenness = betweenness(tnbc_graph),
    closeness = closeness(tnbc_graph),
    clustering_coefficient = transitivity(tnbc_graph, type="local"),
    genes = V(tnbc_graph)$genes
  )
  
  # Combine and identify status
  all_nodes <- rbindlist(list(normal_nodes, tnbc_nodes))
  
  normal_nodes_set <- V(normal_graph)$name
  tnbc_nodes_set <- V(tnbc_graph)$name
  
  all_nodes[, node_status := case_when(
    node_id %in% intersect(normal_nodes_set, tnbc_nodes_set) ~ "common",
    network == "normal" ~ "normal_only",
    network == "tnbc" ~ "tnbc_only"
  )]
  
  return(all_nodes)
}

# 8. Detailed edge metrics
get_detailed_edge_metrics <- function(normal_graph, tnbc_graph, chr) {
  # Get edges for normal network
  normal_edges <- as.data.table(get.edgelist(normal_graph))
  setnames(normal_edges, c("node1", "node2"))
  normal_edges[, `:=`(
    chromosome = chr,
    network = "normal",
    edge_id = paste(pmin(node1, node2), pmax(node1, node2), sep="--"),
    zscore = E(normal_graph)$zscore,
    distance = E(normal_graph)$distance
  )]
  
  # Get edges for TNBC network
  tnbc_edges <- as.data.table(get.edgelist(tnbc_graph))
  setnames(tnbc_edges, c("node1", "node2"))
  tnbc_edges[, `:=`(
    chromosome = chr,
    network = "tnbc",
    edge_id = paste(pmin(node1, node2), pmax(node1, node2), sep="--"),
    zscore = E(tnbc_graph)$zscore,
    distance = E(tnbc_graph)$distance
  )]
  
  # Combine and identify edge status
  all_edges <- rbindlist(list(normal_edges, tnbc_edges))
  
  normal_edge_ids <- normal_edges$edge_id
  tnbc_edge_ids <- tnbc_edges$edge_id
  
  all_edges[, edge_status := case_when(
    edge_id %in% intersect(normal_edge_ids, tnbc_edge_ids) ~ "common",
    network == "normal" ~ "normal_only",
    network == "tnbc" ~ "tnbc_only"
  )]
  
  return(all_edges)
}

# Main analysis function for a single chromosome
analyze_chromosome <- function(normal_graph, tnbc_graph, chr, output_dir) {
  cat(sprintf("\nProcessing chromosome %s...\n", chr))
  
  # Create chromosome-specific directory
  chr_dir <- file.path(output_dir, paste0("chr", chr))
  dir.create(chr_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get all metrics
  cat("- Computing basic metrics...\n")
  basic <- get_basic_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing node comparison metrics...\n")
  nodes <- get_node_comparison_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing edge comparison metrics...\n")
  edges <- get_edge_comparison_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing network properties...\n")
  properties <- get_network_properties(normal_graph, tnbc_graph, chr)
  
  cat("- Computing gene metrics...\n")
  genes <- get_gene_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing complexity metrics...\n")
  complexity <- get_complexity_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing detailed node metrics...\n")
  detailed_nodes <- get_detailed_node_metrics(normal_graph, tnbc_graph, chr)
  
  cat("- Computing detailed edge metrics...\n")
  detailed_edges <- get_detailed_edge_metrics(normal_graph, tnbc_graph, chr)
  
  # Create summary table
  summary_metrics <- Reduce(function(x, y) merge(x, y, by="chromosome"), 
                            list(basic, nodes, edges, properties, genes, complexity))
  
  # Save all results
  cat("- Saving results...\n")
  fwrite(summary_metrics, file=file.path(chr_dir, "summary_metrics.csv"))
  fwrite(detailed_nodes, file=file.path(chr_dir, "detailed_node_metrics.csv"))
  fwrite(detailed_edges, file=file.path(chr_dir, "detailed_edge_metrics.csv"))
  
  cat(sprintf("Chromosome %s completed successfully\n", chr))
  
  return(list(
    summary = summary_metrics,
    detailed_nodes = detailed_nodes,
    detailed_edges = detailed_edges
  ))
}

# Main function to analyze all chromosomes
analyze_all_chromosomes <- function(normal_graphs, tnbc_graphs, output_dir = "network_analysis_results") {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Start time
  start_time <- Sys.time()
  cat(sprintf("\nStarting analysis at %s\n", start_time))
  
  # Process each chromosome
  results <- list()
  
  for (chr in names(normal_graphs)) {
    tryCatch({
      results[[chr]] <- analyze_chromosome(normal_graphs[[chr]], 
                                           tnbc_graphs[[chr]], 
                                           chr, 
                                           output_dir)
    }, error = function(e) {
      cat(sprintf("Error processing chromosome %s: %s\n", chr, e$message))
    })
  }
  
  # Combine all summary results
  cat("\nCombining results...\n")
  final_table <- rbindlist(lapply(results, function(x) x$summary))
  
  # Save combined results
  fwrite(final_table, file=file.path(output_dir, "all_chromosomes_metrics.csv"))
  
  # End time
  end_time <- Sys.time()
  time_taken <- difftime(end_time, start_time, units="mins")
  
  cat(sprintf("\nAnalysis completed at %s\n", end_time))
  cat(sprintf("Total time taken: %.2f minutes\n", time_taken))
  
  return(final_table)
}

# Main execution
main <- function() {
  cat("Starting network analysis...\n")
  
  # Load data
  normal_path ="results/igraph_intra/normal/40000/normal_intra_graphs_q001_genes.Rds"
  tnbc_path="results/igraph_intra/tnbc/40000/tnbc_intra_graphs_q001_genes.Rds"
  
  normal_graphs <- readRDS(normal_path)
  tnbc_graphs <- readRDS(tnbc_path)
  
  # Run analysis
  results <- analyze_all_chromosomes(normal_graphs, tnbc_graphs)
  
  cat("Analysis complete. Results are saved in the 'network_analysis_results' directory.\n")
}

# Run the main function
main()

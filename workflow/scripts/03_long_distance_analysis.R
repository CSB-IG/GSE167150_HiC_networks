library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)

analyze_long_range_coding_network <- function(graph, name = "chr") {
  # 1. Get the long-range coding subnetwork
  coding_edges <- E(graph)[E(graph)$edge_type %in% c("C-C", "C-N", "C-R")]
  coding_distances <- E(graph)$distance[E(graph)$edge_type %in% c("C-C", "C-N", "C-R")]
  distance_threshold <- quantile(coding_distances, 0.99)
  
  selected_edges <- which(
    E(graph)$edge_type %in% c("C-C", "C-N", "C-R") & 
      E(graph)$distance >= distance_threshold
  )
  
  long_range_net <- subgraph.edges(graph, selected_edges, delete.vertices = TRUE)
  
  # 2. Network Statistics
  network_stats <- list(
    name = name,
    threshold_distance = distance_threshold,
    num_nodes = vcount(long_range_net),
    num_edges = ecount(long_range_net),
    network_density = edge_density(long_range_net),
    mean_degree = mean(degree(long_range_net)),
    median_degree = median(degree(long_range_net)),
    clustering_coefficient = transitivity(long_range_net, type = "global"),
    components = components(long_range_net)$no,
    largest_component_size = max(components(long_range_net)$csize)
  )
  
  # 3. Edge Type Analysis
  edge_type_stats <- as.data.frame(table(E(long_range_net)$edge_type))
  colnames(edge_type_stats) <- c("type", "count")
  edge_type_stats$percentage <- edge_type_stats$count / sum(edge_type_stats$count) * 100
  
  # 4. Distance Analysis by Edge Type
  distance_by_type <- data.frame(
    distance = E(long_range_net)$distance,
    type = E(long_range_net)$edge_type
  )
  
  # 5. Node Analysis
  node_metrics <- data.frame(
    node = V(long_range_net)$name,
    degree = degree(long_range_net),
    betweenness = betweenness(long_range_net, normalized = TRUE),
    clustering = transitivity(long_range_net, type = "local"),
    node_type = V(long_range_net)$node_type,
    start = V(long_range_net)$start,
    end = V(long_range_net)$end,
    genes = V(long_range_net)$genes,
    gene_types = V(long_range_net)$gene_types
  )
  
  # 6. Coding Region Analysis
  coding_nodes <- node_metrics[node_metrics$node_type == "C",]
  coding_stats <- list(
    num_coding_regions = nrow(coding_nodes),
    mean_coding_degree = mean(coding_nodes$degree),
    median_coding_degree = median(coding_nodes$degree)
  )
  
  # 7. Gene Content Analysis
  gene_list <- unique(unlist(strsplit(
    paste(V(long_range_net)$genes[V(long_range_net)$node_type == "C"], 
          collapse = ","), 
    ",")))
  gene_list <- gene_list[gene_list != ""]
  
  gene_stats <- list(
    total_genes = length(gene_list),
    genes_per_region = length(gene_list) / coding_stats$num_coding_regions
  )
  
  # Return all results
  results <- list(
    network = long_range_net,
    network_stats = network_stats,
    edge_type_stats = edge_type_stats,
    distance_by_type = distance_by_type,
    node_metrics = node_metrics,
    coding_stats = coding_stats,
    gene_stats = gene_stats
  )
  
  return(results)
}

generate_coding_network_report <- function(analysis_results) {
  cat("=== Long-range Coding Interactions Network Analysis Report ===\n")
  cat(sprintf("Chromosome: %s\n\n", analysis_results$network_stats$name))
  
  # 1. Network Overview
  cat("1. Network Overview\n")
  cat("-----------------\n")
  cat(sprintf("Distance threshold: %d bp\n", 
              round(analysis_results$network_stats$threshold_distance)))
  cat(sprintf("Total nodes: %d\n", analysis_results$network_stats$num_nodes))
  cat(sprintf("Total edges: %d\n", analysis_results$network_stats$num_edges))
  cat(sprintf("Network density: %.3f\n", 
              analysis_results$network_stats$network_density))
  cat("\n")
  
  # 2. Edge Type Distribution
  cat("2. Edge Type Distribution\n")
  cat("-----------------------\n")
  print(analysis_results$edge_type_stats)
  cat("\n")
  
  # 3. Distance Characteristics
  cat("3. Distance Characteristics by Edge Type\n")
  cat("------------------------------------\n")
  dist_summary <- analysis_results$distance_by_type %>%
    group_by(type) %>%
    summarise(
      mean_distance = mean(distance),
      median_distance = median(distance),
      min_distance = min(distance),
      max_distance = max(distance)
    )
  print(dist_summary)
  cat("\n")
  
  # 4. Coding Regions Analysis
  cat("4. Coding Regions Analysis\n")
  cat("------------------------\n")
  cat(sprintf("Number of coding regions: %d\n", 
              analysis_results$coding_stats$num_coding_regions))
  cat(sprintf("Mean degree of coding regions: %.2f\n", 
              analysis_results$coding_stats$mean_coding_degree))
  cat(sprintf("Median degree of coding regions: %.2f\n", 
              analysis_results$coding_stats$median_coding_degree))
  cat("\n")
  
  # 5. Gene Content
  cat("5. Gene Content\n")
  cat("--------------\n")
  cat(sprintf("Total unique genes: %d\n", 
              analysis_results$gene_stats$total_genes))
  cat(sprintf("Average genes per coding region: %.2f\n", 
              analysis_results$gene_stats$genes_per_region))
  cat("\n")
  
  # 6. Top Hub Regions
  cat("6. Top Hub Regions\n")
  cat("----------------\n")
  top_hubs <- head(arrange(analysis_results$node_metrics, 
                           desc(degree)), 5)
  print(dplyr::select(top_hubs, node, node_type, degree, genes))
  cat("\n")
  
  # 7. Network Structure
  cat("7. Network Structure\n")
  cat("------------------\n")
  cat(sprintf("Number of components: %d\n", 
              analysis_results$network_stats$components))
  cat(sprintf("Largest component size: %d nodes\n", 
              analysis_results$network_stats$largest_component_size))
}




# Analyze a single chromosome
chr19_analysis <- analyze_long_range_coding_network(normal_graphs$chr19, "chr19")

# Generate the report
generate_coding_network_report(chr19_analysis)

# To analyze all chromosomes:
all_chr_analyses <- lapply(names(normal_graphs), function(chr) {
  analyze_long_range_coding_network(normal_graphs[[chr]], chr)
})
names(all_chr_analyses) <- names(normal_graphs)



# Define the PDF output file
pdf("output_distance.pdf")

# Redirect printed output to a file
sink("output_distance.pdf")
all_chr_analyses <- lapply(names(normal_graphs), function(chr) {
  analyze_long_range_coding_network(normal_graphs[[chr]], chr)
})

sink()  # Reset sink

# Close the PDF
dev.off()








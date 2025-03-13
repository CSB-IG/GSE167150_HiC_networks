library(igraph)
library(muxViz) # For multilayer network analysis
library(ggplot2)
library(tidyverse)

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


# # 1. Prepare the layers for multilayer analysis
# prepare_multilayer_network <- function(hic_graph, mi_data, chromosome) {
#   # Create gene-to-node mapping
#   gene_to_node <- list()
#   for (i in 1:length(V(hic_graph))) {
#     if (!is.na(V(hic_graph)$gene_ids[i]) && 
#         V(hic_graph)$node_type[i] %in% c("C", "R")) {
#       
#       genes <- unlist(strsplit(V(hic_graph)$gene_ids[i], ";"))
#       for (gene in genes) {
#         if (gene != "") {
#           if (is.null(gene_to_node[[gene]])) {
#             gene_to_node[[gene]] <- c(i)
#           } else {
#             gene_to_node[[gene]] <- c(gene_to_node[[gene]], i)
#           }
#         }
#       }
#     }
#   }
#   
#   # Filter MI data for this chromosome
#   mi_data_chr <- mi_data[mi_data$Chromosome == chromosome,]
#   
#   # Create gene-gene adjacency matrices for both layers
#   # Get unique genes that appear in both networks
#   mi_genes <- unique(c(mi_data_chr$gene1, mi_data_chr$gene2))
#   hic_genes <- names(gene_to_node)
#   common_genes <- intersect(mi_genes, hic_genes)
#   
#   # Create empty adjacency matrices
#   n_genes <- length(common_genes)
#   hic_adj <- matrix(0, n_genes, n_genes)
#   mi_adj <- matrix(0, n_genes, n_genes)
#   rownames(hic_adj) <- colnames(hic_adj) <- common_genes
#   rownames(mi_adj) <- colnames(mi_adj) <- common_genes
#   
#   # Fill MI adjacency matrix
#   for (i in 1:nrow(mi_data_chr)) {
#     gene1 <- mi_data_chr$gene1[i]
#     gene2 <- mi_data_chr$gene2[i]
#     
#     if (gene1 %in% common_genes && gene2 %in% common_genes) {
#       mi_adj[gene1, gene2] <- mi_data_chr$MI[i]
#       mi_adj[gene2, gene1] <- mi_data_chr$MI[i] # Symmetric
#     }
#   }
#   
#   # Fill HiC adjacency matrix (taking max zscore if multiple edges)
#   for (i in 1:length(common_genes)) {
#     for (j in i:length(common_genes)) {
#       gene1 <- common_genes[i]
#       gene2 <- common_genes[j]
#       
#       # Skip if same gene
#       if (gene1 == gene2) next
#       
#       # Get nodes containing these genes
#       nodes1 <- gene_to_node[[gene1]]
#       nodes2 <- gene_to_node[[gene2]]
#       
#       max_zscore <- 0
#       
#       # Check all possible node pairs
#       for (n1 in nodes1) {
#         for (n2 in nodes2) {
#           if (are_adjacent(hic_graph, n1, n2)) {
#             edge_id <- get_edge_ids(hic_graph, c(n1, n2))
#             zscore <- E(hic_graph)$zscore[edge_id]
#             if (!is.na(zscore) && zscore > max_zscore) {
#               max_zscore <- zscore
#             }
#           }
#         }
#       }
#       
#       if (max_zscore > 0) {
#         hic_adj[gene1, gene2] <- max_zscore
#         hic_adj[gene2, gene1] <- max_zscore
#       }
#     }
#   }
#   
#   # Return adjacency matrices and common genes
#   return(list(
#     hic_adj = hic_adj,
#     mi_adj = mi_adj,
#     common_genes = common_genes
#   ))
# }
# 
# # 2. Create multilayer networks for normal and cancer
# create_chromosome_multilayer <- function(normal_hic, normal_mi, tnbc_hic, tnbc_mi, chr) {
#   # Prepare networks
#   normal_layers <- prepare_multilayer_network(normal_hic, normal_mi, chr)
#   tnbc_layers <- prepare_multilayer_network(tnbc_hic, tnbc_mi, chr)
#   
#   # Find genes common to all networks
#   common_genes <- intersect(normal_layers$common_genes, tnbc_layers$common_genes)
#   
#   # Filter adjacency matrices to common genes
#   normal_hic_adj <- normal_layers$hic_adj[common_genes, common_genes]
#   normal_mi_adj <- normal_layers$mi_adj[common_genes, common_genes]
#   tnbc_hic_adj <- tnbc_layers$hic_adj[common_genes, common_genes]
#   tnbc_mi_adj <- tnbc_layers$mi_adj[common_genes, common_genes]
#   
#   # Format for muxViz (if using) or return adjacency matrices
#   return(list(
#     normal = list(hic = normal_hic_adj, mi = normal_mi_adj),
#     tnbc = list(hic = tnbc_hic_adj, mi = tnbc_mi_adj),
#     common_genes = common_genes
#   ))
# }
# 
# # 3. Compute layer correlations and overlaps
# analyze_layer_relationships <- function(multilayer_net) {
#   # For each condition (normal, TNBC)
#   results <- list()
#   
#   for (condition in c("normal", "tnbc")) {
#     # Get adjacency matrices
#     hic_adj <- multilayer_net[[condition]]$hic
#     mi_adj <- multilayer_net[[condition]]$mi
#     
#     # Convert to edge lists for analysis
#     hic_edges <- which(hic_adj > 0, arr.ind = TRUE)
#     mi_edges <- which(mi_adj > 0, arr.ind = TRUE)
#     
#     # Calculate edge overlap
#     hic_edge_set <- paste(pmin(hic_edges[,1], hic_edges[,2]), 
#                           pmax(hic_edges[,1], hic_edges[,2]), sep="_")
#     mi_edge_set <- paste(pmin(mi_edges[,1], mi_edges[,2]), 
#                          pmax(mi_edges[,1], mi_edges[,2]), sep="_")
#     
#     # Calculate Jaccard similarity between edge sets
#     overlap_edges <- intersect(hic_edge_set, mi_edge_set)
#     jaccard_sim <- length(overlap_edges) / length(union(hic_edge_set, mi_edge_set))
#     
#     # Calculate correlation between edge weights (for overlapping edges)
#     overlap_indices <- match(overlap_edges, hic_edge_set)
#     hic_weights <- sapply(overlap_indices, function(i) {
#       if(is.na(i)) return(NA)
#       row <- hic_edges[i,1]
#       col <- hic_edges[i,2]
#       return(hic_adj[row, col])
#     })
#     
#     overlap_indices <- match(overlap_edges, mi_edge_set)
#     mi_weights <- sapply(overlap_indices, function(i) {
#       if(is.na(i)) return(NA)
#       row <- mi_edges[i,1]
#       col <- mi_edges[i,2]
#       return(mi_adj[row, col])
#     })
#     
#     weight_correlation <- cor(hic_weights, mi_weights, method="spearman", use="complete.obs")
#     
#     # Store results
#     results[[condition]] <- list(
#       jaccard_similarity = jaccard_sim,
#       weight_correlation = weight_correlation,
#       overlap_size = length(overlap_edges),
#       hic_edge_count = length(hic_edge_set),
#       mi_edge_count = length(mi_edge_set),
#       hic_weights = hic_weights,
#       mi_weights = mi_weights
#     )
#   }
#   
#   return(results)
# }
# 
# # 4. Compare normal vs cancer multilayer networks
# compare_conditions <- function(multilayer_results) {
#   # Compare Jaccard similarities
#   normal_jaccard <- multilayer_results$normal$jaccard_similarity
#   tnbc_jaccard <- multilayer_results$tnbc$jaccard_similarity
#   
#   # Compare weight correlations
#   normal_correlation <- multilayer_results$normal$weight_correlation
#   tnbc_correlation <- multilayer_results$tnbc$weight_correlation
#   
#   # Visualize comparison
#   comparison_df <- data.frame(
#     Metric = c("Jaccard Similarity", "Weight Correlation"),
#     Normal = c(normal_jaccard, normal_correlation),
#     TNBC = c(tnbc_jaccard, tnbc_correlation)
#   )
#   
#   comparison_long <- pivot_longer(comparison_df, cols=c("Normal", "TNBC"), 
#                                   names_to="Condition", values_to="Value")
#   
#   ggplot(comparison_long, aes(x=Metric, y=Value, fill=Condition)) +
#     geom_bar(stat="identity", position="dodge") +
#     theme_minimal() +
#     labs(title="Comparison of Layer Relationships: Normal vs TNBC",
#          y="Value") +
#     scale_fill_brewer(palette="Set1")
# }
# 
# 
# # 5. Identify genes with significantly altered relationships
# find_altered_relationships <- function(multilayer_net) {
#   common_genes <- multilayer_net$common_genes
#   normal_hic <- multilayer_net$normal$hic
#   normal_mi <- multilayer_net$normal$mi
#   tnbc_hic <- multilayer_net$tnbc$hic
#   tnbc_mi <- multilayer_net$tnbc$mi
#   
#   # Calculate gene-wise metrics
#   gene_results <- data.frame(
#     gene = common_genes,
#     normal_hic_degree = colSums(normal_hic > 0),
#     normal_mi_degree = colSums(normal_mi > 0),
#     tnbc_hic_degree = colSums(tnbc_hic > 0),
#     tnbc_mi_degree = colSums(tnbc_mi > 0)
#   )
#   
#   # Calculate ratios of MI to HiC connectivity
#   gene_results$normal_ratio <- gene_results$normal_mi_degree / 
#     (gene_results$normal_hic_degree + 0.1) # Add 0.1 to avoid division by zero
#   gene_results$tnbc_ratio <- gene_results$tnbc_mi_degree / 
#     (gene_results$tnbc_hic_degree + 0.1)
#   
#   # Calculate change in ratio
#   gene_results$ratio_change <- gene_results$tnbc_ratio - gene_results$normal_ratio
#   
#   # Identify genes with biggest changes
#   gene_results <- gene_results[order(-abs(gene_results$ratio_change)),]
#   
#   return(gene_results)
# }
# 
# # 6. Create multilayer network visualization
# visualize_multilayer <- function(multilayer_net, condition, top_genes=30) {
#   # Get adjacency matrices
#   hic_adj <- multilayer_net[[condition]]$hic
#   mi_adj <- multilayer_net[[condition]]$mi
#   common_genes <- multilayer_net$common_genes
#   
#   # Select top genes by total degree
#   total_degree <- colSums(hic_adj > 0) + colSums(mi_adj > 0)
#   top_indices <- order(-total_degree)[1:min(top_genes, length(common_genes))]
#   
#   # Subset matrices
#   hic_sub <- hic_adj[top_indices, top_indices]
#   mi_sub <- mi_adj[top_indices, top_indices]
#   selected_genes <- common_genes[top_indices]
#   
#   # Create igraph objects
#   hic_graph <- graph_from_adjacency_matrix(hic_sub, weighted=TRUE, mode="undirected")
#   mi_graph <- graph_from_adjacency_matrix(mi_sub, weighted=TRUE, mode="undirected")
#   
#   # Set vertex names
#   V(hic_graph)$name <- V(mi_graph)$name <- selected_genes
#   
#   # Scale edge weights for visualization
#   E(hic_graph)$width <- scales::rescale(E(hic_graph)$weight, to=c(0.5, 3))
#   E(mi_graph)$width <- scales::rescale(E(mi_graph)$weight, to=c(0.5, 3))
#   
#   # Set colors for different layers
#   E(hic_graph)$color <- "darkorchid4"
#   E(mi_graph)$color <- "darkgreen"
#   
#   # Use same layout for both
#   layout <- layout_with_fr(hic_graph)
#   
#   # Create combined graph for visualization
#   combined <- hic_graph
#   for (i in 1:ecount(mi_graph)) {
#     e <- ends(mi_graph, i)
#     if (!are_adjacent(combined, e[1], e[2])) {
#       combined <- add_edges(combined, e, color="darkgreen", width=E(mi_graph)$width[i])
#     }
#   }
#   
#   # Plot
#   plot(combined, layout=layout, vertex.label.cex=0.7,
#        main=paste("Multilayer Network -", condition, "condition"),
#        edge.curved=0.2)
#   legend("bottomright", c("HiC", "MI"), lty=1, col=c("darkorchid4", "darkgreen"))
# }
# 
# # 7. Execute full analysis
# run_multilayer_analysis <- function(normal_graphs, mi_normal, tnbc_graphs, mi_tnbc) {
#   # Get chromosomes
#   chromosomes <- names(normal_graphs)
#   
#   all_results <- list()
#   
#   for (chr in chromosomes) {
#     cat("Processing", chr, "...\n")
#     
#     # Create multilayer network
#     multilayer <- create_chromosome_multilayer(
#       normal_graphs[[chr]], 
#       mi_normal[[chr]], 
#       tnbc_graphs[[chr]], 
#       mi_tnbc[[chr]], 
#       chr
#     )
#     
#     # Analyze layer relationships
#     layer_results <- analyze_layer_relationships(multilayer)
#     
#     # Find genes with altered relationships
#     altered_genes <- find_altered_relationships(multilayer)
#     
#     # Store results
#     all_results[[chr]] <- list(
#       multilayer = multilayer,
#       layer_relationships = layer_results,
#       altered_genes = altered_genes
#     )
#     
#     # Visualize (for selected chromosomes)
#     if (chr %in% c("chr1", "chr17")) {
#       pdf(paste0("multilayer_", chr, "_normal.pdf"))
#       visualize_multilayer(multilayer, "normal")
#       dev.off()
#       
#       pdf(paste0("multilayer_", chr, "_tnbc.pdf"))
#       visualize_multilayer(multilayer, "tnbc")
#       dev.off()
#       
#       # Compare conditions
#       pdf(paste0("comparison_", chr, ".pdf"))
#       compare_conditions(layer_results)
#       dev.off()
#     }
#   }
#   
#   return(all_results)
# }
# 
# # Run the analysis
# results <- run_multilayer_analysis(normal_graphs, mi_normal, tnbc_graphs, mi_tnbc)
# 
# saveRDS(results, file="results/multilayer/multilayer_results.Rds")
# 
# # run_multilayer_analysis(normal_graphs[21:22], mi_normal[21:22],
# #                         tnbc_graphs[21:22], mi_tnbc[21:22])

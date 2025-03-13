###########################################################################
# mi_and_hic.R
#
# analysis of mi and hic networks
# Resolution is 40kb
#
###########################################################################
#setwd("/efs/users/hreyes/tcga_brca_SE/")

# libraries
library(igraph)
library(data.table)

#
#
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

# how many genes do the Hi-C Networks have per node
coding <- tcga_annotation$gene_id[tcga_annotation$gene_type == "protein_coding"] 

# lapply(
#   V(normal_graphs$chr1)$gene_ids[V(normal_graphs$chr1)$node_type != "N"], 
#   function(i) {
#     strsplit(i, ";") %>%
#       unlist() -> gnodes
#   }
# ) -> n_genes_data
# 
# unlist(lapply(n_genes_data, function(n) {
#   genes <- n %in% coding
#   
#   sum(genes)
# })) -> n_genes_data
# 
# ng <- table(n_genes_data[n_genes_data != 0])
# 
# names(ng) <- paste0(names(ng), "\n(", ng, " nodes)")
# 
# barplot(ng, las=1, xlab = "Number of Coding Genes per Node", ylab="Nodes")


lapply(chromosomes, function(c) {
  lapply(V(normal_graphs[[c]])$gene_ids[V(normal_graphs[[c]])$node_type != "N"],
         function(i) {
           strsplit(i, ";") %>%
             unlist()
         }) %>%
    lapply(., function(n) { 
      n %in% coding %>%
        sum()
      }) %>%
    unlist() -> n_genes
  
  table(n_genes[n_genes != 0]) %>% 
    as.data.frame() -> df
  
  colnames(df) <- c("Genes", "Nodes")
  
  data.frame(df, Chromosome=c)
}) -> genes_per_node_df

genes_per_node_df <- do.call(rbind, genes_per_node_df)

genes_per_node_df$Chromosome <- factor(genes_per_node_df$Chromosome, levels = chromosomes)

svg(file.path("results/nodes_hic_genes", "barplot_genes_per_node.svg"), width = 10, height = 10)

ggplot(genes_per_node_df, aes(Genes, Nodes)) + geom_bar(stat = "identity") +
  facet_wrap(~Chromosome, scales = "free") +
  xlab("Genes per node") + ylab("Count") +
  ggtitle("Number of Coding Genes within each Hi-C 40kb Node")

dev.off()

#######################################
# add MI as attribute to hi-c networks

hic_graph= tnbc_graphs$chr1
mi_network = mi_tnbc$chr1


assign_mi_to_hic_network <- function(hic_graph, mi_network) {
  # Filter edges of type "C-C"
  cc_edges <- which(E(hic_graph)$edge_type == "C-C")
  
  # Process C-C edges
  mi_values <- lapply(cc_edges, function(edge_idx) {
    # Get the two nodes for this edge
    edge_nodes <- ends(hic_graph, edge_idx)
    
    # Get gene IDs for each node
    genes1 <- strsplit(V(hic_graph)$gene_ids[match(edge_nodes[1], V(hic_graph)$name)], ";")[[1]]
    genes2 <- strsplit(V(hic_graph)$gene_ids[match(edge_nodes[2], V(hic_graph)$name)], ";")[[1]]
    
    # If both nodes have only one gene
    if (length(genes1) == 1 && length(genes2) == 1) {
      # Look up MI value
      mi_row <- mi_network[
        (mi_network$gene1 == genes1 & mi_network$gene2 == genes2) |
          (mi_network$gene1 == genes2 & mi_network$gene2 == genes1),
      ]
      
      # If MI value exists, return it, otherwise return NA
      if (nrow(mi_row) > 0) {
        return(list(
          mi = mi_row$MI[1],
          pvalue = mi_row$pvalue[1]
        ))
      } else {
        return(list(
          mi = NA_real_,
          pvalue = NA_real_
        ))
      }
    } 
    # If multiple genes in either node
    else {
      # Generate all gene pairs
      gene_pairs <- expand.grid(gene1 = genes1, gene2 = genes2, stringsAsFactors = FALSE)
      
      # Look up MI values for all pairs
      gene_pairs$mi_value <- apply(gene_pairs, 1, function(pair) {
        mi_row <- mi_network[
          (mi_network$gene1 == pair['gene1'] & mi_network$gene2 == pair['gene2']) |
            (mi_network$gene1 == pair['gene2'] & mi_network$gene2 == pair['gene1']),
        ]
        
        return(ifelse(nrow(mi_row) > 0, mi_row$MI[1], NA_real_))
      })
      
      gene_pairs$pvalue <- apply(gene_pairs, 1, function(pair) {
        mi_row <- mi_network[
          (mi_network$gene1 == pair['gene1'] & mi_network$gene2 == pair['gene2']) |
            (mi_network$gene1 == pair['gene2'] & mi_network$gene2 == pair['gene1']),
        ]
        
        return(ifelse(nrow(mi_row) > 0, mi_row$pvalue[1], NA_real_))
      })
      
      # Remove rows with NA MI values
      gene_pairs <- gene_pairs[!is.na(gene_pairs$mi_value), ]
      
      # If any valid pairs exist, sort and select the first
      if (nrow(gene_pairs) > 0) {
        sorted_pairs <- gene_pairs[order(gene_pairs$pvalue, -gene_pairs$mi_value), ]
        return(list(
          mi = sorted_pairs$mi_value[1],
          pvalue = sorted_pairs$pvalue[1]
        ))
      } else {
        return(list(
          mi = NA_real_,
          pvalue = NA_real_
        ))
      }
    }
  })
  
  # Initialize edge attributes
  E(hic_graph)$MI <- rep(NA_real_, ecount(hic_graph))
  E(hic_graph)$MI_pvalue <- rep(NA_real_, ecount(hic_graph))
  
  # Assign values to C-C edges
  E(hic_graph)$MI[cc_edges] <- sapply(mi_values, `[[`, "mi")
  E(hic_graph)$MI_pvalue[cc_edges] <- sapply(mi_values, `[[`, "pvalue")
  
  return(hic_graph)
}

# ### attempt 
# test <- assign_mi_to_hic_network(
#   tnbc_graphs$chr21, 
#   mi_tnbc$chr21
# )
# 
# NORMAL
lapply(chromosomes, function(c) {
  o <- assign_mi_to_hic_network(
    normal_graphs[[c]],
    mi_normal[[c]]
  )
}) -> normal_hic_mi

# TNBC
lapply(chromosomes, function(c) {
  o <- assign_mi_to_hic_network(
    tnbc_graphs[[c]],
    mi_tnbc[[c]]
  )
}) -> tnbc_hic_mi

##################################################
lapply(chromosomes, function(c) {
  
  hic_graph= normal_graphs[[c]]
  mi_network = mi_normal[[c]]
  
  # Filter edges based on edge_type attribute
  selected_edges <- which(E(hic_graph)$edge_type == "C-C")
  
  # Get the nodes connected by these edges
  node1 <- ends(hic_graph, selected_edges)[,1]
  node2 <- ends(hic_graph, selected_edges)[,2]
  
  # Extract attributes for the nodes
  node1_gene_ids <- V(hic_graph)[node1]$gene_ids
  node2_gene_ids <- V(hic_graph)[node2]$gene_ids
  chromosome <- V(hic_graph)[node1]$chr  # Assuming both nodes share the same chromosome
  zscore <- E(hic_graph)[selected_edges]$zscore
  
  # Create a data frame
  result_df <- data.frame(
    node1 = node1,
    node2 = node2,
    node1_gene_ids = node1_gene_ids,
    node2_gene_ids = node2_gene_ids,
    zscore = zscore,
    chromosome = chromosome
  )
  
  # View the result
  #head(result_df)
  
  multi1 <- str_count(result_df$node1_gene_ids, ";")
  multi2 <- str_count(result_df$node2_gene_ids, ";")
  
  #table(multi1, multi2)
  
  result_df <- result_df[which(multi1 == 0 & multi2 == 0), ]
  
  merge_interaction_data <- function(result_df, mi_network) {
    # Ensure gene pairs are comparable regardless of order
    mi_network$key <- paste0(pmin(mi_network$gene1, mi_network$gene2), "_", 
                             pmax(mi_network$gene1, mi_network$gene2))
    
    result_df$key <- paste0(pmin(result_df$node1_gene_ids, result_df$node2_gene_ids), "_", 
                            pmax(result_df$node1_gene_ids, result_df$node2_gene_ids))
    
    # Merge based on the key
    merged_df <- merge(result_df, mi_network[, c("key", "MI", "pvalue")], by = "key", all.x = TRUE)
    
    # Remove the helper key column
    #merged_df$key <- NULL
    
    return(merged_df)
    
  }
  
  # Apply the function
  result_df <- merge_interaction_data(result_df, mi_network)
  
  return(result_df)
}) -> mi_hic_normal

mi_hic_normal <- do.call(rbind, mi_hic_normal)
mi_hic_normal <- data.frame(mi_hic_normal, Phenotype="Normal")


lapply(chromosomes, function(c) {

  hic_graph= tnbc_graphs[[c]]
  mi_network = mi_tnbc[[c]]
  
  # Filter edges based on edge_type attribute
  selected_edges <- which(E(hic_graph)$edge_type == "C-C")
  
  # Get the nodes connected by these edges
  node1 <- ends(hic_graph, selected_edges)[,1]
  node2 <- ends(hic_graph, selected_edges)[,2]
  
  # Extract attributes for the nodes
  node1_gene_ids <- V(hic_graph)[node1]$gene_ids
  node2_gene_ids <- V(hic_graph)[node2]$gene_ids
  chromosome <- V(hic_graph)[node1]$chr  # Assuming both nodes share the same chromosome
  zscore <- E(hic_graph)[selected_edges]$zscore
  
  # Create a data frame
  result_df <- data.frame(
    node1 = node1,
    node2 = node2,
    node1_gene_ids = node1_gene_ids,
    node2_gene_ids = node2_gene_ids,
    zscore = zscore,
    chromosome = chromosome
    )
  
  # View the result
  #head(result_df)
  
  multi1 <- str_count(result_df$node1_gene_ids, ";")
  multi2 <- str_count(result_df$node2_gene_ids, ";")
  
  #table(multi1, multi2)
  
  result_df <- result_df[which(multi1 == 0 & multi2 == 0), ]
  
  merge_interaction_data <- function(result_df, mi_network) {
    # Ensure gene pairs are comparable regardless of order
    mi_network$key <- paste0(pmin(mi_network$gene1, mi_network$gene2), "_", 
                             pmax(mi_network$gene1, mi_network$gene2))
    
    result_df$key <- paste0(pmin(result_df$node1_gene_ids, result_df$node2_gene_ids), "_", 
                            pmax(result_df$node1_gene_ids, result_df$node2_gene_ids))
    
    # Merge based on the key
    merged_df <- merge(result_df, mi_network[, c("key", "MI", "pvalue")], by = "key", all.x = TRUE)
    
    # Remove the helper key column
    #merged_df$key <- NULL
    
    return(merged_df)
    
  }
  
  # Apply the function
  result_df <- merge_interaction_data(result_df, mi_network)
  
  return(result_df)
}) -> mi_hic_tnbc

mi_hic_tnbc <- do.call(rbind, mi_hic_tnbc)
mi_hic_tnbc <- data.frame(mi_hic_tnbc, Phenotype="TNBC")


mi_hi_all <- rbind(mi_hic_normal, mi_hic_tnbc)
mi_hi_all$chromosome <- factor(mi_hic_all$chromosome, levels=chromosomes)



##### scatterplots

# # plot 
# ggplot(mi_hi_all, aes(zscore, MI, col=Phenotype)) + geom_point() +
#   facet_wrap(~chromosome, scales="free") + theme_minimal()
# 

svg(file.path("results/nodes_hic_genes", "scatterplot_zscore_mi.svg"), width = 10, height = 10)

ggplot(mi_hi_all, aes(x = zscore, y = MI, color = Phenotype)) +
  geom_point(alpha = 0.4, size = 2) +  # Scatter points
#  geom_smooth(method = "lm", se = FALSE) +  # Linear trend lines per phenotype
  theme_minimal() + facet_wrap(~chromosome, scales = "free") +
  labs(title = "Relationship between Hi-C count z-score and Mutual Information",
       x = "Edge Z-score",
       y = "Edge Mutual Information (MI)") +
  theme(legend.position = "top")  +
  scale_color_manual(values = c("Normal" = "#008080", "TNBC" = "darkred"))  # Adjust colors as needed
dev.off()


#################################################################################
#################################################################################
#################################################################################
#' Comprehensive Multilayer Network Analysis for Hi-C and MI Networks
#' 
#' This script implements a multilayer network analysis framework to study
#' the relationship between chromatin conformation (Hi-C) and gene expression
#' correlation (MI) networks in normal and cancer conditions.

library(igraph)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(RColorBrewer)

#' Create a multilayer network from Hi-C and MI data
#' 
#' @param hic_graph An igraph object containing the Hi-C network
#' @param mi_data Data frame with MI values between gene pairs
#' @param chromosome Chromosome identifier (e.g., "chr1")
#' @return A list containing the multilayer network and analysis results
create_multilayer_network <- function(hic_graph, mi_data, chromosome) {
  # 1. Map genes to Hi-C nodes
  gene_to_node <- list()
  for (i in 1:vcount(hic_graph)) {
    if (!is.na(V(hic_graph)$gene_ids[i]) && V(hic_graph)$node_type[i] %in% c("C", "R")) {
      genes <- unlist(strsplit(V(hic_graph)$gene_ids[i], ";"))
      for (gene in genes) {
        if (gene != "") {
          if (is.null(gene_to_node[[gene]])) {
            gene_to_node[[gene]] <- c(i)
          } else {
            gene_to_node[[gene]] <- c(gene_to_node[[gene]], i)
          }
        }
      }
    }
  }
  
  # 2. Filter MI data for this chromosome
  mi_chr <- mi_data[mi_data$Chromosome == chromosome, ]
  
  # 3. Create gene set intersection
  mi_genes <- unique(c(mi_chr$gene1, mi_chr$gene2))
  hic_genes <- names(gene_to_node)
  common_genes <- intersect(mi_genes, hic_genes)
  
  # 4. Create gene-gene adjacency matrices
  n_genes <- length(common_genes)
  gene_indices <- setNames(1:n_genes, common_genes)
  
  # Initialize sparse matrices
  hic_adj <- Matrix(0, n_genes, n_genes, sparse = TRUE)
  mi_adj <- Matrix(0, n_genes, n_genes, sparse = TRUE)
  
  # 5. Fill MI adjacency matrix with MI values
  for (i in 1:nrow(mi_chr)) {
    gene1 <- mi_chr$gene1[i]
    gene2 <- mi_chr$gene2[i]
    
    if (gene1 %in% common_genes && gene2 %in% common_genes) {
      idx1 <- gene_indices[gene1]
      idx2 <- gene_indices[gene2]
      mi_adj[idx1, idx2] <- mi_chr$MI[i]
      mi_adj[idx2, idx1] <- mi_chr$MI[i]  # Symmetric
    }
  }
  
  # 6. Fill Hi-C adjacency matrix
  for (i in 1:(n_genes-1)) {
    for (j in (i+1):n_genes) {
      gene1 <- common_genes[i]
      gene2 <- common_genes[j]
      
      # Get nodes containing these genes
      nodes1 <- gene_to_node[[gene1]]
      nodes2 <- gene_to_node[[gene2]]
      
      max_zscore <- 0
      
      # Check all possible node pairs
      for (n1 in nodes1) {
        for (n2 in nodes2) {
          if (are.connected(hic_graph, n1, n2)) {
            edge_id <- get.edge.ids(hic_graph, c(n1, n2))
            zscore <- E(hic_graph)$zscore[edge_id]
            if (!is.na(zscore) && zscore > max_zscore) {
              max_zscore <- zscore
            }
          }
        }
      }
      
      if (max_zscore > 0) {
        hic_adj[i, j] <- max_zscore
        hic_adj[j, i] <- max_zscore
      }
    }
  }
  
  # 7. Create multilayer network representation
  multilayer <- list(
    hic_layer = hic_adj,
    mi_layer = mi_adj,
    genes = common_genes,
    gene_indices = gene_indices
  )
  
  return(multilayer)
}

#' Calculate multilayer network metrics
#' 
#' @param multilayer A multilayer network object
#' @return A data frame of gene-level multilayer metrics
calculate_multilayer_metrics <- function(multilayer) {
  hic_adj <- multilayer$hic_layer
  mi_adj <- multilayer$mi_layer
  genes <- multilayer$genes
  
  # 1. Degree-based metrics
  hic_degree <- rowSums(hic_adj > 0)
  mi_degree <- rowSums(mi_adj > 0)
  
  # 2. Strength-based metrics (weighted sum)
  hic_strength <- rowSums(hic_adj)
  mi_strength <- rowSums(mi_adj)
  
  # 3. Cross-layer metrics
  # Normalized adjacency matrices
  hic_norm <- hic_adj / max(hic_adj)
  mi_norm <- mi_adj / max(mi_adj)
  
  # Multilayer participation coefficient
  participation <- 1 - ((hic_degree/(hic_degree + mi_degree))^2 + 
                          (mi_degree/(hic_degree + mi_degree))^2)
  participation[is.nan(participation)] <- 0
  
  # 4. Layer agreement
  # Convert to igraph objects
  hic_graph <- graph_from_adjacency_matrix(hic_adj, weighted = TRUE, mode = "undirected")
  mi_graph <- graph_from_adjacency_matrix(mi_adj, weighted = TRUE, mode = "undirected")
  
  # Calculate eigenvector centrality
  hic_eigen <- tryCatch({
    eigen_centrality(hic_graph)$vector
  }, error = function(e) {
    rep(0, vcount(hic_graph))
  })
  
  mi_eigen <- tryCatch({
    eigen_centrality(mi_graph)$vector
  }, error = function(e) {
    rep(0, vcount(mi_graph))
  })
  
  # Layer agreement
  layer_agreement <- abs(scale(hic_eigen) - scale(mi_eigen))
  
  # 5. Cross-layer assortativity
  # Measure if genes with high HiC connectivity also have high MI connectivity
  cross_assortativity <- cor(hic_degree, mi_degree, method = "spearman")
  
  # Compile metrics
  metrics <- data.frame(
    gene = genes,
    hic_degree = hic_degree,
    mi_degree = mi_degree,
    hic_strength = hic_strength,
    mi_strength = mi_strength,
    total_degree = hic_degree + mi_degree,
    participation_coef = participation,
    layer_agreement = layer_agreement,
    hic_eigen = hic_eigen,
    mi_eigen = mi_eigen
  )
  
  # Add global metrics as attributes
  attr(metrics, "cross_assortativity") <- cross_assortativity
  
  return(metrics)
}

#' Compare multilayer networks between conditions
#' 
#' @param normal_multilayer Multilayer network for normal condition
#' @param cancer_multilayer Multilayer network for cancer condition
#' @return Comparison results and differential metrics
compare_multilayer_networks <- function(normal_multilayer, cancer_multilayer) {
  # 1. Ensure common gene set
  normal_genes <- normal_multilayer$genes
  cancer_genes <- cancer_multilayer$genes
  common_genes <- intersect(normal_genes, cancer_genes)
  
  # Get indices for common genes
  normal_indices <- match(common_genes, normal_genes)
  cancer_indices <- match(common_genes, cancer_genes)
  
  # 2. Extract relevant submatrices
  normal_hic <- normal_multilayer$hic_layer[normal_indices, normal_indices]
  normal_mi <- normal_multilayer$mi_layer[normal_indices, normal_indices]
  cancer_hic <- cancer_multilayer$hic_layer[cancer_indices, cancer_indices]
  cancer_mi <- cancer_multilayer$mi_layer[cancer_indices, cancer_indices]
  
  # 3. Calculate differential matrices
  diff_hic <- cancer_hic - normal_hic
  diff_mi <- cancer_mi - normal_mi
  
  # 4. Calculate metrics for both conditions
  normal_metrics <- calculate_multilayer_metrics(list(
    hic_layer = normal_hic,
    mi_layer = normal_mi,
    genes = common_genes,
    gene_indices = setNames(1:length(common_genes), common_genes)
  ))
  
  cancer_metrics <- calculate_multilayer_metrics(list(
    hic_layer = cancer_hic,
    mi_layer = cancer_mi,
    genes = common_genes,
    gene_indices = setNames(1:length(common_genes), common_genes)
  ))
  
  # 5. Calculate differential metrics
  diff_metrics <- data.frame(
    gene = common_genes,
    diff_hic_degree = cancer_metrics$hic_degree - normal_metrics$hic_degree,
    diff_mi_degree = cancer_metrics$mi_degree - normal_metrics$mi_degree,
    diff_participation = cancer_metrics$participation_coef - normal_metrics$participation_coef,
    diff_layer_agreement = cancer_metrics$layer_agreement - normal_metrics$layer_agreement,
    fold_change_hic = log2((cancer_metrics$hic_strength + 0.1) / (normal_metrics$hic_strength + 0.1)),
    fold_change_mi = log2((cancer_metrics$mi_strength + 0.1) / (normal_metrics$mi_strength + 0.1))
  )
  
  # 6. Calculate global network differences
  global_diff <- list(
    cross_assortativity_normal = attr(normal_metrics, "cross_assortativity"),
    cross_assortativity_cancer = attr(cancer_metrics, "cross_assortativity"),
    hic_density_normal = sum(normal_hic > 0) / (length(common_genes) * (length(common_genes) - 1)),
    hic_density_cancer = sum(cancer_hic > 0) / (length(common_genes) * (length(common_genes) - 1)),
    mi_density_normal = sum(normal_mi > 0) / (length(common_genes) * (length(common_genes) - 1)),
    mi_density_cancer = sum(cancer_mi > 0) / (length(common_genes) * (length(common_genes) - 1))
  )
  
  # 7. Identify significantly altered genes
  # Z-score of differential metrics
  diff_metrics$z_combined <- scale(abs(diff_metrics$diff_hic_degree) + 
                                     abs(diff_metrics$diff_mi_degree))
  
  # Top altered genes
  top_altered_genes <- diff_metrics[order(-diff_metrics$z_combined), ]
  
  return(list(
    normal_metrics = normal_metrics,
    cancer_metrics = cancer_metrics,
    diff_metrics = diff_metrics,
    global_diff = global_diff,
    top_altered_genes = top_altered_genes
  ))
}

#' Visualize multilayer networks
#' 
#' @param multilayer A multilayer network object
#' @param condition Label for the condition (e.g., "Normal" or "Cancer")
#' @param top_n Number of top genes to include in visualization
#' @return A ggplot object
visualize_multilayer_network <- function(multilayer, condition, top_n = 50) {
  # 1. Get top nodes by total degree
  metrics <- calculate_multilayer_metrics(multilayer)
  top_genes <- metrics[order(-metrics$total_degree), ][1:min(top_n, nrow(metrics)), "gene"]
  
  # 2. Subset matrices to top genes
  gene_indices <- match(top_genes, multilayer$genes)
  hic_sub <- multilayer$hic_layer[gene_indices, gene_indices]
  mi_sub <- multilayer$mi_layer[gene_indices, gene_indices]
  
  # 3. Create edge lists for visualization
  hic_edges <- which(hic_sub > 0, arr.ind = TRUE)
  mi_edges <- which(mi_sub > 0, arr.ind = TRUE)
  
  hic_df <- data.frame(
    from = top_genes[hic_edges[, 1]],
    to = top_genes[hic_edges[, 2]],
    weight = hic_sub[hic_edges],
    type = "HiC"
  )
  
  mi_df <- data.frame(
    from = top_genes[mi_edges[, 1]],
    to = top_genes[mi_edges[, 2]],
    weight = mi_sub[mi_edges],
    type = "MI"
  )
  
  # 4. Combine edges
  all_edges <- rbind(hic_df, mi_df)
  
  # 5. Create graph
  g <- graph_from_data_frame(all_edges, directed = FALSE, 
                             vertices = data.frame(name = top_genes))
  
  # 6. Compute layout
  set.seed(42)  # For reproducibility
  layout <- layout_with_fr(g)
  
  # 7. Create plot data
  # For edges
  edge_data <- as_data_frame(g)
  edge_data$x1 <- layout[edge_data$from, 1]
  edge_data$y1 <- layout[edge_data$from, 2]
  edge_data$x2 <- layout[edge_data$to, 1]
  edge_data$y2 <- layout[edge_data$to, 2]
  
  # For vertices
  vertex_data <- data.frame(
    name = V(g)$name,
    x = layout[, 1],
    y = layout[, 2],
    hic_degree = metrics$hic_degree[match(V(g)$name, metrics$gene)],
    mi_degree = metrics$mi_degree[match(V(g)$name, metrics$gene)]
  )
  
  # 8. Create plot
  p <- ggplot() +
    # Add edges
    geom_segment(data = edge_data, aes(x = x1, y = y1, xend = x2, yend = y2, 
                                       color = type, alpha = weight),
                 size = 0.5) +
    # Add vertices
    geom_point(data = vertex_data, aes(x = x, y = y, size = hic_degree + mi_degree),
               color = "darkblue", alpha = 0.7) +
    # Add labels for important nodes
    geom_text(data = vertex_data[order(-vertex_data$hic_degree - vertex_data$mi_degree), ][1:10, ],
              aes(x = x, y = y, label = name), size = 3, hjust = -0.3) +
    # Styling
    scale_color_manual(values = c("HiC" = "blue", "MI" = "red")) +
    scale_alpha(range = c(0.1, 0.8)) +
    scale_size(range = c(1, 6)) +
    theme_void() +
    labs(title = paste("Multilayer Network -", condition),
         subtitle = paste0(top_n, " genes with highest total degree")) +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Identify multilayer communities
#'
#' @param multilayer A multilayer network object
#' @return Community detection results
detect_multilayer_communities <- function(multilayer) {
  # 1. Create aggregate network (sum of both layers with weights)
  hic_norm <- multilayer$hic_layer / max(multilayer$hic_layer)
  mi_norm <- multilayer$mi_layer / max(multilayer$mi_layer)
  
  agg_adj <- hic_norm + mi_norm
  
  # 2. Create graph
  agg_graph <- graph_from_adjacency_matrix(agg_adj, weighted = TRUE, mode = "undirected")
  V(agg_graph)$name <- multilayer$genes
  
  # 3. Run community detection algorithms
  # Louvain method
  louvain <- cluster_louvain(agg_graph)
  
  # Infomap
  infomap <- cluster_infomap(agg_graph)
  
  # 4. Evaluate community quality
  louvain_modularity <- modularity(louvain)
  infomap_modularity <- modularity(infomap)
  
  # 5. Choose best method
  if (louvain_modularity >= infomap_modularity) {
    communities <- louvain
    method <- "louvain"
  } else {
    communities <- infomap
    method <- "infomap"
  }
  
  # 6. Analyze community structure
  community_sizes <- sizes(communities)
  
  # 7. Calculate layer-specific modularity
  hic_graph <- graph_from_adjacency_matrix(multilayer$hic_layer, weighted = TRUE, 
                                           mode = "undirected")
  mi_graph <- graph_from_adjacency_matrix(multilayer$mi_layer, weighted = TRUE, 
                                          mode = "undirected")
  
  hic_modularity <- modularity(hic_graph, membership(communities))
  mi_modularity <- modularity(mi_graph, membership(communities))
  
  # 8. Create community data frame
  community_df <- data.frame(
    gene = multilayer$genes,
    community = membership(communities)
  )
  
  # Add statistics
  attr(community_df, "method") <- method
  attr(community_df, "modularity") <- modularity(communities)
  attr(community_df, "hic_modularity") <- hic_modularity
  attr(community_df, "mi_modularity") <- mi_modularity
  attr(community_df, "community_sizes") <- community_sizes
  
  return(community_df)
}

#' Run comprehensive multilayer analysis pipeline
#'
#' @param normal_graphs List of normal Hi-C networks by chromosome
#' @param mi_normal Normal MI data
#' @param tnbc_graphs List of cancer Hi-C networks by chromosome 
#' @param mi_tnbc Cancer MI data
#' @return Results of multilayer analysis
run_multilayer_analysis <- function(normal_graphs, mi_normal, tnbc_graphs, mi_tnbc) {
  # Get chromosomes
  chromosomes <- names(normal_graphs)
  
  results <- list()
  
  for (chr in chromosomes) {
    cat("Processing", chr, "...\n")
    
    # 1. Create multilayer networks
    normal_multilayer <- create_multilayer_network(
      normal_graphs[[chr]], mi_normal, chr
    )
    
    cancer_multilayer <- create_multilayer_network(
      tnbc_graphs[[chr]], mi_tnbc, chr
    )
    
    # 2. Compare networks
    comparison <- compare_multilayer_networks(
      normal_multilayer, cancer_multilayer
    )
    
    # 3. Detect communities
    normal_communities <- detect_multilayer_communities(normal_multilayer)
    cancer_communities <- detect_multilayer_communities(cancer_multilayer)
    
    # 4. Store results
    results[[chr]] <- list(
      normal_multilayer = normal_multilayer,
      cancer_multilayer = cancer_multilayer,
      comparison = comparison,
      normal_communities = normal_communities,
      cancer_communities = cancer_communities
    )
    
    # 5. Generate key visualizations
    vis_normal <- visualize_multilayer_network(normal_multilayer, "Normal")
    vis_cancer <- visualize_multilayer_network(cancer_multilayer, "Cancer")
    
    # Save plots
    ggsave(paste0("multilayer_", chr, "_normal.pdf"), vis_normal, width = 10, height = 8)
    ggsave(paste0("multilayer_", chr, "_cancer.pdf"), vis_cancer, width = 10, height = 8)
  }
  
  return(results)
}

#' Generate summary report for multilayer analysis
#'
#' @param results Results from run_multilayer_analysis
#' @return Summary data frame
generate_summary_report <- function(results) {
  # Initialize summary data frame
  summary_df <- data.frame(
    chromosome = character(),
    normal_genes = integer(),
    cancer_genes = integer(),
    common_genes = integer(),
    normal_hic_density = numeric(),
    normal_mi_density = numeric(),
    cancer_hic_density = numeric(),
    cancer_mi_density = numeric(),
    normal_cross_assortativity = numeric(),
    cancer_cross_assortativity = numeric(),
    normal_communities = integer(),
    cancer_communities = integer()
  )
  
  # Top altered genes across all chromosomes
  all_top_genes <- data.frame()
  
  for (chr in names(results)) {
    res <- results[[chr]]
    
    # Get key metrics
    normal_genes <- length(res$normal_multilayer$genes)
    cancer_genes <- length(res$cancer_multilayer$genes)
    common_genes <- length(intersect(res$normal_multilayer$genes, res$cancer_multilayer$genes))
    
    # Add to summary
    summary_df <- rbind(summary_df, data.frame(
      chromosome = chr,
      normal_genes = normal_genes,
      cancer_genes = cancer_genes,
      common_genes = common_genes,
      normal_hic_density = res$comparison$global_diff$hic_density_normal,
      normal_mi_density = res$comparison$global_diff$mi_density_normal,
      cancer_hic_density = res$comparison$global_diff$hic_density_cancer,
      cancer_mi_density = res$comparison$global_diff$mi_density_cancer,
      normal_cross_assortativity = res$comparison$global_diff$cross_assortativity_normal,
      cancer_cross_assortativity = res$comparison$global_diff$cross_assortativity_cancer,
      normal_communities = length(attr(res$normal_communities, "community_sizes")),
      cancer_communities = length(attr(res$cancer_communities, "community_sizes"))
    ))
    
    # Add top genes
    top_genes <- res$comparison$top_altered_genes[1:min(10, nrow(res$comparison$top_altered_genes)), ]
    top_genes$chromosome <- chr
    all_top_genes <- rbind(all_top_genes, top_genes)
  }
  
  # Rank genes across all chromosomes
  all_top_genes <- all_top_genes[order(-all_top_genes$z_combined), ]
  
  return(list(
    summary = summary_df,
    top_altered_genes = all_top_genes
  ))
}

# Example usage:
ml_results <- run_multilayer_analysis(normal_graphs, mi_normal, tnbc_graphs, mi_tnbc)
# summary <- generate_summary_report(results)

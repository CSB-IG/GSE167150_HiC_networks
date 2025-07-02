###########################################################################
#
# Summary Statistics of TNBC and Normal Networks
# Resolution is 40kb
#
###########################################################################
# setwd("/efs/users/hreyes/GSE167150_hiedge.git")
#
# 1. Read in annotated hi-c networks of significant contacts
# 2. compute topologicla properties
# 3. visualize
#
# libraries
library(igraph)
library(data.table)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(scales)
library(tidyr)
library(dplyr)
library(reshape2)
library(GenomicRanges)
library(karyoploteR)
library(ggplotify)
library(ggrepel)

######### set up variables
normal_path ="results/igraph_intra/normal/40000/normal_intra_graphs_q001_genes.Rds"
tnbc_path="results/igraph_intra/tnbc/40000/tnbc_intra_graphs_q001_genes.Rds"

outdir_plots = "results/igraph_intra_plots"

my_resolution = 40000

#my_significance = 0.001

# Ensure chromosomes are in the correct order
chromosomes <- paste0("chr", c(seq(1:22), "X"))  

######### Load the data
normal_graphs <- readRDS(normal_path)
tnbc_graphs <- readRDS(tnbc_path)

# annotation
tcga_annotation <- readRDS("results/igraph_intra/tcga_annotation.Rds")

######### FUNCTIONS 
# function to count overlap of genes between annotation and networks
analyze_chromosome_genes <- function(tcga_annotation, normal_networks, tnbc_networks, chromosomes) {
  # Initialize a list to store results that we'll convert to a dataframe
  results <- list()
  
  # For each chromosome in our vector
  for (chr in chromosomes) {
    # Get all genes for this chromosome from TCGA annotation
    chr_genes <- tcga_annotation[seqnames(tcga_annotation) == chr]
    
    # Count genes by type for this chromosome
    gene_types <- table(chr_genes$gene_type)
    
    # For each gene type found in this chromosome
    for (gene_type in names(gene_types)) {
      # Get genes of this specific type
      type_specific_genes <- chr_genes[chr_genes$gene_type == gene_type]$gene_id
      total_genes <- length(type_specific_genes)
      
      # Get genes present in normal network for this chromosome
      normal_network_genes <- unique(unlist(strsplit(
        V(normal_networks[[chr]])$gene_ids, ";"
      )))
      present_normal <- sum(type_specific_genes %in% normal_network_genes)
      
      # Get genes present in TNBC network for this chromosome
      tnbc_network_genes <- unique(unlist(strsplit(
        V(tnbc_networks[[chr]])$gene_ids, ";"
      )))
      present_tnbc <- sum(type_specific_genes %in% tnbc_network_genes)
      
      # Calculate genes present in both networks
      genes_in_normal <- type_specific_genes[type_specific_genes %in% normal_network_genes]
      genes_in_tnbc <- type_specific_genes[type_specific_genes %in% tnbc_network_genes]
      present_both <- length(intersect(genes_in_normal, genes_in_tnbc))
      
      # Store results for this chromosome and gene type
      results[[length(results) + 1]] <- data.frame(
        chr = chr,
        total_genes = total_genes,
        type = gene_type,
        present_normal = present_normal,
        present_tnbc = present_tnbc,
        present_both = present_both
      )
    }
  }
  
  # Combine all results into a single dataframe
  final_table <- do.call(rbind, results)
  
  # Reset row names for cleaner output
  rownames(final_table) <- NULL
  
  return(final_table)
}

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

# 2. nodes to granges
# Function to extract node regions from all chromosome networks
extract_all_node_regions <- function(network_list) {
  # Combine node attributes from all networks
  all_starts <- unlist(lapply(network_list, function(g) V(g)$start))
  all_ends <- unlist(lapply(network_list, function(g) V(g)$end))
  all_chrs <- unlist(lapply(network_list, function(g) V(g)$chr))
  all_names <- unlist(lapply(network_list, function(g) V(g)$name))
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = all_chrs,
    ranges = IRanges(start = all_starts, end = all_ends),
    names = all_names
  )
  
  return(gr)
}


# 3. count Interactions
count_edge_type <- function(c, g_list, ph) {
  E(g_list[[c]])$edge_type %>%
    table() %>%
    data.frame(check.names = FALSE, 
               Chromosome=c,
               Phenotype=ph) -> df
  colnames(df)[1] <- "Interaction"
  return(df)
}


# 4. jaccard index nodes
# Function to compute Jaccard Index for a given chromosome
compute_jaccard <- function(normal_graph, tnbc_graph, chrom) {
  # Extract start and end coordinates as numeric vectors
  normal_start <- as.numeric(V(normal_graph[[chrom]])$start)
  normal_end <- as.numeric(V(normal_graph[[chrom]])$end)
  
  tnbc_start <- as.numeric(V(tnbc_graph[[chrom]])$start)
  tnbc_end <- as.numeric(V(tnbc_graph[[chrom]])$end)
  
  # Remove NA values
  normal_start <- normal_start[!is.na(normal_start)]
  normal_end <- normal_end[!is.na(normal_end)]
  tnbc_start <- tnbc_start[!is.na(tnbc_start)]
  tnbc_end <- tnbc_end[!is.na(tnbc_end)]
  
  # Create GenomicRanges objects
  normal_ranges <- GRanges(seqnames = chrom, 
                           ranges = IRanges(start = normal_start, end = normal_end))
  
  tnbc_ranges <- GRanges(seqnames = chrom, 
                         ranges = IRanges(start = tnbc_start, end = tnbc_end))
  
  # Compute intersection (overlapping regions)
  intersected_ranges <- GenomicRanges::intersect(normal_ranges, tnbc_ranges)
  intersection_bp <- sum(width(intersected_ranges))  # Total bp in overlap
  
  # Compute union (total unique genomic regions)
  union_ranges <- GenomicRanges::union(normal_ranges, tnbc_ranges)
  union_bp <- sum(width(union_ranges))  # Total bp covered
  
  # Compute Jaccard Index
  jaccard_index <- intersection_bp / union_bp
  
  # Return Jaccard Index
  return(jaccard_index)
}


# 5. Jaccard index Edges
get_ji_edges_function <- function(chr, g_list1, g_list2, l) {
  g1 <- g_list1[[chr]]
  e1 <- get.edgelist(g1, names=FALSE)
  
  # set1 <- paste(paste(V(g1)$start[e1[,1]], V(g1)$end[e1[,1]], sep="-"),
  #       paste(V(g1)$start[e1[,2]], V(g1)$end[e1[,2]], sep="-"), 
  #       sep = "AND"
  # )
  set1 <- data.frame(edges=paste(paste(V(g1)$start[e1[,1]],
                                       V(g1)$start[e1[,2]], sep="-")),
                     edge_type=E(g1)$edge_type)
  
  set1 <- split(set1, set1$edge_type)
  
  g2 <- g_list2[[chr]]
  e2 <- get.edgelist(g2, names=FALSE)
  
  # set2 <- paste(paste(V(g2)$start[e2[,1]], V(g2)$end[e2[,1]], sep="-"),
  #               paste(V(g2)$start[e2[,2]], V(g2)$end[e2[,2]], sep="-"), 
  #               sep = "AND"
  # )
  set2 <- data.frame(edges=paste(paste(V(g2)$start[e2[,1]],
                                       V(g2)$start[e2[,2]], sep="-")),
                     edge_type=E(g2)$edge_type)
  
  set2 <- split(set2, set2$edge_type)
  
  # Compute Jaccard index
  unlist(lapply(l, function(x) { 
    intersection <- length(intersect(set1[[x]]$edges, set2[[x]]$edges))
    union <- length(union(set1[[x]]$edges, set2[[x]]$edges))
    jaccard_index <- intersection / union
  })) -> JIE
  
  data.frame(Chromosome=chr, Interactions=l, JI_Edges=JIE)
  
}

# 6. Degree and weighted degree
# Function to calculate DEGREE
compute_degree_df <- function(graph_list, phenotype) {
  degree_data <- do.call(rbind, lapply(names(graph_list), function(chr) {
    g <- graph_list[[chr]]
    
    data.frame(
      node = V(g)$name,
      chr = chr,
      w_degree_zscore = strength(g, weights = E(g)$zscore),
      w_degree_count = strength(g, weights = E(g)$count),
      degree = degree(g, mode = "all"),
      node_type = V(g)$node_type,
      phenotype = phenotype  # "Normal" or "Cancer"
    )
  }))
  
  return(degree_data)
}

# 7. zscore distribution plot
plot_zscore_distribution_function <- function(df, chromosome, outdir_plots) {
  df_filtered <- df %>% filter(chromosome == !!chromosome)
  
  # Define percentiles
  prob <- seq(0, 1, 0.0001)
  
  # Compute quantiles for each phenotype
  dat <- df_filtered %>% 
    group_by(phenotype) %>% 
    summarise(
      zscore = list(quantile(zscore, probs = prob, na.rm = TRUE)),
      Percentile = list(prob * 100)
    ) %>% 
    unnest(cols = c(zscore, Percentile))
  
  # Determine y-axis limits dynamically
  y_min <- min(dat$zscore, na.rm = TRUE)
  y_max <- max(dat$zscore, na.rm = TRUE)
  
  # Define colors for phenotypes
  phenotype_colors <- c("Normal" = "#008080", "TNBC" = "darkred")
  
  # Create percentile-based plots
  p1 <- ggplot(dat, aes(Percentile, zscore, colour = phenotype)) +
    geom_line() +
    scale_color_manual(values = phenotype_colors) +
    theme_bw() +
    ylim(y_min, y_max) + ylab("Z-Score") +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
    ggtitle("Lower")
  
  p2 <- ggplot(dat, aes(Percentile, zscore, colour = phenotype)) +
    geom_line() +
    scale_color_manual(values = phenotype_colors) +
    theme_bw() +
    ylim(y_min, y_max) + ylab("Z-Score") +
    scale_x_continuous(limits = c(90, 100), breaks = seq(90, 100, by = 1)) +
    ggtitle("Upper")
  
  # Create the boxplot (without outliers)
  p3 <- ggplot(df_filtered, aes(x = phenotype, y = zscore, fill = phenotype)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_manual(values = phenotype_colors) +
    theme_bw() +
    ylab("Z-Score (outliers = FALSE)") + xlab("Phenotype") +
    ggtitle(paste("Z-score Distribution", chromosome)) +
    theme(legend.position = "none")
  
  # Arrange all three plots
  combined_plot <- p3 + (p1 + p2 +
                           plot_layout(axes = "collect", guides = "collect")  
                         & theme(legend.position = 'bottom')) +
    plot_layout(widths = c(2.5, 5))
  
  # Save the plot
  svg(file.path(outdir_plots, paste0("zscore_distribution_", chromosome, ".svg")),
      width = 8.5, height = 4.3)
  print(combined_plot)
  dev.off()
}


# 8. function to obtain the intersection of nodes 
intersect_nodes_function <- function(ph1_l, ph2_l, c, nty="C",
                                     attrs=c("name", "chr", "start", "end", "genes")) {
  data.frame(name=V(ph1_l[[c]])$name, 
             node_type=V(ph1_l[[c]])$node_type,
             regions=paste(V(ph1_l[[c]])$start, V(ph1_l[[c]])$end, sep="-")
  ) %>%
    filter(node_type == !!nty) -> n1
  
  data.frame(name=V(ph2_l[[c]])$name, 
             node_type=V(ph2_l[[c]])$node_type,
             regions=paste(V(ph2_l[[c]])$start, V(ph2_l[[c]])$end, sep="-")
  ) %>%
    filter(node_type == !!nty) -> n2
  
  #list(n1, n2)
  n <- n1[n1$regions %in% intersect(n1$regions, n2$regions), "name"]
  
  #  return(n)
  
  # obtain attributes
  lapply(attrs, function(a) {
    get.vertex.attribute(graph = ph1_l[[c]], 
                         name = a, 
                         index = n)
    
  }) %>%
    do.call(cbind, .) -> out
  
  colnames(out) <- attrs
  
  return(out)
  
}

# 9. Function to find n points farthest from y=x line
find_diagonal_outliers <- function(data, x_col, y_col, n = 4) {
  # Calculate distance from y=x line
  data$distance_from_diagonal <- abs(data[[y_col]] - data[[x_col]]) / sqrt(2)
  
  # Add direction indicator (up or down from diagonal)
  data$direction <- ifelse(data[[y_col]] > data[[x_col]], "up", "down")
  
  # Get indices of n largest distances
  outlier_indices <- order(data$distance_from_diagonal, decreasing = TRUE)[1:n]
  
  return(data[outlier_indices, ])
}


# 10. calculate edge dissimilarity for nodes
node_dissimilarity_function <- function(chr) {
  normal_net <- normal_graphs[[chr]]
  cancer_net <- tnbc_graphs[[chr]]
  
  # Get common nodes between the two networks
  common_nodes <- intersect(V(normal_net)$name, V(cancer_net)$name)
  
  # Calculate Jaccard dissimilarity for each common node
  jaccard_dissimilarity <- sapply(common_nodes, function(node) {
    # Get neighbors in each network
    normal_neighbors <- neighbors(normal_net, node)$name
    cancer_neighbors <- neighbors(cancer_net, node)$name
    
    # Calculate Jaccard dissimilarity
    intersection_size <- length(intersect(normal_neighbors, cancer_neighbors))
    union_size <- length(union(normal_neighbors, cancer_neighbors))
    
    # Handle potential division by zero if node has no neighbors in either network
    if (union_size == 0) return(0)
    
    # Return dissimilarity (1 - similarity)
    return(1 - intersection_size / union_size)
  })
  
  
  # Create a data frame of results
  results <- data.frame(
    node = common_nodes,
    dissimilarity = jaccard_dissimilarity,
    chromosome = chr,
    start = V(normal_net)[V(normal_net)$name %in% common_nodes]$start,
    node_type = V(normal_net)[V(normal_net)$name %in% common_nodes]$node_type,
    genes = V(normal_net)[V(normal_net)$name %in% common_nodes]$genes
  )
  
  # Sort by dissimilarity to find the most changed nodes
  results <- results[order(results$dissimilarity), ]
  
  return(results)
  
}









######### MAIN 
#
##### how many genes are present
# Generate the table
gene_coverage_table <- analyze_chromosome_genes(
  tcga_annotation,
  normal_networks = normal_graphs,
  tnbc_networks = tnbc_graphs,
  chromosomes
)

# add percentages
gene_coverage_table$percent_normal <- (gene_coverage_table$present_normal / gene_coverage_table$total_genes) * 100
gene_coverage_table$percent_tnbc <- (gene_coverage_table$present_tnbc / gene_coverage_table$total_genes) * 100
gene_coverage_table$percent_both <- (gene_coverage_table$present_both / gene_coverage_table$total_genes) * 100

# export
write.csv(gene_coverage_table, 
          file = "results/igraph_intra/tcga_annotation_covered_genes.csv", quote = F, row.names = F)

##### gene coverage heatmap
gene_coverage_table %>%
  subset(type == "protein_coding") %>%
  rename("chr" = "Chromosome", 
         "present_both" = "Hi-C Nodes\n(Normal ∩ TNBC)",
         "total_genes" = "GRCh38\nAnnotation") %>%
  dplyr::select(-type) %>%
  mutate(percent_total=100) %>%
  pivot_longer(
    cols = c("Hi-C Nodes\n(Normal ∩ TNBC)", "GRCh38\nAnnotation"),
    names_to = "Nodes",
    values_to = "Genes"
  ) %>%
  # Add corresponding percentages
  mutate(`% of\nGenes` = case_when(
    Nodes == "GRCh38\nAnnotation" ~ percent_total,
    Nodes == "Hi-C Nodes\n(Normal ∩ TNBC)" ~ percent_both
  ))  %>%
  dplyr::select(Chromosome, Nodes, Genes, `% of\nGenes`) %>%
  mutate(Chromosome = factor(Chromosome, levels=rev(chromosomes)),
         Nodes = factor(Nodes, levels=c("Hi-C Nodes\n(Normal ∩ TNBC)",
                                        "GRCh38\nAnnotation"))) %>%
  ggplot(aes(Nodes, Chromosome, fill=`% of\nGenes`)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient(low = "grey", high = "steelblue") +
  geom_text(aes(label = scales::comma(Genes)), color = "white", size = 4) +
  theme_minimal() + xlab("\nNumber of Genes") + guides(fill="none") -> pheat_total


# ideogram of genome and coverage by the hic networks

# Filter protein_coding genes
protein_coding_genes <- tcga_annotation[mcols(tcga_annotation)$gene_type == "protein_coding"]

# continuous regions normal and tnbc
normal_node_regions <- extract_all_node_regions(normal_graphs)
tnbc_node_regions <- extract_all_node_regions(tnbc_graphs)

# genome ideogram
as.ggplot(expression(
  
  pp <- getDefaultPlotParams(plot.type=2),
  pp$data1height <- 350,
  #pp$topmargin <- 200
  
  kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes=chromosomes, plot.params = pp),
  #kpAddBaseNumbers(kp, chromosomes="chr1")
  #kpDataBackground(kp, data.panel = "ideogram", color = "#FFFFFFAA")
  kp <- kpPlotDensity(kp, data.panel = 1, protein_coding_genes,
                      window.size = 1000000, col="steelblue", border=NA),
  
  
  kp <- kpPlotRegions(kp, data=GenomicRanges::reduce(normal_node_regions),
                      r0=0, r1=0.45, data.panel = 2, col = "#008080", border=NA, ),
  kp <- kpPlotRegions(kp, data=GenomicRanges::reduce(tnbc_node_regions),
                      r0=0.6, r1=1, data.panel = 2, col = "darkred", border=NA),
  
  
  legend("bottomright", legend = c("Coding Genes Density", "Normal Breast Hi-C Nodes", "TNBC Hi-C Nodes"),
         col = c("steelblue", "#008080", "darkred"), #ncol=3,
         pch = c(NA, 15, 15), lty = c(1, NA, NA), lwd = 2.5, bty = "n",
         pt.cex = 1.5)
)) -> pkar


svg(file.path(outdir_plots, "supp1_heatmap_ideogram_nodes.svg"), width = 11, height = 8)
pheat_total + pkar +
  plot_layout(widths = c(2, 8), axes="collect")
dev.off()


###########  barplot how many nodes (categorical variable is type of node by content)
# Normal nodes data
lapply(normal_graphs, function(x) { table(V(x)$node_type)}) %>%
  do.call(rbind, .) %>% data.frame(check.names = FALSE) %>%
  mutate(Phenotype="Normal Breast",
         Chromosome=rownames(.)) %>%
  melt() %>% 
  rename("variable" = "Node Type", "value" = "Number of Nodes") %>%
  mutate(Chromosome = factor(Chromosome, levels=chromosomes),
         `Node Type` = factor(`Node Type`, 
                              levels=c("R",
                                       "C",
                                       "N"))) -> normal_bar_data


# TNBC nodes data
lapply(tnbc_graphs, function(x) { table(V(x)$node_type)}) %>%
  do.call(rbind, .) %>% data.frame(check.names = FALSE) %>%
  mutate(Phenotype="TNBC",
         Chromosome=rownames(.)) %>%
  melt() %>% 
  rename("variable" = "Node Type", "value" = "Number of Nodes") %>%
  mutate(Chromosome = factor(Chromosome, levels=chromosomes),
         `Node Type` = factor(`Node Type`, 
                              levels=c("R",
                                       "C",
                                       "N"))) -> tnbc_bar_data


# rbind both
both_bar_data <- rbind(normal_bar_data, tnbc_bar_data)

ggplot(both_bar_data, aes(x = Phenotype, y = `Number of Nodes`, fill = `Node Type`)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Chromosome) +
  theme_bw() +
  theme(
     axis.text.x = element_text(angle = 45, hjust = 1),
     legend.position = "top"
  ) +
  ggtitle("Genomic Feature within Hi-C Nodes") + xlab("\nChromosme") +
  facet_grid(~ Chromosome, labeller = labeller(Chromosome = function(x) gsub("chr", "", x))) +
  scale_fill_brewer(palette = "Set2", name = "", 
                    labels = c("R"="miRNA/lncRNA", "C"="Coding Gene", "N"= "ncDNA")) -> barp_node_n

svg(file.path(outdir_plots, "supp2_barplot_nodes.svg"), width = 11, height = 5.5)
barp_node_n
dev.off()

###########  barplot how many edges (categorical variable is type of interaction)
#Extract number of edges per chromosome for normal and TNBC networks
normal_edges <- sapply(normal_graphs, function(g) gsize(g))
tnbc_edges <- sapply(tnbc_graphs, function(g) gsize(g))

# Create a data frame for plotting
df <- data.frame(
  Chromosome = names(normal_edges),
  Normal = normal_edges,
  TNBC = tnbc_edges
)

# Reshape data for ggplot
df_melted <- melt(df, id.vars = "Chromosome", variable.name = "Phenotype", value.name = "Edges")
df_melted$Chromosome <- factor(df_melted$Chromosome, levels=chromosomes)

# Plot
ggplot(df_melted, aes(x = Chromosome, y = Edges, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Edges per Intrachromosomal Network", x = "Chromosome", y = "Number of Edges") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.title=element_blank()) +
  scale_y_continuous(labels = comma) +  # Adds commas to y-axis numbers
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) -> p1

svg(file.path(outdir_plots, "supp_barplot_edges.svg"), width = 11, height = 5.5)
p1
dev.off()


# Jaccard index NODES
jaccard_results_nodes <- sapply(chromosomes, function(chr) {
  compute_jaccard(normal_graphs, tnbc_graphs, chr)
})

# Print results
data.frame(nodes_JI=round(jaccard_results_nodes, digits = 2),
           JI="Jaccard Index\n(Nodes)     ",
           Chromosome=factor(chromosomes, levels = chromosomes)) -> JI_df

ggplot(JI_df, aes(Chromosome, y = JI, fill=nodes_JI)) +
  geom_tile(color="white", lwd=1.25, linetype=1) +
  scale_fill_gradient(low = "grey", high = "grey50") +
  geom_text(aes(label = nodes_JI), color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(color="black", size=7)) +
  guides(fill="none") + ylab("") +
  labs(x="Chromosome") -> p3


# CHANGE by chromosome
df$log2FC <- log2(df$TNBC / df$Normal)
df$Chromosome <- factor(df$Chromosome, levels=chromosomes)

ggplot(df, aes(x = Chromosome, y = log2FC)) +
  geom_bar(stat = "identity", fill = "blue2") +
  theme_minimal() +
  labs(x = "Chromosome", y = "TNBC / Normal\n(Log2 FC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = seq(-0.25,-0.75,-0.25), linetype="dotted", lwd=0.5, col="grey") -> p2

svg(file.path(outdir_plots, "fig1_edges_nodes_jaccard.svg"), width = 9, height = 7)

p1 / plot_spacer() / p2 / plot_spacer() / p3 +
  plot_layout(axes = 'collect', axis_titles = 'collect', heights = c(5, -0.3, 1.15, -0.2, 0.65))

dev.off()



########## Number of edges by type (supp)
rbind(
  do.call(rbind, lapply(chromosomes, count_edge_type, g_list=normal_graphs, ph="Normal")),
  do.call(rbind, lapply(chromosomes, count_edge_type, g_list=tnbc_graphs, ph="TNBC"))
) -> edge_types_df

edge_types_df$Chromosome <- factor(edge_types_df$Chromosome, levels = chromosomes)

dplyr::recode(edge_types_df$Interaction,
              "C-C" = "Coding - Coding", "C-R" = "Coding - miRNA/lncRNA",
              "C-N" = "Coding - ncDNA", "R-R" = "miRNA/lncRNA - miRNA/lncRNA", 
              "N-R" = "ncDNA - miRNA/lncRNA", "N-N" = "ncDNA - ncDNA"
              ) %>%
  factor(levels = c("Coding - Coding","Coding - miRNA/lncRNA",
                    "Coding - ncDNA", "miRNA/lncRNA - miRNA/lncRNA", 
                    "ncDNA - miRNA/lncRNA", "ncDNA - ncDNA")) -> edge_types_df$edges_hic

svg(file.path(outdir_plots, "supp_edges_types_bars_long.svg"), width = 19, height = 12)

ggplot(edge_types_df, aes(fill=Phenotype, y=Freq, x=edges_hic)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
   xlab("\nHi-C Interactions (Edges)") + ggtitle("Edge Types by Chromosome") +
  facet_wrap(~Chromosome, ncol = 6, scales = "free_y")

dev.off()


########## Jaccard of Edges by type of Edge

# chr="chr21"
# g_list1=normal_graphs
# g_list2=tnbc_graphs
# 
# 

lapply(chromosomes, get_ji_edges_function,
       g_list1 = normal_graphs, 
       g_list2 = tnbc_graphs, 
       l=rev(c("C-C", "C-R", "C-N", "R-R", "N-R", "N-N"))
) %>%
  do.call(rbind, .) -> jaccard_results_edges

jaccard_results_edges$Chromosome <- factor(jaccard_results_edges$Chromosome, levels = chromosomes)

jaccard_results_edges$Interactions <- factor(jaccard_results_edges$Interactions, 
                                             levels = rev(c("C-C", "C-R", "C-N", "R-R", "N-R", "N-N")))

svg(file.path(outdir_plots, "fig1_hmap_edges_jaccard.svg"), width = 9, height = 3.5)

ggplot(jaccard_results_edges, aes(Chromosome, Interactions, fill=JI_Edges)) +
  geom_tile(color = "white",
            lwd = 0.25,
            linetype = 1) +
  scale_fill_viridis_c(name="Jaccard\nIndex", option = "plasma") +
  #scale_fill_viridis_c(name="Jaccard\nIndex", option = "plasma", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "left") +
  ylab("Interaction (edge) Type") +
  ggtitle("Pairwise Jaccard Similarity of Hi-C Interactions Between Normal and TNBC")

dev.off()


############ degree and hubs
# Compute degrees for normal and cancer networks
normal_degree_df <- compute_degree_df(normal_graphs, "Normal")
cancer_degree_df <- compute_degree_df(tnbc_graphs, "TNBC")

# Combine both datasets
degree_df <- rbind(normal_degree_df, cancer_degree_df)
degree_df$node_type <- factor(degree_df$node_type, levels=c("C", "R", "N"))
degree_df$chr <- factor(degree_df$chr, levels=chromosomes)


svg(file.path(outdir_plots, "fig1_vln_node_degree.svg"), width = 12, height = 8)

ggplot(degree_df, aes(x = node_type, y = degree, fill = phenotype)) +
  geom_violin(position = position_dodge(width = 0.75), trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, position = position_dodge(0.75), 
               outlier.alpha = 0.75, outlier.color = "darkgrey", outlier.fill = NA, outlier.shape = 1) +  # Adds a boxplot inside the violin
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) +  
#  scale_fill_manual(values = c("Normal" = "#1b9e77", "TNBC" = "#d95f02")) + 
  scale_x_discrete(labels=c("C"="Coding", "R"="miRNA/\nlncRNA", "N"="ncDNA")) +
  theme_minimal() +
  theme(legend.position = "top") +
#  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Node Degree Distribution by Node Type",
       x = "\nNode Type",
       y = "Degree",
       fill = "Phenotype") +
  facet_wrap(~chr, scales="free_y")

dev.off()

######### Z score distribution
# plot the distribution (or density of zscore for the edges, and edge types)
#
# normal
do.call(rbind, lapply(chromosomes, function(c) {
  g <- normal_graphs[[c]]
  
  data.frame(
             zscore = E(g)$zscore,
             chromosome = c,
             edge_type = E(g)$edge_type,
             distance=E(g)$distance,
             phenotype="Normal")
})) -> zedges_normal 

# tnbc
do.call(rbind, lapply(chromosomes, function(c) {
  g <- tnbc_graphs[[c]]
  
  data.frame(
    zscore = E(g)$zscore,
    chromosome = c,
    edge_type = E(g)$edge_type,
    distance=E(g)$distance,
    phenotype="TNBC")
})) -> zedges_tnbc 

zedges_all <- rbind(zedges_normal, zedges_tnbc)
rm(zedges_normal, zedges_tnbc)
zedges_all$chromosome <- factor(zedges_all$chromosome, levels=chromosomes)
zedges_all$edge_type <- factor(zedges_all$edge_type,
                               levels = rev(c("C-C", "C-R", "C-N", "R-R", "N-R", "N-N")))


zedges_coding <- zedges_all[zedges_all$edge_type %in% c("C-C", "C-R", "C-N"), ]

# plot each chromosomes zscore distribution 
lapply(chromosomes, function(chr) plot_zscore_distribution_function(zedges_coding, chr, outdir_plots))


########### WEIGHTED DEGREE
svg(file.path(outdir_plots, "fig1_vln_node_degree_WEIGHTED_ZSCORE.svg"), width = 12, height = 8)

ggplot(degree_df, aes(x = node_type, y = w_degree_zscore, fill = phenotype)) +
  geom_violin(position = position_dodge(width = 0.75), trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, position = position_dodge(0.75), 
               outlier.alpha = 0.75, outlier.color = "darkgrey", outlier.fill = NA, outlier.shape = 1) +  # Adds a boxplot inside the violin
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) +  
  #  scale_fill_manual(values = c("Normal" = "#1b9e77", "TNBC" = "#d95f02")) + 
  scale_x_discrete(labels=c("C"="Coding", "R"="miRNA/\nlncRNA", "N"="ncDNA")) +
  theme_minimal() +
  theme(legend.position = "top") +
  #  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Node Weighted Degree Distribution (Z-Score) by Node Type",
       x = "\nNode Type",
       y = "Weighted Degree (Z-Score)",
       fill = "Phenotype") +
  facet_wrap(~chr, scales="free_y")

dev.off()

############ HUBS for coding nodes
# Nodes that maintain their high degree in cancer versus those that lose connections
# a hub in normal, lost in cancer
# a hub in cancer that was not there in normal

# rbind phenos
# plot X and Y degree and weighted degree
# facet wrap by chromosome
# find the x (5?) most extreme values per quadrant and do ggtext with genes

# get common nodes for each chromosome
common_nodes_all_chr <- lapply(chromosomes, intersect_nodes_function, 
                               ph1_l = normal_graphs, ph2_l = tnbc_graphs)

names(common_nodes_all_chr) <- chromosomes

# get degrees and weighted degrees for both phenotypes
lapply(chromosomes, function(c) {
  n <- common_nodes_all_chr[[c]][,"name"]
  
  data.frame(
    common_nodes_all_chr[[c]], 
    degree_normal = degree(normal_graphs[[c]], v = n),
    w_degree_normal = strength(normal_graphs[[c]], 
                                weights = E(normal_graphs[[c]])$zscore, vids = n),
    degree_tnbc = degree(tnbc_graphs[[c]], v = n),
    w_degree_tnbc = strength(tnbc_graphs[[c]], 
                              weights = E(tnbc_graphs[[c]])$zscore, vids = n)
  )

}) -> degree_data_common_nodes

#names(degree_data_common_nodes) <- chromosomes

degree_data_common_nodes <- do.call(rbind, degree_data_common_nodes)

degree_data_common_nodes$chr <- factor(degree_data_common_nodes$chr, 
                                       levels = chromosomes)

# plot(degree_data_common_nodes[["chrX"]]$degree_normal,
#      degree_data_common_nodes[["chrX"]]$degree_tnbc)
# 
# abline(1,1, col="blue")

# abline(lm(degree_data_common_nodes[["chr21"]]$w_degree_normal ~
#    degree_data_common_nodes[["chr21"]]$w_degree_tnbc), col="blue")

# Create outlier dataset first
outlier_data <- degree_data_common_nodes %>% 
  group_by(chr) %>% 
  group_modify(~find_diagonal_outliers(.x, "w_degree_normal", "w_degree_tnbc"))

# just label coding
unlist(lapply(outlier_data$genes, function(x) {
  s <- unlist(strsplit(x, ";"))
  s <- s[s %in% tcga_annotation$gene_name[tcga_annotation$gene_type == "protein_coding"]]
  
  paste(s, collapse=";")

})) -> outlier_data$genes



svg(file.path(outdir_plots, "scatterplot_weighted_degree_small.svg"),  width = 18, height = 16)

# Create the enhanced plot
ggplot(degree_data_common_nodes, aes(x = w_degree_normal, y = w_degree_tnbc)) +
  # Add diagonal y=x line first
  geom_abline(slope = 1, intercept = 0, color = "#2c3e50", 
              linetype = "dashed", alpha = 0.5) +
  
  # Add regular points
  geom_point(data = . %>% anti_join(outlier_data, by = "name"),  # non-outlier points
             size = 1, alpha = 0.6, color = "#4a4e4d") +
  
  # Add colored outlier points
  geom_point(data = outlier_data,
             aes(color = direction),
             size = 2) +  # Make outlier points slightly larger
  
  # Add labels for outliers with colored boxes
  geom_label_repel(
    data = outlier_data,
    aes(label = genes,
        fill = direction),  # Color based on up/down
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    label.padding = unit(0.15, "lines"),
    segment.color = "#666666",
    alpha = 0.7,  # Semi-transparent boxes
    color = "black"  # Text color
  ) +
  
  # Define colors for up/down
  scale_fill_manual(values = c(
    "up" = "#ff9999",    # Light red for upregulated
    "down" = "#99ccff"   # Light blue for downregulated
  )) +
  
  # Match point colors to label colors
  scale_color_manual(values = c(
    "up" = "#ff4444",    # Darker red for points
    "down" = "#4477ff"   # Darker blue for points
  )) +
  
  # Customize theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#eeeeee"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none"  # Hide the legend
  ) +
  
  # Add labels
  labs(
    x = "Weighted Degree Z-score (Normal)",
    y = "Weighted Degree Z-score (TNBC)",
    title = "Changes in Hi-C Interaction Patterns between Normal and TNBC Samples"
  ) +
  
  # Facet by chromosome
  facet_wrap(~chr, scales = "free")

dev.off()


### examples
# PTEN and PXDC1

pxdc1_node = degree_data_common_nodes[grepl("PXDC1", degree_data_common_nodes$genes), ]
pxdc1_node = degree_data_common_nodes[grepl("PXDC1", degree_data_common_nodes$genes), ]



############ SIMILARITY INDEX


############ DISTANCE

calculate_distance_zscores <- function(graph) {
  # Check if the graph has a distance attribute
  if (!"distance" %in% edge_attr_names(graph)) {
    stop("Graph must have a 'distance' edge attribute")
  }
  
  # Get the distance values
  distances <- E(graph)$distance
  
  # Calculate z-scores
  distance_mean <- mean(distances, na.rm = TRUE)
  distance_sd <- sd(distances, na.rm = TRUE)
  zscores <- (distances - distance_mean) / distance_sd
  
  # Add z-scores as a new edge attribute
  E(graph)$distance_zscore <- zscores
  
  # Return the modified graph
  return(graph)
}

#hist(E(calculate_distance_zscores(normal_graphs$chr1))$distance_zscore)

normal_graphs <- lapply(normal_graphs, calculate_distance_zscores)
tnbc_graphs <- lapply(tnbc_graphs, calculate_distance_zscores)

# Function get the distances
get_dists <- function(c, glist, ph) {
  g = glist[[c]]
  
  data.frame(dist=E(g)$distance,
             zdist=E(g)$distance_zscore,
             edge_type=E(g)$edge_type,
             chromosome=c,
             phenotype=ph)
}

# get them all
rbind(
  do.call(rbind, lapply(chromosomes, get_dists, glist=normal_graphs, ph="Normal")),
  do.call(rbind, lapply(chromosomes, get_dists, glist=tnbc_graphs, ph="TNBC"))
) -> distances_df


distances_df_coding <- distances_df[grepl("C", distances_df$edge_type),  ]

distances_df_coding$chromosome <- factor(distances_df_coding$chromosome, 
                                         levels = chromosomes)

distances_df_coding$intype <- paste(distances_df_coding$phenotype, distances_df_coding$edge_type)
distances_df_coding$intype <- factor(distances_df_coding$intype, 
                                     levels = c("Normal C-C", "Normal C-R", "Normal C-N",
                                                "TNBC C-C", "TNBC C-R", "TNBC C-N"))



svg(file.path(outdir_plots, "fig1_hist_edge_GENOMICDIST.svg"), width = 10, height = 10)

# ggplot(distances_df_coding, aes(x = dist/1000000, fill = phenotype)) +
#   geom_histogram(bins = 30, alpha = 0.45, position = "identity") +
#   labs(title = "Edge Distance by Chromosome",
#        x = "\nGenomic Distance (Mb)",
#        y = "Count") +
#   theme_minimal() +
# #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) + 
#   #theme_minimal() + 
#   facet_wrap(~chromosome, scales="free")

ggplot(distances_df_coding, aes(x = dist/1000000, fill = phenotype, color = phenotype)) +
  geom_histogram(bins = 30, alpha = 0.3, position = "identity") +  # Set alpha for transparency
  labs(title = "Edge Distance by Chromosome",
       x = "\nGenomic Distance (Mb)",
       y = "Count",
       fill = "Phenotype",
       color = "Phenotype") + 
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) +  # Custom fill colors
  scale_color_manual(values = c("Normal" = "#008080", "TNBC" = "darkred")) +  # Matching outline colors
  facet_wrap(~chromosome, scales = "free")

dev.off()

# density for supplementary
svg(file.path(outdir_plots, "fig1_hist_edge_GENOMICDIST_TYPE.svg"),  width = 10.5, height = 7)
ggplot(distances_df_coding, aes(x = dist/1000000, color = intype, linetype = intype)) +
  geom_density(size = 1, lwd=1) +
  scale_linetype_manual(values = c("Normal C-C" = "solid", "Normal C-R" = "dashed", "Normal C-N" = "dotted",
                                   "TNBC C-C" = "solid", "TNBC C-R" = "dashed", "TNBC C-N" = "dotted")) +
  scale_color_manual(values = c("Normal C-C" = "#008080", "Normal C-R" = "#66b2b2", "Normal C-N" = "#b2d8d8",
                                "TNBC C-C" = "darkred", "TNBC C-R" = "#b22222", "TNBC C-N" = "#e9967a")) +
  labs(title = "Edge Genomic Distance by Chromosome",
       x = "\nDistance (Mb)",
       y = "Density") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top") +
  facet_wrap(~chromosome, scales = "free", ncol=4)

dev.off()

######### set up a 1% threshold for long range interactions 
# threshold based edge filtering
# Function to filter network keeping top 1% of edges by distance
filter_network_by_distance <- function(graph, percentile = 0.99) {
  # Get all distance values
  distances <- E(graph)$distance
  
  # Calculate the threshold for top 1%
  distance_threshold <- quantile(distances, percentile)
  
  # Create subgraph keeping only edges above threshold
  # Keep all vertices but only edges that meet our criteria
  subgraph_from_edges(graph, 
                 which(E(graph)$distance >= distance_threshold), 
                 delete.vertices = TRUE)
}

# Apply to all chromosomes
filt_normal_graphs <- lapply(normal_graphs, filter_network_by_distance)
filt_tnbc_graphs <- lapply(tnbc_graphs, filter_network_by_distance)


################## Similarity INDEX
# which nodes retain their interactions, which nodes lose them
# JACCARD DISSIMILARITY 
#Dissimilarity(node) = 1 - |Neighbors_normal ∩ Neighbors_cancer| / |Neighbors_normal ∪ Neighbors_cancer|

# calculate node dissimilarity
jaccard_dissimilarity_all <- lapply(chromosomes, node_dissimilarity_function)
jaccard_dissimilarity_all <- do.call(rbind, jaccard_dissimilarity_all)

jaccard_dissimilarity_all$chromosome <- factor(jaccard_dissimilarity_all$chromosome, 
                                               levels = chromosomes)


jaccard_dissimilarity_genes <- jaccard_dissimilarity_all[jaccard_dissimilarity_all$node_type %in% c("C", "R"), ]

# Filter nodes with dissimilarity index below 0.25
low_dissimilarity_nodes <- jaccard_dissimilarity_genes %>%
  filter(dissimilarity < 0.25)

# update geneID labels with current annotation names
#low_dissimilarity_nodes$genes[grepl("AC114808.1", low_dissimilarity_nodes$genes)] <- "ENSG00000235403"

# Plot
svg(file.path(outdir_plots, "box_node_dissimilarity.svg"),  width = 11, height = 6)

ggplot(jaccard_dissimilarity_genes, aes(chromosome, dissimilarity)) + 
  geom_boxplot(notch = TRUE, fill = "grey85") +
  geom_text_repel(data = low_dissimilarity_nodes, aes(label = genes), 
                  size = 3, box.padding = 0.5, col= "darkblue",
                  point.padding = 0.3, max.overlaps = Inf) +
  ggtitle("Rewiring of Hi-C Interactions by Node") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  ylab("Jaccard Dissimilarity Index\n") + 
  xlab("\nChromosome")

dev.off()




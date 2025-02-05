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
library(scales)
library(tidyr)
library(dplyr)
library(reshape2)
library(GenomicRanges)
library(karyoploteR)
library(ggplotify)

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


# 3. create node type
create_node_type <- function(graph) {
  # Get the gene types vector
  gene_types <- V(graph)$gene_types
  
  # Create node_type attribute
  V(graph)$node_type <- sapply(gene_types, function(gt) {
    # Check if gene_type is NA
    if (is.na(gt)) {
      return("No transcript in node")
    }
    
    # Check if gene_type contains "protein_coding"
    if (grepl("protein_coding", gt)) {
      return("Coding transcript in node")
    } else {
      return("miRNA/lncRNA in node")
    }
  })
  
  return(graph)
}

# Function to process all networks in a list
process_networks <- function(network_list) {
  # Apply the create_node_type function to each network in the list
  processed_networks <- lapply(network_list, create_node_type)
  return(processed_networks)
}

# alt function for real networks
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
add_ntype <- function(network_list) {
  # Apply the create_node_type function to each network in the list
  processed_networks <- lapply(network_list, really_create_node_type)
  return(processed_networks)
}


# 4.function to add edge types
add_edge_types <- function(graph) {
  edge_list <- get.edgelist(graph, names=FALSE)  # Get numeric indices
  type1 <- V(graph)$node_type[edge_list[,1]]
  type2 <- V(graph)$node_type[edge_list[,2]]
  
  edge_types <- paste(pmin(type1, type2), pmax(type1, type2), sep="-")
  E(graph)$edge_type <- edge_types
  return(graph)
}

# test <- add_edge_types(test)
# table(E(test)$edge_type)







# 4. Jaccard nodes (regions)
calculate_region_jaccard <- function(normal_graph, tnbc_graph) {
  # Extract node attributes to create GRanges objects
  normal_regions <- GRanges(
    seqnames = V(normal_graph)$chr,
    ranges = IRanges(
      start = V(normal_graph)$start,
      end = V(normal_graph)$end
    )
  )
  
  tnbc_regions <- GRanges(
    seqnames = V(tnbc_graph)$chr,
    ranges = IRanges(
      start = V(tnbc_graph)$start,
      end = V(tnbc_graph)$end
    )
  )
  
  # Calculate intersection
  intersection_width <- sum(width(intersect(normal_regions, tnbc_regions)))
  
  # Calculate union
  union_width <- sum(width(union(normal_regions, tnbc_regions)))
  
  # Calculate Jaccard index
  jaccard_index <- intersection_width / union_width
  
  return(jaccard_index)
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


write.csv(gene_coverage_table, 
          file = "results/igraph_intra/tcga_annotation_covered_genes.csv", quote = F, row.names = F)

# get the heatmaps 
gene_coverage_table %>%
  subset(type == "protein_coding") %>%
  rename("chr" = "Chromosome", 
         "present_normal" = "Normal",
         "present_tnbc" = "TNBC",
         "present_both" = "Overlap") %>%
  dplyr::select(-type) %>%
  mutate(percent_total=100) -> plot_data

# total genes heatmap
plot_data %>%
  dplyr::select(Chromosome, total_genes, percent_total) %>%
  mutate(all_genes = "Total Number of\ngenes (100%)",
         Chromosome = factor(Chromosome, levels=rev(chromosomes))) %>%
  #ggplot(aes(all_genes, Chromosome, fill=percent_total)) +
  ggplot(aes(all_genes, Chromosome)) +
  geom_tile(color = "white",
          lwd = 1.5,
          linetype = 1, fill="#56B1F7") +
  geom_text(aes(label = scales::comma(total_genes)), color = "white", size = 4) +
  theme_minimal() + xlab("") -> pheat_total


# genes present in normal or tnbc hic networks heatmap
plot_data %>%
  pivot_longer(
    cols = c(Normal, TNBC, Overlap),
    names_to = "Nodes",
    values_to = "Genes"
  ) %>%
  # Add corresponding percentages
  mutate(`% of\nGenes` = case_when(
    Nodes == "Normal" ~ percent_normal,
    Nodes == "TNBC" ~ percent_tnbc,
    Nodes == "Overlap" ~ percent_both
  )) %>%
  dplyr::select(Chromosome, Nodes, Genes, `% of\nGenes`) %>%
  mutate(Chromosome = factor(Chromosome, levels=rev(chromosomes)),
         Nodes = factor(Nodes, levels=c("Normal", "TNBC", "Overlap"))) %>%
  ggplot(aes(Nodes, Chromosome, fill=`% of\nGenes`)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(label = scales::comma(Genes)), color = "white", size = 4) +
  theme_minimal() + xlab("Hi-C Nodes") -> pheat_genes


# hp <- pheat_genes + pheat_total +
#   plot_layout(widths = c(3, 1), axes="collect")


##### regions ideogram
# ideogram of genome and coverage by the hic networks

# Filter protein_coding genes
protein_coding_genes <- tcga_annotation[mcols(tcga_annotation)$gene_type == "protein_coding"]

# continuous regions normal and tnbc
normal_node_regions <- extract_all_node_regions(normal_graphs)
tnbc_node_regions <- extract_all_node_regions(tnbc_graphs)

#svg(file.path(outdir_plots, "supp1_ideogram_hic_nodes.svg"), width = 10, height = 13)
# Prepare gene density track

as.ggplot(expression(

pp <- getDefaultPlotParams(plot.type=2),
pp$data1height <- 350,
#pp$topmargin <- 200

kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes=chromosomes, plot.params = pp),
#kpAddBaseNumbers(kp, chromosomes="chr1")
#kpDataBackground(kp, data.panel = "ideogram", color = "#FFFFFFAA")
kp <- kpPlotDensity(kp, data.panel = 1, protein_coding_genes,
                    window.size = 1000000, col="cornflowerblue", border=NA),


kp <- kpPlotRegions(kp, data=GenomicRanges::reduce(normal_node_regions),
              r0=0, r1=0.45, data.panel = 2, col = "#008080", border=NA, ),
kp <- kpPlotRegions(kp, data=GenomicRanges::reduce(tnbc_node_regions),
              r0=0.6, r1=1, data.panel = 2, col = "darkred", border=NA),


legend("bottomright", legend = c("Coding Genes Density", "Normal Breast Hi-C Nodes", "TNBC Hi-C Nodes"),
       col = c("cornflowerblue", "#008080", "darkred"), #ncol=3,
       pch = c(NA, 15, 15), lty = c(1, NA, NA), lwd = 2.5, bty = "n",
       pt.cex = 1.5)
)) -> pkar

#dev.off()

svg(file.path(outdir_plots, "supp1_heatmap_ideogram_nodes.svg"), width = 11, height = 8)
pheat_genes + pheat_total + pkar +
  plot_layout(widths = c(3, 1, 9), axes="collect")
dev.off()


# barplot how many nodes (categorical variable is type of node by content)
# Normal nodes data
process_networks(normal_graphs) %>% 
  lapply(., function(x) { table(V(x)$node_type)}) %>%
  do.call(rbind, .) %>% data.frame(check.names = FALSE) %>%
  mutate(Phenotype="Normal Breast",
         Chromosome=rownames(.)) %>%
  melt() %>% 
  rename("variable" = "Node Type", "value" = "Number of Nodes") %>%
  mutate(Chromosome = factor(Chromosome, levels=chromosomes),
         `Node Type` = factor(`Node Type`, 
                              levels=c("miRNA/lncRNA in node",
                                       "Coding transcript in node",
                                       "No transcript in node"))) -> normal_bar_data

# 
# ggplot(normal_bar_data, aes(x = Chromosome, y = `Number of Nodes`, fill = `Node Type`)) +
#   geom_bar(stat = "identity") +
#   labs(
#     title = "Number of Hi-C Nodes by Chromosome",
#     x = "Chromosome",
#     y = "Number of Nodes"
#   ) +
#   theme_minimal() +
#   # Rotate x-axis labels if needed
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right"
#   ) +
#   # Use a colorblind-friendly palette
#   scale_fill_brewer(palette = "Set2")

# TNBC nodes data
process_networks(tnbc_graphs) %>% 
  lapply(., function(x) { table(V(x)$node_type)}) %>%
  do.call(rbind, .) %>% data.frame(check.names = FALSE) %>%
  mutate(Phenotype="TNBC",
         Chromosome=rownames(.)) %>%
  melt() %>% 
  rename("variable" = "Node Type", "value" = "Number of Nodes") %>%
  mutate(Chromosome = factor(Chromosome, levels=chromosomes),
         `Node Type` = factor(`Node Type`, 
                              levels=c("miRNA/lncRNA in node",
                                       "Coding transcript in node",
                                       "No transcript in node"))) -> tnbc_bar_data


# rbind both
both_bar_data <- rbind(normal_bar_data, tnbc_bar_data)

ggplot(both_bar_data, aes(x = Phenotype, y = `Number of Nodes`, fill = `Node Type`)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Chromosome) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) + ggtitle("Intrachromosomal Networks") +
  facet_grid(~ Chromosome, labeller = labeller(Chromosome = function(x) gsub("chr", "", x))) +
  scale_fill_brewer(palette = "Set2") -> barp_node_n

svg(file.path(outdir_plots, "supp2_barplot_nodes.svg"), width = 12, height = 5.5)
barp_node_n
dev.off()



########### nodes jaccard index vector for each chromosome

########## interaction types
# we actually need the node types
normal_graphs <- add_ntype(normal_graphs)
tnbc_graphs <- add_ntype(tnbc_graphs)

# add edge type as attribute



# Process all graphs
normal_graphs <- lapply(normal_graphs, add_edge_types)
tnbc_graphs <- lapply(tnbc_graphs, add_edge_types)




# SAVE THE ATTRIBUTED NETS
#saveRDS(edge_types_normal, file="edge_types_normal.Rds")


#################################################################################
#################################################################################

# Set seed for reproducibility
set.seed(123)

# Define number of nodes per type
num_ncDNA <- 500
num_gene <- 300
num_ncRNA <- 200

# Generate synthetic node degree data
create_data <- function(condition, chromosome) {
  data.frame(
    node_type = rep(c("ncDNA", "gene", "ncRNA"), times = c(num_ncDNA, num_gene, num_ncRNA)),
    degree = c(
      rpois(num_ncDNA, lambda = ifelse(condition == "TNBC", 15, 13)),  # ncDNA nodes
      rpois(num_gene, lambda = ifelse(condition == "TNBC", 20, 18)),    # Gene nodes
      rpois(num_ncRNA, lambda = ifelse(condition == "TNBC", 10, 9))     # ncRNA nodes
    ),
    phenotype = condition,
    chromosome = chromosome
  )
}

# Create datasets for Normal and TNBC across chromosomes
chromosomes <- paste0("chr", c(1:22, "X"))
data_list <- lapply(chromosomes, function(chr) {
  rbind(create_data("Normal", chr), create_data("TNBC", chr))
})

data <- do.call(rbind, data_list)

data$chromosome <- factor(data$chromosome, levels=chromosomes)


# Plot
p <- ggplot(data, aes(x = node_type, y = degree, fill = phenotype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~chromosome, scales = "free") +
  theme_minimal() +
  labs(title = "Node Degree Distribution by Node Type, Condition, and Chromosome",
       x = "Node Type", y = "Node Degree") +
  scale_fill_manual(values = c("Normal" = "#008080", "TNBC" = "darkred"))

print(p)


######
library(ggplot2)
library(igraph)
library(stringi)

# Set seed for reproducibility
set.seed(123)

# Generate random Ensembl gene IDs
random_gene_id <- function() {
  paste0("ENSG", stri_rand_strings(1, 11, pattern = "[0-9]"))
}

lost_hub <- random_gene_id()
connections_normal <- data.frame(
  from = rep(lost_hub, 10),
  to = replicate(10, random_gene_id()),
  weight = runif(10, 0.1, 1)
)

g_normal <- graph_from_data_frame(connections_normal, directed = FALSE)
E(g_normal)$color <- scales::col_numeric("Blues", domain = NULL)(E(g_normal)$weight)

# Plot the lost hub in Normal
plot(g_normal, edge.width = E(g_normal)$weight * 5, edge.color = E(g_normal)$color, 
     main = paste("Lost Hub in Normal: ", lost_hub))

# Create a reduced connection network in Cancer
connections_cancer <- data.frame(
  from = rep(lost_hub, 5),  # Fewer connections
  to = replicate(5, random_gene_id()),
  weight = runif(5, 0.1, 1)
)

g_cancer <- graph_from_data_frame(connections_cancer, directed = FALSE)
E(g_cancer)$color <- scales::col_numeric("Reds", domain = NULL)(E(g_cancer)$weight)

# Plot the lost hub in Cancer
plot(g_cancer, edge.width = E(g_cancer)$weight * 5, edge.color = E(g_cancer)$color, 
     main = paste("Lost Hub in Cancer: ", lost_hub))





#################################################################################
#################################################################################
#################################################################################
## degree histogram
lapply(normal_graphs, function(x) {
  degree(x)
}) -> normal_degrees 

names(normal_degrees) <- chromosomes


lapply(tnbc_graphs, function(x) {
  degree(x)
}) -> tnbc_degrees 

names(tnbc_degrees) <- chromosomes

lapply(chromosomes, function(c) {
  rbind(
    data.frame(degree=normal_degrees[[c]], Phenotype="Normal", Chromosome=c),
    data.frame(degree=tnbc_degrees[[c]], Phenotype="TNBC", Chromosome=c)
  )
}) %>%
  do.call(rbind, .) -> degrees_df

degrees_df$Chromosome <- factor(degrees_df$Chromosome, levels=chromosomes)

# 
# rbind(data.frame(degree=normal_degrees$chr1, Phenotype="Normal", Chromosome="chr1"),
#       data.frame(degree=tnbc_degrees$chr1, Phenotype="TNBC", Chromosome="chr1")
# )
 
ggplot(degrees_df, aes(x=degree)) + geom_histogram(fill=NA, col="black") +
  facet_grid(Phenotype~Chromosome, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## 



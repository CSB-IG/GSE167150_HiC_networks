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

######### POINTS ISSUE
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

test21 <- zedges_coding[zedges_coding$chromosome == "chr21", ]
# boxplot(test21$zscore~test21$phenotype, horizontal=TRUE)
# boxplot(test21$zscore~test21$phenotype, horizontal=TRUE, outline=F, notch=T)

####### plot to do lower and higher percentiles
# Define parameters
low_percentile <- c(0, 10)     # Adjust range for the lower tail
high_percentile <- c(90, 100)  # Adjust range for the upper tail

# Define colors for phenotypes
phenotype_colors <- c("Normal" = "#008080", "TNBC" = "darkred")

# Define percentiles
prob <- seq(0, 1, 0.0001)

# Compute quantiles for each phenotype
dat <- test21 %>% 
  group_by(phenotype) %>% 
  summarise(
    zscore = list(quantile(zscore, probs = prob, na.rm = TRUE)),
    Percentile = list(prob * 100)
  ) %>% 
  unnest(cols = c(zscore, Percentile))

# Determine y-axis limits dynamically
y_min <- min(dat$zscore, na.rm = TRUE)
y_max <- max(dat$zscore, na.rm = TRUE)

# Create the percentile-based plots
p1 <- ggplot(dat, aes(Percentile, zscore, colour = phenotype)) +
  geom_line() +
  scale_color_manual(values = phenotype_colors) +
  theme_bw() + #theme(legend.position = "top") +
  ylim(y_min, y_max) + ylab("Z-Score") +
  scale_x_continuous(limits = low_percentile, breaks = seq(low_percentile[1], low_percentile[2], by = 1)) +
  #ggtitle(paste0("Lower ", low_percentile[2], "% Percentile"))
  ggtitle("Lower")

p2 <- ggplot(dat, aes(Percentile, zscore, colour = phenotype)) +
  geom_line() +
  scale_color_manual(values = phenotype_colors) +
  theme_bw() + #theme(legend.position = "top") +
  ylim(y_min, y_max) + ylab("Z-Score") +
  scale_x_continuous(limits = high_percentile, breaks = seq(high_percentile[1], high_percentile[2], by = 1)) +
  #ggtitle(paste0("Upper ", high_percentile[2], "% Percentile"))
  ggtitle("Upper")

# Create the boxplot (without outliers)
p3 <- ggplot(test21, aes(x = phenotype, y = zscore, fill = phenotype)) +
  geom_boxplot(outliers = FALSE) +  # Remove outliers
  scale_fill_manual(values = phenotype_colors) +
  theme_bw() +
  ylab("Z-Score\t(outliers = FALSE)") + xlab("Phenotype") +
  ggtitle(paste("Z-score Distribution", "Chr21")) +
#  labs(title = "Z-Score Distribution by Phenotype", x = "", y = "Z-Score") +
  theme(legend.position = "none")  # Hide legend for boxplot

# Arrange all three plots
#p3 + plot_spacer() + p1 + p2 +
svg(file.path(outdir_plots, paste0("zscore_distribution_", "chr21", ".svg")), 
    width = 8.5, height = 4.3)

p3 + (p1 + p2 +
  plot_layout(axes="collect", guides = "collect")  & theme(legend.position = 'bottom')) +
  plot_layout(widths = c(2.5, 5)) 

dev.off()    





# # Define colors for phenotypes
# phenotype_colors <- c("Normal" = "#008080", "TNBC" = "darkred")
# 
# # Define percentiles
# prob <- seq(0, 1, 0.0001)
# 
# # Compute quantiles for each phenotype
# dat <- test21 %>% 
#   group_by(phenotype) %>% 
#   summarise(
#     zscore = list(quantile(zscore, probs = prob, na.rm = TRUE)),
#     Percentile = list(prob * 100)
#   ) %>% 
#   unnest(cols = c(zscore, Percentile))
# 
# # Create a single plot with full percentile range
# ggplot(dat, aes(Percentile, zscore, colour = phenotype)) +
#   geom_line() +
#   scale_color_manual(values = phenotype_colors) +
#   theme_bw() +
#   labs(title = "Quantile Plot of Z-Scores", 
#        x = "Percentile", 
#        y = "Z-Score") +
# #  xlim(90, 100) +
#   theme(legend.title = element_blank())











# WEIGHTED DEGREE
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
  labs(title = "Node Degree Distribution by Node Type",
       x = "\nNode Type",
       y = "Weighted Degree (zscore)",
       fill = "Phenotype") +
  facet_wrap(~chr, scales="free_y")

dev.off()

############ HUBS for coding nodes

# a hub in normal, lost in cancer

# a hub in cancer that was not there in normal








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

##### scatterplots
library(ggplot2)
library(dplyr)
library(tidyr)

# Set random seed for reproducibility
set.seed(123)

# Define the number of genes per chromosome
genes_per_chr <- c(1900, 1245, 1000, 751, 885, 1047, 919, 683, 778, 727, 
                   1310, 1033, 321, 610, 596, 854, 1183, 269, 1470, 546, 234, 444, 760)
chromosomes <- c(1:22, "X")  # Include Chromosome X

# Function to generate clustering coefficient data
generate_data <- function(chr, n) {
  x <- runif(n, 0, 1)  # Simulating clustering coefficients in normal networks
  y <- x + rnorm(n, 0, 0.1)  # TNBC sample with variation around diagonal
  
  # Convert chromosome number safely
  chr_num <- suppressWarnings(as.numeric(chr))  # Will be NA for "X"
  
  # Introduce chromosome-specific variations (skip if NA)
  if (!is.na(chr_num)) {
    if (chr_num %% 3 == 0) {
      y <- x + abs(rnorm(n, 0, 0.2))  # More points above the diagonal
    } else if (chr_num %% 4 == 0) {
      y <- x + rnorm(n, 0, 0.05)  # Less variation, more clustering on diagonal
    }
  }
  
  data.frame(Chromosome = paste0("Chr", chr), Normal = x, TNBC = y)
}

# Generate data for all chromosomes
mock_data <- bind_rows(mapply(generate_data, chromosomes, genes_per_chr, SIMPLIFY = FALSE))
mock_data$Chromosome <- factor(mock_data$Chromosome, 
                               levels = paste0("Chr", c(seq(1,22), "X")))

# Plot using ggplot2
ggplot(mock_data, aes(x = TNBC, y = Normal)) +
  geom_point(alpha = 0.5, color = "blue", size = 0.5) +
  facet_wrap(~Chromosome, scales = "free", ncol = 6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Clustering Coefficient (Normal)", 
       y = "Clustering Coefficient (TNBC)", 
       title = "Scatterplots of Clustering Coefficients per Chromosome") +
  theme_minimal()




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



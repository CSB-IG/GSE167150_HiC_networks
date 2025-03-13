# Function to generate search URLs for different databases
create_search_urls <- function(gene_name) {
  # Clean gene name (in case it contains multiple genes)
  genes <- strsplit(gene_name, ";")[[1]]
  
  # Create URLs for each gene
  lapply(genes, function(gene) {
    list(
      gene = gene,
      pubmed = sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s+AND+%s", 
                       gene, URLencode("triple negative breast cancer")),
      genecards = sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene),
      gwas = sprintf("https://www.ebi.ac.uk/gwas/search?query=%s", gene),
      clinicaltrials = sprintf("https://clinicaltrials.gov/search?term=%s+AND+triple+negative+breast+cancer", 
                               URLencode(gene))
    )
  })
}

# Function to get genes from our outlier analysis
get_outlier_genes <- function(data, x_col, y_col, n = 3) {
  # Use our previous outlier detection
  outliers <- data %>%
    group_by(chr) %>%
    group_modify(~find_diagonal_outliers(.x, x_col, y_col, n)) %>%
    ungroup()
  
  # Extract gene information with direction
  gene_info <- outliers %>%
    dplyr::select(chr, genes, direction, distance_from_diagonal) %>%
    arrange(desc(distance_from_diagonal))
  
  return(gene_info)
}

# Get the outlier genes
outlier_genes <- get_outlier_genes(degree_data_common_nodes, 
                                   "w_degree_normal", 
                                   "w_degree_tnbc")

# Create a dataframe with search URLs
generate_search_table <- function(gene_info) {
  # Split multiple genes and create individual rows
  gene_data <- gene_info %>%
    mutate(gene_list = strsplit(genes, ";")) %>%
    unnest(gene_list) %>%
    mutate(
      pubmed_url = sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s+AND+%s", 
                           gene_list, URLencode("triple negative breast cancer")),
      genecards_url = sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", 
                              gene_list),
      gwas_url = sprintf("https://www.ebi.ac.uk/gwas/search?query=%s", 
                         gene_list),
      clinicaltrials_url = sprintf("https://clinicaltrials.gov/search?term=%s+AND+triple+negative+breast+cancer", 
                                   URLencode(gene_list))
    )
  
  return(gene_data)
}

# Create the search table
search_table <- generate_search_table(outlier_genes)

# Write to CSV (optional)
write.csv(search_table, "tnbc_outlier_genes_search.csv", row.names = FALSE)

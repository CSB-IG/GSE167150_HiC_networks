library(igraph)
library(ggraph)
library(circlize)



# GRCh38 chromosome sizes
GRCH38_SIZES <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559,
  chr4 = 190214555, chr5 = 181538259, chr6 = 170805979,
  chr7 = 159345973, chr8 = 145138636, chr9 = 138394717,
  chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
  chr16 = 90338345, chr17 = 83257441, chr18 = 80373285,
  chr19 = 58617616, chr20 = 64444167, chr21 = 46709983,
  chr22 = 50818468, chrX = 156040895, chrY = 57227415
)


set.seed(999)
n = 1000
df = data.frame(sectors = sample(letters[1:8], n, replace = TRUE),
                x = rnorm(n), y = runif(n))


library(circlize)

circos.par("start.degree" = 90)
circos.initializeWithIdeogram(chromosome.index = "chr3")






extract_node_interactions <- function(graph, node_of_interest, ph) {
  # Find the vertex ID of the node based on its name
  node_id <- which(V(graph)$name == node_of_interest)
  
  # If node is not found, return NULL
  if (length(node_id) == 0) {
    warning("Node not found in the graph.")
    return(NULL)
  }
  
  # Get edges involving the node
  edges <- E(graph)[.inc(node_id)]
  
  # Get neighboring node IDs (first neighbors)
  neighbor_ids <- neighbors(graph, node_id)
  
  # Extract edge attributes and get the neighbor nodes
  edge_data <- data.frame(
    node1 = rep(node_of_interest, length(neighbor_ids)),  # Repeated node_of_interest
    node2 = V(graph)[neighbor_ids]$name,  # Neighbor names
    node2_genes = V(graph)[neighbor_ids]$genes,  # Genes of the neighbors (node2)
    node2_start = V(graph)[neighbor_ids]$start,
#    node2_midpoint = V(graph)[neighbor_ids]$midpoint,
    node2_end = V(graph)[neighbor_ids]$end,
    node2_type = V(graph)[neighbor_ids]$node_type,
    zscore = E(graph)[edges]$zscore,  # Edge z-score
#    distance = E(graph)[edges]$distance,  # Edge distance
    phenotype = ph
  )
  
  # just label coding
  unlist(lapply(edge_data$node2_genes, function(x) {
    s <- unlist(strsplit(x, ";"))
    s <- s[s %in% tcga_annotation$gene_name[tcga_annotation$gene_type == "protein_coding"]]
    
    paste(s, collapse=";")
    
  })) -> edge_data$node2_genes
  
  
  return(edge_data)
}


# input is gene name, common nodes, list1, list2
gene_name = "PTEN"

# get the data from common nodes object (node ID, node chr)
# PTEN normal
# PTEN tnbc
pov < - degree_data_common_nodes[grepl("PTEN", degree_data_common_nodes$genes), ]
# edit pov to show only coding

pov_node = "51353"
pov_chr = "chr10"

# call the extract node interactions function on both lists
# add a node type parameter c("C", "N", "R")
# call function and rbind

pten_interactions <- rbind(extract_node_interactions(normal_graphs$chr10, pten_node, ph="Normal"),
                           extract_node_interactions(tnbc_graphs$chr10, pten_node, ph="TNBC"))


# we could give back a list where the first element is POV and the second are ints


# use this list to plot FOR EACH PHENOTYPE (how to make comparable??)

# chromosome layout starting at midnight

# POV
# pov segment position start to end
# pov label is genes
# pov segment needs a color according to degree (degree_normal, degree_tnbc)

# NODE2
# node2 position start to end
# node2 label is genes

# chords are drawn according to ints to each node2
# zscore is parameter to plot chords colors (and width?)












############ PLOTTING FUN
# Function to create chromosome interaction diagrams with ideograms
# Function to create chromosome interaction diagrams with ideograms
create_chromosome_diagrams <- function(pov_data, interactions_data, output_file = NULL) {
  # Load required packages
  require(circlize)
  require(colorRamp2)
  require(grid)
  
  # GRCh38 chromosome sizes
  GRCH38_SIZES <- c(
    chr1 = 248956422, chr2 = 242193529, chr3 = 198295559,
    chr4 = 190214555, chr5 = 181538259, chr6 = 170805979,
    chr7 = 159345973, chr8 = 145138636, chr9 = 138394717,
    chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
    chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
    chr16 = 90338345, chr17 = 83257441, chr18 = 80373285,
    chr19 = 58617616, chr20 = 64444167, chr21 = 46709983,
    chr22 = 50818468, chrX = 156040895, chrY = 57227415
  )
  
  # Extract the POV information
  pov <- pov_data
  ints <- interactions_data
  
  # Determine which chromosome we're working with
  chr_name <- pov$chr[1]
  
  # Prepare POV node information
  pov_node <- data.frame(
    name = pov$name,
    chr = pov$chr,
    start = pov$start,
    end = pov$end,
    genes = pov$genes,
    degree_normal = pov$degree_normal,
    w_degree_normal = pov$w_degree_normal,
    degree_tnbc = pov$degree_tnbc,
    w_degree_tnbc = pov$w_degree_tnbc,
    type = "POV"
  )
  
  # Get chromosome length
  chr_length <- GRCH38_SIZES[chr_name]
  
  # Ensure all coordinates are within chromosome boundaries
  ensure_valid_coords <- function(start, end, chr_length) {
    # Make sure start is valid
    start <- max(1, min(as.numeric(start), chr_length))
    # Make sure end is valid and greater than start
    end <- max(start, min(as.numeric(end), chr_length))
    return(list(start = start, end = end))
  }
  
  # Fix POV coordinates if needed
  valid_coords <- ensure_valid_coords(pov_node$start, pov_node$end, chr_length)
  pov_node$start <- valid_coords$start
  pov_node$end <- valid_coords$end
  
  # Extract interaction nodes with validated coordinates
  int_nodes <- data.frame(
    name = character(),
    chr = character(),
    start = numeric(),
    end = numeric(),
    genes = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each interaction, validating coordinates
  for (i in 1:nrow(ints)) {
    # Get valid coordinates
    valid_coords <- ensure_valid_coords(ints$node2_start[i], ints$node2_end[i], chr_length)
    
    # Only add if not already in the dataframe
    if (!ints$node2[i] %in% int_nodes$name) {
      int_nodes <- rbind(int_nodes, data.frame(
        name = ints$node2[i],
        chr = chr_name,
        start = valid_coords$start,
        end = valid_coords$end,
        genes = ints$node2_genes[i],
        type = "interaction",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create color scales for degree values
  # Get the range of degree values for both phenotypes to create a comparable scale
  max_degree_normal <- max(pov$degree_normal)
  max_degree_tnbc <- max(pov$degree_tnbc)
  max_degree <- max(max_degree_normal, max_degree_tnbc)
  
  degree_color_scale <- colorRamp2(
    c(0, max_degree/2, max_degree),
    c("#FFFFFF", "#FFD700", "#FF0000")  # White to gold to red
  )
  
  # Get the max zscore for consistent chord width scaling
  max_zscore <- max(ints$zscore)
  
  # Split interactions by phenotype
  normal_ints <- ints[ints$phenotype == "Normal", ]
  tnbc_ints <- ints[ints$phenotype == "TNBC", ]
  
  # Function to create a single diagram
  create_diagram <- function(phenotype, current_page) {
    # Select interactions for the current phenotype
    if (phenotype == "Normal") {
      current_ints <- normal_ints
      degree_column <- "degree_normal"
    } else {
      current_ints <- tnbc_ints
      degree_column <- "degree_tnbc"
    }
    
    # Set up the circlize environment
    circos.clear()
    
    # Set global parameters
    circos.par(
      gap.degree = 2,
      cell.padding = c(0, 0, 0, 0),
      start.degree = 90  # Start at midnight (12 o'clock)
    )
    
    # Initialize with ideogram
    circos.initializeWithIdeogram(
      chromosome.index = chr_name,
      plotType = c("ideogram", "labels")
    )
    
    # Add track for highlighting nodes
    circos.track(
      ylim = c(0, 1),
      panel.fun = function(x, y) {
        # Draw POV node
        pov_rect_col <- degree_color_scale(pov_node[[degree_column]])
        circos.rect(
          pov_node$start, 0,
          pov_node$end, 1,
          col = pov_rect_col,
          border = "black"
        )
        
        # Label POV gene with an arrow pointing to it
        text_pos <- (pov_node$start + pov_node$end) / 2
        circos.text(
          text_pos, 0.7,
          paste0("POV: ", pov_node$genes), 
          cex = 0.9, 
          facing = "inside", 
          niceFacing = TRUE,
          col = "black",
          font = 2
        )
        
        # Draw interaction nodes
        for (i in 1:nrow(int_nodes)) {
          node_in_phenotype <- any(current_ints$node2 == int_nodes$name[i])
          rect_col <- ifelse(node_in_phenotype, "#3CB371", "#AAAAAA")  # Green if present, grey if not
          
          circos.rect(
            int_nodes$start[i], 0,
            int_nodes$end[i], 1,
            col = rect_col,
            border = "black"
          )
          
          # Check if text fits within the rectangle
          width <- int_nodes$end[i] - int_nodes$start[i]
          node_center <- (int_nodes$start[i] + int_nodes$end[i]) / 2
          
          # Extract gene names (some might have multiple genes separated by semicolons)
          genes <- int_nodes$genes[i]
          
          # If there are multiple genes, try to abbreviate
          if (grepl(";", genes)) {
            # Split genes
            gene_list <- strsplit(genes, ";")[[1]]
            # If more than 2 genes, keep first and indicate there are more
            if (length(gene_list) > 2) {
              genes <- paste0(gene_list[1], "+", length(gene_list)-1)
            }
          }
          
          # Only label nodes in current phenotype
          if (node_in_phenotype) {
            # Try to place text at an appropriate position and size
            text_size <- ifelse(width > chr_length / 60, 0.7, 0.6)
            
            tryCatch({
              circos.text(
                node_center, 0.5,
                genes, 
                cex = text_size, 
                facing = "inside", 
                niceFacing = TRUE,
                col = "black"
              )
            }, error = function(e) {
              # If text placement fails, try a simpler approach or skip
              message(paste("Note: Could not place text for gene:", genes))
            })
          }
        }
      },
      track.height = 0.15,
      bg.border = NA
    )
    
    # Draw links (chords) from POV to interaction nodes
    pov_center <- (pov_node$start + pov_node$end) / 2
    
    for (i in 1:nrow(current_ints)) {
      # Get interaction node details
      int_node_idx <- which(int_nodes$name == current_ints$node2[i])
      
      if (length(int_node_idx) > 0) {
        # Get valid coordinates
        int_start <- int_nodes$start[int_node_idx]
        int_end <- int_nodes$end[int_node_idx]
        int_center <- (int_start + int_end) / 2
        
        # Calculate link width based on zscore and with a minimum width
        width_factor <- current_ints$zscore[i] / max_zscore
        link_width <- max(0.005, min(0.02, width_factor * 0.02))  # Scale between 0.005 and 0.02
        
        # Use darker colors for higher zscores
        alpha_level <- 0.4 + (width_factor * 0.4)  # Scale between 0.4 and 0.8
        link_color <- alpha("#0000FF", alpha_level)
        
        # Draw the link with error handling
        tryCatch({
          circos.link(
            chr_name, pov_center,
            chr_name, int_center,
            h.ratio = 0.5,
            w = link_width, 
            col = link_color
          )
        }, error = function(e) {
          message(paste("Note: Could not draw link to node:", int_nodes$name[int_node_idx]))
        })
      }
    }
    
    # Add a title
    title <- paste("Chromosome", sub("chr", "", chr_name), "Interactions -", phenotype, "Phenotype")
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    
    # Add a legend for the degree color scale
    create_legend()
  }
  
  # Function to create legends
  create_legend <- function() {
    # Create color legend for node degree
    pushViewport(viewport(x = 0.85, y = 0.2, width = 0.15, height = 0.3))
    
    # Add gradient legend
    n_steps <- 100
    for (i in 1:n_steps) {
      y_pos <- (i - 0.5) / n_steps
      degree_val <- (i - 1) / (n_steps - 1) * max_degree
      col <- degree_color_scale(degree_val)
      
      grid.rect(x = 0.5, y = y_pos, width = 0.8, height = 1/n_steps, 
                gp = gpar(fill = col, col = NA))
    }
    
    # Add labels
    grid.text("Node Degree", x = 0.5, y = 1.05, gp = gpar(fontsize = 9, fontface = "bold"))
    grid.text(paste0("High (", round(max_degree), ")"), x = 0.5, y = 0.95, gp = gpar(fontsize = 8))
    grid.text("Low (0)", x = 0.5, y = 0.05, gp = gpar(fontsize = 8))
    
    popViewport()
    
    # Create legend for node types
    pushViewport(viewport(x = 0.85, y = 0.55, width = 0.15, height = 0.15))
    
    # POV node
    grid.rect(x = 0.2, y = 0.8, width = 0.2, height = 0.1, 
              gp = gpar(fill = "#FF0000", col = "black"))
    grid.text("POV node", x = 0.6, y = 0.8, just = "left", gp = gpar(fontsize = 8))
    
    # Interaction node
    grid.rect(x = 0.2, y = 0.6, width = 0.2, height = 0.1, 
              gp = gpar(fill = "#3CB371", col = "black"))
    grid.text("Interaction node", x = 0.6, y = 0.6, just = "left", gp = gpar(fontsize = 8))
    
    # Inactive node
    grid.rect(x = 0.2, y = 0.4, width = 0.2, height = 0.1, 
              gp = gpar(fill = "#AAAAAA", col = "black"))
    grid.text("Inactive in phenotype", x = 0.6, y = 0.4, just = "left", gp = gpar(fontsize = 8))
    
    # Link width
    grid.text("Link width ∝ zscore", x = 0.5, y = 0.2, gp = gpar(fontsize = 8))
    
    popViewport()
  }
  
  # If output file is provided, create a PDF with both diagrams
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 16)
    
    # Create layout for two diagrams
    layout(matrix(1:2, nrow = 2))
    
    # Normal phenotype
    create_diagram("Normal", 1)
    
    # TNBC phenotype
    create_diagram("TNBC", 2)
    
    dev.off()
    
    return(paste("Plots saved to", output_file))
  } else {
    # Create interactive plots (one after another)
    # Normal phenotype
    create_diagram("Normal", 1)
    message("Press [Enter] for TNBC phenotype diagram...")
    invisible(readline())
    
    # TNBC phenotype
    create_diagram("TNBC", 2)
    
    return("Plots displayed in interactive mode")
  }
}


create_chromosome_diagrams(pten_list$pov, pten_list$ints,
                           "pten_interactions.pdf")


































# PXDC1

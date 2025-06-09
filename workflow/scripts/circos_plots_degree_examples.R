library(igraph)
library(circlize)

# function to get a gene's pov interactions
find_gene_connections <- function(graph_list, chromosome, gene_name) {
  # Get the graph for the specified chromosome
  graph <- graph_list[[chromosome]]
  
  if (is.null(graph)) {
    stop(paste("Chromosome", chromosome, "not found in the graph list"))
  }
  
  # Find vertices containing the gene name
  vertex_indices <- which(sapply(V(graph)$genes, function(genes) {
    if (is.null(genes) || is.na(genes) || genes == "") {
      return(FALSE)
    }
    gene_list <- strsplit(genes, ";")[[1]]
    return(gene_name %in% gene_list)
  }))
  
  if (length(vertex_indices) == 0) {
    message(paste("Gene", gene_name, "not found in chromosome", chromosome))
    return(NULL)
  }
  
  result <- list()
  for (i in seq_along(vertex_indices)) {
    vertex_index <- vertex_indices[i]
    vertex <- V(graph)[vertex_index]
    vertex_name <- V(graph)[vertex_index]$name
    
    # Get all neighbors connected to this vertex
    neighbors <- neighbors(graph, vertex)
    
    # Get all edges connected to this vertex
    edges <- incident(graph, vertex, mode = "all")
    
    # Create a subgraph with the vertex and its neighbors
    subgraph <- induced_subgraph(graph, c(vertex, neighbors))
    
    # Compile information about the vertex
    vertex_info <- vertex_attr(graph, index = vertex)
    
    # Compile information about connections
    edge_ends <- ends(graph, edges)
    connections <- data.frame(
      from = edge_ends[,1],
      to = edge_ends[,2],
      stringsAsFactors = FALSE
    )
    
    # Add connected vertex information and attributes
    connections$connected_vertex <- ifelse(
      connections$from == vertex_name, 
      connections$to, 
      connections$from
    )
    
    # Add connected vertex attributes (including chromosome, positions, and gene types)
    connections$connected_chr <- sapply(connections$connected_vertex, function(v_name) {
      v_index <- which(V(graph)$name == v_name)
      if (length(v_index) > 0 && v_name != vertex_name) {
        return(V(graph)[v_index]$chr)
      } else {
        return(NA)
      }
    })
    
    connections$connected_start <- sapply(connections$connected_vertex, function(v_name) {
      v_index <- which(V(graph)$name == v_name)
      if (length(v_index) > 0 && v_name != vertex_name) {
        return(V(graph)[v_index]$start)
      } else {
        return(NA)
      }
    })
    
    connections$connected_end <- sapply(connections$connected_vertex, function(v_name) {
      v_index <- which(V(graph)$name == v_name)
      if (length(v_index) > 0 && v_name != vertex_name) {
        return(V(graph)[v_index]$end)
      } else {
        return(NA)
      }
    })
    
    connections$connected_genes <- sapply(connections$connected_vertex, function(v_name) {
      v_index <- which(V(graph)$name == v_name)
      if (length(v_index) > 0 && v_name != vertex_name) {
        return(V(graph)[v_index]$genes)
      } else {
        return(NA)
      }
    })
    
    connections$connected_gene_types <- sapply(connections$connected_vertex, function(v_name) {
      v_index <- which(V(graph)$name == v_name)
      if (length(v_index) > 0 && v_name != vertex_name) {
        return(V(graph)[v_index]$gene_types)
      } else {
        return(NA)
      }
    })
    
    # Add current vertex chromosome and positions
    connections$vertex_chr <- chromosome
    connections$vertex_start <- vertex_info$start
    connections$vertex_end <- vertex_info$end
    connections$vertex_gene_types <- vertex_info$gene_types
    
    # Add edge attributes to connections dataframe
    for (attr_name in edge_attr_names(graph)) {
      connections[[attr_name]] <- edge_attr(graph, attr_name, edges)
    }
    
    result[[i]] <- list(
      vertex_info = vertex_info,
      connections = connections,
      subgraph = subgraph
    )
  }
  
  if (length(result) == 1) {
    return(result[[1]])
  } else {
    message(paste("Found", length(result), "vertices containing gene", gene_name))
    return(result)
  }
}


plot_interactions_pov <- function(pov_normal_list, pov_tnbc_list, topN=100, pov, tcga_annotation, show_non_coding=FALSE) {
  
  # stop if these things are not identical 
  stopifnot(identical(pov_normal_list$vertex_info, pov_tnbc_list$vertex_info))
  
  # STEP 1 build from_bed with pov information
  from_bed = data.frame(chr=pov_normal_list$vertex_info$chr,
                        start=pov_normal_list$vertex_info$start,
                        end=pov_normal_list$vertex_info$end)
  
  # STEP 2 EXTRACT THE TOP INTERACTIONS BY ZSCORE
  top_normal <- pov_normal_list$connections[
    order(pov_normal_list$connections$zscore, decreasing = T)[1:topN], 
  ]
  
  top_tnbc <- pov_tnbc_list$connections[
    order(pov_tnbc_list$connections$zscore, decreasing = T)[1:topN], 
  ]
  
  # STEP 3 calculate colors based on the zscore value
  all_zscores <- c(top_normal$zscore, top_tnbc$zscore)
  
  # Find global min and max
  min_zscore <- min(all_zscores)
  max_zscore <- max(all_zscores)
  
  library(viridis)
  
  # Create a blue gradient function
  get_blue_color <- function(zscore) {
    normalized <- (zscore - min_zscore) / (max_zscore - min_zscore)
    viridis::plasma(100)[ceiling(normalized * 99) + 1]
  }
  
  top_normal$color <- sapply(top_normal$zscore, get_blue_color)
  top_tnbc$color <- sapply(top_tnbc$zscore, get_blue_color)
  
  # Create a legend dataframe
  legend_breaks <- seq(min_zscore, max_zscore, length.out = 6)
  legend_colors <- sapply(legend_breaks, get_blue_color)
  
  legend_df <- data.frame(
    zscore = round(legend_breaks, 2),
    color = legend_colors
  )
  
  # STEP 4 extract to_bed objects for circos plotting
  to_bed_normal <- top_normal[,c("connected_chr", "connected_start", "connected_end", "color")]
  colnames(to_bed_normal) <- c("chr", "start", "end", "col")
  
  to_bed_tnbc <- top_tnbc[,c("connected_chr", "connected_start", "connected_end", "color")]
  colnames(to_bed_tnbc) <- c("chr", "start", "end", "col")
  
  # Convert tcga_annotation to data.frame for easier lookup
  tcga_df <- as.data.frame(tcga_annotation)
  
  # Helper function to get gene color based on type
  get_gene_color <- function(gene_name, pov_gene) {
    if (gene_name == pov_gene) {
      return("blue")
    }
    
    gene_info <- tcga_df[tcga_df$gene_name == gene_name, ]
    
    if (nrow(gene_info) > 0) {
      gene_type <- gene_info$gene_type[1]
      if (gene_type == "protein_coding") {
        return("black")
      } else {
        return("grey")
      }
    } else {
      return("grey")
    }
  }
  
  # Helper function to process genes
  process_genes <- function(connections_df, pov_info, pov_gene) {
    labels_bed <- data.frame()
    
    connections_with_genes <- connections_df[!is.na(connections_df$connected_genes), ]
    
    if (nrow(connections_with_genes) > 0) {
      for (i in 1:nrow(connections_with_genes)) {
        genes <- unlist(strsplit(connections_with_genes$connected_genes[i], ";"))
        
        for (gene in genes) {
          gene_color <- get_gene_color(gene, pov_gene)
          
          labels_bed <- rbind(labels_bed, data.frame(
            chr = connections_with_genes$connected_chr[i],
            start = connections_with_genes$connected_start[i],
            end = connections_with_genes$connected_end[i],
            value1 = gene,
            color = gene_color,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    pov_genes <- unlist(strsplit(pov_info$genes, ";"))
    
    for (gene in pov_genes) {
      gene_color <- get_gene_color(gene, pov_gene)
      
      labels_bed <- rbind(labels_bed, data.frame(
        chr = pov_info$chr,
        start = pov_info$start,
        end = pov_info$end,
        value1 = gene,
        color = gene_color,
        stringsAsFactors = FALSE
      ))
    }
    
    return(labels_bed)
  }
  
  # Create label beds for both conditions
  labels_bed_normal <- process_genes(top_normal, pov_normal_list$vertex_info, pov)
  labels_bed_tnbc <- process_genes(top_tnbc, pov_tnbc_list$vertex_info, pov)
  
  # Filter out grey genes if show_non_coding is FALSE
  if (!show_non_coding) {
    labels_bed_normal <- labels_bed_normal[labels_bed_normal$color != "grey", ]
    labels_bed_tnbc <- labels_bed_tnbc[labels_bed_tnbc$color != "grey", ]
  }
  
  collapse_genes <- function(labels_bed) {
    labels_bed$coord_id <- paste(labels_bed$chr, labels_bed$start, labels_bed$end, sep="_")
    
    collapsed <- aggregate(value1 ~ coord_id + chr + start + end, 
                           data = labels_bed, 
                           FUN = function(x) paste(x, collapse = ";"))
    
    collapsed$color <- sapply(1:nrow(collapsed), function(i) {
      coord <- collapsed$coord_id[i]
      colors <- labels_bed$color[labels_bed$coord_id == coord]
      if ("blue" %in% colors) {
        return("blue")
      } else if ("black" %in% colors) {
        return("black")
      } else {
        return("grey")
      }
    })
    
    collapsed$coord_id <- NULL
    
    return(collapsed)
  }
  
  # Collapse both label beds
  labels_bed_normal_collapsed <- collapse_genes(labels_bed_normal)
  labels_bed_tnbc_collapsed <- collapse_genes(labels_bed_tnbc)
  
  # Function to create a single circos plot
  create_single_circos <- function(labels_bed, to_bed, from_bed, topN, title_text) {
    function() {
      circos.clear()
      circos.par("start.degree" = 90)
      
      circos.initializeWithIdeogram(plotType = NULL, species = "hg38", chromosome.index = from_bed$chr)
      
      circos.genomicLabels(labels_bed, labels.column = 4, side = "outside",
                           cex = 0.6,
                           padding =  0.001, connection_height = 0.2, line_lwd = 0.45,
                           labels_height = min(c(convert_height(0.5, "cm"))),
                           col = labels_bed$color, niceFacing = TRUE, line_col = "darkgrey")
      
      circos.genomicIdeogram(species = "hg38")
      
      lapply(seq(1:topN), function(i){
        circos.genomicLink(from_bed, to_bed[i,], col = to_bed[i,"col"], 
                           border = NA)
      })
      
      circos.genomicAxis()
      text(0, 0, title_text, cex = 1.5, font = 2)
    }
  }
  
  # Create plot functions
  normal_plot <- create_single_circos(labels_bed_normal_collapsed, to_bed_normal, from_bed, topN, "Normal")
  tnbc_plot <- create_single_circos(labels_bed_tnbc_collapsed, to_bed_tnbc, from_bed, topN, "TNBC")
  
  # Create legend plot function
  legend_plot <- function() {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    
    # Create gradient legend
    n_colors <- nrow(legend_df)
    rect_height <- 0.8 / n_colors
    
    for (i in 1:n_colors) {
      rect(0.3, 0.1 + (i-1) * rect_height, 
           0.5, 0.1 + i * rect_height, 
           col = legend_df$color[i], 
           border = NA)
    }
    
    # Add z-score labels
    text(0.55, seq(0.1, 0.9, length.out = n_colors), 
         labels = sprintf("%.2f", legend_df$zscore),
         adj = 0, cex = 0.8)
    
    text(0.4, 0.95, "Z-score", cex = 1.2, font = 2)
    
    # Add gene color legend
    text(0.4, 0.02, "Gene colors:", cex = 0.8, font = 2)
    points(0.2, -0.05, pch = 19, col = "blue", cex = 1.2)
    text(0.25, -0.05, "POV gene", adj = 0, cex = 0.7)
    points(0.45, -0.05, pch = 19, col = "black", cex = 1.2)
    text(0.5, -0.05, "Protein coding", adj = 0, cex = 0.7)
    if (show_non_coding) {
      points(0.7, -0.05, pch = 19, col = "grey", cex = 1.2)
      text(0.75, -0.05, "Non-coding", adj = 0, cex = 0.7)
    }
  }
  
  # Return everything
  return(list(
    normal_plot = normal_plot,
    tnbc_plot = tnbc_plot,
    legend_plot = legend_plot,
    legend_df = legend_df,
    data = list(
      labels_normal = labels_bed_normal_collapsed,
      labels_tnbc = labels_bed_tnbc_collapsed,
      from_bed = from_bed,
      to_bed_normal = to_bed_normal,
      to_bed_tnbc = to_bed_tnbc
    )
  ))
}



#### diagram function 
# 
# pov_normal_list = pten_normal
# pov_tnbc_list = pten_tnbc
# pov="PTEN"
# 
# plot_interactions_pov <- function(pov_normal_list, pov_tnbc_list, topN=100, pov, tcga_annotation, show_non_coding=FALSE) {
#   
#   # stop if these things are not identical 
#   stopifnot(identical(pov_normal_list$vertex_info, pov_tnbc_list$vertex_info))
#   
#   # STEP 1 build from_bed with pov information, it only has three fields, get this from the vertex_info
#   from_bed = data.frame(chr=pov_normal_list$vertex_info$chr,
#                         start=pov_normal_list$vertex_info$start,
#                         end=pov_normal_list$vertex_info$end)
#   
#   # STEP 2 EXTRACT THE TOP 100 INTERACTIONS BY ZSCORE FROM EACH "to" DATAFRAME
#   # so 
#   # to_bed_normal is
#   top_normal <- pov_normal_list$connections[
#     order(pov_normal_list$connections$zscore, decreasing = T)[1:topN], 
#     ]
#   
#   # to_bed_tnbc is
#   top_tnbc <- pov_tnbc_list$connections[
#     order(pov_tnbc_list$connections$zscore, decreasing = T)[1:topN], 
#     ]
#   
#   # STEP 3 calculate colors based on the zscore value and add as a column to each to_bed dataframe
#   all_zscores <- c(top_normal$zscore, top_tnbc$zscore)
#   
#   # Find global min and max
#   min_zscore <- min(all_zscores)
#   max_zscore <- max(all_zscores)
#   
#   # Using the viridis package for a perceptually uniform color scale
#   # This is colorblind-friendly and prints well in grayscale
#   library(viridis)
#   
#   # Create a blue gradient function
#   get_blue_color <- function(zscore) {
#     # Normalize to 0-1 range
#     normalized <- (zscore - min_zscore) / (max_zscore - min_zscore)
#     
#     # Get color from blue palette
#     # Using 'plasma' scale for high visual distinction between values
#     viridis::plasma(100)[ceiling(normalized * 99) + 1]
#   }
#   
#   top_normal$color <- sapply(top_normal$zscore, get_blue_color)
#   top_tnbc$color <- sapply(top_tnbc$zscore, get_blue_color)
#   
#   # Create a legend dataframe with regular intervals for reference
#   legend_breaks <- seq(min_zscore, max_zscore, length.out = 6)
#   legend_colors <- sapply(legend_breaks, get_blue_color)
#   
#   #### OUT Print for reference
#   legend_df <- data.frame(
#     zscore = round(legend_breaks, 2),
#     color = legend_colors
#   )
#   #print(legend_df)
#   
#   # STEP 4 extract to_bed objects for circos plotting
#   to_bed_normal <- top_normal[,c("connected_chr", "connected_start", "connected_end", "color")]
#   colnames(to_bed_normal) <- c("chr", "start", "end", "col")
#   
#   to_bed_tnbc <- top_tnbc[,c("connected_chr", "connected_start", "connected_end", "color")]
#   colnames(to_bed_tnbc) <- c("chr", "start", "end", "col")
#   
#   # STEP 4 LABELS TAKE EACH top DATAFRAME AND EXTRACT THE GENE NAMES AND COORDINATES
#   # Helper function to process genes
#   # Convert tcga_annotation to data.frame for easier lookup
#   tcga_df <- as.data.frame(tcga_annotation)
#   
#   # Helper function to get gene color based on type
#   get_gene_color <- function(gene_name, pov_gene) {
#     if (gene_name == pov_gene) {
#       return("blue")
#     }
#     
#     # Look up gene type in tcga_annotation
#     gene_info <- tcga_df[tcga_df$gene_name == gene_name, ]
#     
#     if (nrow(gene_info) > 0) {
#       gene_type <- gene_info$gene_type[1]
#       if (gene_type == "protein_coding") {
#         return("black")
#       } else {
#         return("grey")
#       }
#     } else {
#       # Default to grey if gene not found
#       return("grey")
#     }
#   }
#   
#   # Helper function to process genes
#   process_genes <- function(connections_df, pov_info, pov_gene) {
#     labels_bed <- data.frame()
#     
#     # Process connected genes (remove NA rows first)
#     connections_with_genes <- connections_df[!is.na(connections_df$connected_genes), ]
#     
#     # Process each row with genes
#     if (nrow(connections_with_genes) > 0) {
#       for (i in 1:nrow(connections_with_genes)) {
#         genes <- unlist(strsplit(connections_with_genes$connected_genes[i], ";"))
#         
#         for (gene in genes) {
#           gene_color <- get_gene_color(gene, pov_gene)
#           
#           labels_bed <- rbind(labels_bed, data.frame(
#             chr = connections_with_genes$connected_chr[i],
#             start = connections_with_genes$connected_start[i],
#             end = connections_with_genes$connected_end[i],
#             value1 = gene,
#             color = gene_color,
#             stringsAsFactors = FALSE
#           ))
#         }
#       }
#     }
#     
#     # Add POV genes
#     pov_genes <- unlist(strsplit(pov_info$genes, ";"))
#     
#     for (gene in pov_genes) {
#       gene_color <- get_gene_color(gene, pov_gene)
#       
#       labels_bed <- rbind(labels_bed, data.frame(
#         chr = pov_info$chr,
#         start = pov_info$start,
#         end = pov_info$end,
#         value1 = gene,
#         color = gene_color,
#         stringsAsFactors = FALSE
#       ))
#     }
#     
#     return(labels_bed)
#   }
#   
#   # Create label beds for both conditions
#   labels_bed_normal <- process_genes(top_normal, pov_normal_list$vertex_info, pov)
#   labels_bed_tnbc <- process_genes(top_tnbc, pov_tnbc_list$vertex_info, pov)
#   
#   # Filter out grey genes if show_non_coding is FALSE
#   if (!show_non_coding) {
#     labels_bed_normal <- labels_bed_normal[labels_bed_normal$color != "grey", ]
#     labels_bed_tnbc <- labels_bed_tnbc[labels_bed_tnbc$color != "grey", ]
#   }
#   
#   collapse_genes <- function(labels_bed) {
#     # Create a unique identifier for each coordinate
#     labels_bed$coord_id <- paste(labels_bed$chr, labels_bed$start, labels_bed$end, sep="_")
#     
#     # Group by coordinates and collapse
#     collapsed <- aggregate(value1 ~ coord_id + chr + start + end, 
#                            data = labels_bed, 
#                            FUN = function(x) paste(x, collapse = ";"))
#     
#     # Determine color for collapsed entries
#     # If any gene is blue (POV), the whole group is blue
#     # Otherwise, if any is black (protein_coding), the group is black
#     # Otherwise, grey
#     collapsed$color <- sapply(1:nrow(collapsed), function(i) {
#       coord <- collapsed$coord_id[i]
#       colors <- labels_bed$color[labels_bed$coord_id == coord]
#       if ("blue" %in% colors) {
#         return("blue")
#       } else if ("black" %in% colors) {
#         return("black")
#       } else {
#         return("grey")
#       }
#     })
#     
#     # Remove the coord_id column
#     collapsed$coord_id <- NULL
#     
#     return(collapsed)
#   }
#   
#   # Collapse both label beds
#   labels_bed_normal_collapsed <- collapse_genes(labels_bed_normal)
#   labels_bed_tnbc_collapsed <- collapse_genes(labels_bed_tnbc)
#   
#   ### NORMAL PLOT
#   circos.par("start.degree" = 90)
#   
#   circos.initializeWithIdeogram(plotType = NULL, species = "hg38", chromosome.index = "chr10")
#   
#   
#   # this takes in a bed to label the genes
#   circos.genomicLabels(labels_bed_normal_collapsed, labels.column = 4, side = "outside",
#                        cex = 0.65,
#                        padding =  0.001, connection_height = 0.2,line_lwd = 0.45,
#                        labels_height = min(c(convert_height(0.5, "cm"))),
#                        col = labels_bed_normal_collapsed$color, niceFacing = TRUE, line_col = "darkgrey")
#   
#   
#   circos.genomicIdeogram(species = "hg38")
#   
#   lapply(seq(1:topN), function(i){
#     circos.genomicLink(from_bed, to_bed_normal[i,], col = to_bed_normal[i,"col"], 
#                        border = NA)
#   })
#   
#   
#   circos.genomicAxis()
#   
#   
#   ### TNBC PLOT
#   #circos.clear()
#   circos.par("start.degree" = 90)
#   
#   circos.initializeWithIdeogram(plotType = NULL, species = "hg38", chromosome.index = "chr10")
#   
#   
#   # this takes in a bed to label the genes
#   circos.genomicLabels(labels_bed_tnbc_collapsed, labels.column = 4, side = "outside",
#                        cex = 0.65,
#                        padding =  0.001, connection_height = 0.2,line_lwd = 0.45,
#                        labels_height = min(c(convert_height(0.5, "cm"))),
#                        col = labels_bed_tnbc_collapsed$color, niceFacing = TRUE, line_col = "darkgrey")
#   
#   circos.genomicIdeogram(species = "hg38")
#   
#   lapply(seq(1:topN), function(i){
#     circos.genomicLink(from_bed, to_bed_tnbc[i,], col = to_bed_tnbc[i,"col"], 
#                        border = NA)
#   })
#   
#   
#   circos.genomicAxis()
#   
#   
#   
# }





# set.seed(123)
# 
# bed1 = generateRandomBed(nr = 100)
# bed1 = bed1[sample(nrow(bed1), 20), ]
# #bed1 <- head(bed1)
# 
# bed2 = generateRandomBed(nr = 100)
# bed2 = bed2[sample(nrow(bed2), 20), ]
# #bed2 <- head(bed2)
# 
# bed1$chr <- "chr1"
# bed2$chr <- "chr1"
# 
# bed1$value1 <- NULL
# 
# bed1$end <- bed1$start+100000
#
# circos.par("start.degree" = 90)
# 
# circos.initializeWithIdeogram(species = "hg38", chromosome = "chr1")
# 
# circos.genomicLink(bed1[1,], bed2[5,], col = rand_color(nrow(bed1[1,]), transparency = 0.5), 
#                    border = NA)
# 
# # somehow we need to take it outside, like track wise
# circos.labels(sectors = bed1$chr, x = bed1$start, 
#               labels = paste("gene", seq(1:20), sep="_"),
#               side = "outside")
#
## second attempt
# bed<-data.frame(chr="chr1", start=24270817, end=24417739, value1="myLittleGene")
# 
# bed$start <- bed1$start[1]
# bed$end <- bed1$end[1]
# 
# bed_lab <- bed
# bed_lab$start <- bed2$start[6]
# bed_lab$end <- bed2$end[6]
# bed_lab$value1 <- "targetG"
# 
# bed <- rbind(bed, bed)
# bed$value1[2] <- "GeneBig"
#
# should thickness be related to the hi-c interaction strength?
# size of the chromosome divided by something times the hic count??
# labels should always have the middle point
# 
# 
# circos.clear()
# circos.par("start.degree" = 90)
# 
# circos.initializeWithIdeogram(plotType = NULL, species = "hg38", chromosome.index = "chr1")
# 
# 
# # this takes in a bed to label the genes
# # we can indicate a color vector in the bed for lncRNa and coding genes
# circos.genomicLabels(rbind(bed, bed_lab), labels.column = 4, side = "outside",
#                      cex = 0.75,
#                      padding =  0.001, connection_height = 0.05,
#                      labels_height = min(c(convert_height(0.5, "cm"))),
#                      col = c("black", "red"), niceFacing = TRUE)
# 
# circos.genomicIdeogram(species = "hg38")
# 
# # there is an lapply here for the second bed (the pov stays the same)
# circos.genomicLink(bed1[1,], bed2[6,], col = rand_color(nrow(bed1[1,]), transparency = 0.5), 
#                    border = NA)
# 
# circos.genomicLink(bed1[1,], bed2[7,], col = rand_color(nrow(bed1[1,]), transparency = 0.5), 
#                    border = NA)
# 
# circos.genomicLink(bed1[1,], bed2[13,], col = rand_color(nrow(bed1[1,]), transparency = 0.5), 
#                    border = NA)
# 
# circos.genomicLink(bed1[1,], bed2[12,], col = bed2[12,"col"], 
#                    border = NA)
# 
# circos.genomicAxis()

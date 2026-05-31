#' @title Convert Taxonomy to Phylogenetic Tree
#'
#' @description Converts a taxonomy data frame (with columns from root level to species level) into an 'ape' \code{phylo} object.
#'
#' @param df A data frame where rows are species and columns represent taxonomic levels ordered from root (left) to tip (right). The last column contains the tip/species labels.
#'
#' @return An object of class \code{phylo}.
#'
#' @export
taxToPhylo <- function(df) {
  # Reorder: Root (Left) to Tips (Right), df should have columns: Level1, Level2, ..., Tips
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  for(i in 1:n_cols) {
    df[[i]] <- as.character(df[[i]])
  }

  # 1. Get unique node names for all internal levels
  # We exclude the last column (tips) from internal nodes
  internal_levels <- df[, -n_cols, drop = FALSE]
  
  # Create unique identifiers for every node at every level
  # This prevents "Genus_A" in "Family_X" being confused with "Genus_A" in "Family_Y"
  for(i in 2:ncol(internal_levels)) {
    internal_levels[, i] <- paste(internal_levels[, i-1], internal_levels[, i], sep = "_")
  }
  
  # Get all unique internal nodes
  all_internal_nodes <- c("root", unique(as.vector(as.matrix(internal_levels))))
  num_internal <- length(all_internal_nodes)
  
  # 2. Assign IDs
  # Tips: 1 to n_rows
  # Internal Nodes: (n_rows + 1) to (n_rows + num_internal)
  node_labels <- c(df[[n_cols]], all_internal_nodes)
  node_ids <- seq_along(node_labels)
  names(node_ids) <- node_labels
  
  # 3. Build Edges
  edges <- matrix(NA, ncol = 2, nrow = 0)
  
  # Edge: Parent (Column i) -> Child (Column i+1)
  for(i in 0:(n_cols - 1)) {
    if(i == 0) {
      parent_names <- "root"
    } else {
      parent_names <- internal_levels[, i]
    }
    
    if(i == (n_cols - 1)) {
      child_names <- df[[n_cols]] # The actual OTU names
    } else {
      child_names <- internal_levels[, i+1]
    }
    
    # Get unique parent-child pairs at this level
    pairs <- unique(data.frame(p = parent_names, c = child_names))
    
    # Map names to IDs
    level_edges <- cbind(node_ids[pairs$p], node_ids[pairs$c])
    edges <- rbind(edges, level_edges)
  }
  
  # 4. Assemble 'phylo' object
  # Root must be the node that is never a child
  all_children <- edges[, 2]
  all_parents <- unique(edges[, 1])
  root_node <- all_parents[!(all_parents %in% all_children)]
  
  obj <- list(
    edge = matrix(as.integer(edges), ncol = 2),
    tip.label = df[[n_cols]],
    node.label = all_internal_nodes,
    Nnode = num_internal,
    edge.length = rep(1, nrow(edges))
  )
  class(obj) <- "phylo"
  return(obj)
}

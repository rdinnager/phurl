#' Prepare data for analysing with phurl models
#'
#' @param tree Phylogenetic tree of class \code{\link[ape]{phylo}} or \code{\link[tibble:tbl_df-class]{tibble}}
#' (e.g. using \code{\link[tidytree]{as_tibble}})
#' @param traits Trait data. Can be a named \code{vector}, a \code{matrix} or \code{data.frame} with named rows, where the rownames
#' are the species names, or a \code{data.frame} or \code{\link[tibble:tbl_df-class]{tibble}} with a column specifying the
#' species name. If species names are not referenced in either the rownames or a column, data will be matched
#' to the tree using order only (e.g. assuming data rows are in the same order as \code{tree$tip.label}).
#' @param species_column Character value specifying which column of \code{traits} refers to the species
#' names, if they are not specified in the rownames. If not \code{NULL}, this will override any rownames in
#' \code{traits}
#' @param scale_edges Should the edges be scaled to a maximum of 1? This is highly recommended for regularised
#' models. Default TRUE.
#' @param centre Should the resulting edge path predictor matrix be centred to so that each column has a mean
#' of 0? This is recommended to improve convergence of the models. Default TRUE.
#' @param standardise_traits Should traits values be standardised by subtracting the mean and dividing by
#' them standard deviation? Recommended. Default TRUE.
#'
#' @return A object of class \code{\link{phurl}}
#' @export
#'
#' @examples
phurl_prepare_data <- function(tree, traits, species_column = NULL, 
                               scale_edges = TRUE, centre = TRUE,
                               transform_traits = NULL,
                               standardise_traits = TRUE) {
  temp_tree <- tree
  max_depth <- max(ape::node.depth.edgelength(temp_tree))
  if(scale_edges) {
    temp_tree$edge.length <- temp_tree$edge.length / max_depth
  }
  temp_tree$root.edge <- 1
  tree_dat <- tidytree::as_tibble(temp_tree) %>%
    dplyr::filter(node <= ape::Ntip(temp_tree))

  edge_mat <- RRphylo::makeL(temp_tree)
  tip_df <- t(apply(edge_mat, 1, function(x) {x[x != 0][-1] <- rev(cumsum(rev(x[x != 0][-1]))); x})) %>%
    dplyr::as_tibble() %>%
    stats::setNames(paste0("node_", names(.))) %>%
    dplyr::mutate(species = rownames(edge_mat)) %>%
    dplyr::select(species, dplyr::everything())
  
  root_node <- list(root_value = dplyr::sym(paste0("node_", length(temp_tree$tip.label) + 1)))

  first_splits_edges <- temp_tree$edge[temp_tree$edge[ , 1] == (length(temp_tree$tip.label) + 1), 2]
  first_splits_name <- paste0("first_split_", first_splits_edges)
  first_splits = list()
  first_splits[[first_splits_name[1]]] <- dplyr::sym(paste0("node_", first_splits_edges[1]))
  first_splits[[first_splits_name[2]]] <- dplyr::sym(paste0("node_", first_splits_edges[2]))
  
  tip_df <- tip_df %>%
    dplyr::rename(!!!root_node) %>%
    dplyr::rename(!!!first_splits)
  
  
  node_mat <- RRphylo::makeL1(temp_tree)
  node_df <- t(apply(node_mat, 1, function(x) {x[x != 0][-1] <- rev(cumsum(rev(x[x != 0][-1]))); x})) %>%
    dplyr::as_tibble() %>%
    stats::setNames(paste0("node_", names(.))) %>%
    dplyr::mutate(node = paste0("node_", rownames(node_mat))) %>%
    dplyr::select(node, dplyr::everything())
  
  node_df <- node_df %>%
    dplyr::rename(!!!root_node) %>%
    dplyr::rename(!!!first_splits)
  
  if(centre) {
    tip_df <- tip_df %>%
      dplyr::mutate_at(dplyr::vars(-species, -root_value),
                       ~ (. - mean(., na.rm = TRUE)))
    
    node_df <- node_df %>%
      dplyr::mutate_at(dplyr::vars(-node, -root_value),
                       ~ (. - mean(., na.rm = TRUE)))
  }
  
  
  nodes_on_path_to_tips <- apply(edge_mat, 1, function(x) paste0("node_", names(which(x > 0))))
  names(nodes_on_path_to_tips) <- paste0("node_", names(nodes_on_path_to_tips))
  nodes_on_path_to_nodes <- apply(node_mat, 1, function(x) paste0("node_", names(which(x > 0))))
  names(nodes_on_path_to_nodes) <- paste0("node_", names(nodes_on_path_to_nodes))
  
  if(!is.null(species_column)){
    tree_df <- tree_dat %>%
      dplyr::left_join(traits, by = c("label" = species_column)) %>%
      dplyr::rename(species = label)
  } else {
    if(is.vector(traits)) {
      if(is.null(names(traits))) {
        trait_df <- dplyr::tibble(species = temp_tree$tip.label,
                                  trait = traits)
        warning("traits object has no names, assuming same order as tree$tip.label")
      } else {
        trait_df <- dplyr::tibble(species = names(traits),
                                  trait = traits)
      }
    } else {
       if(!is.null(rownames(traits))) {
        trait_df <- traits %>%
          as.data.frame() %>%
          tibble::rownames_to_column("species")
       } else {
         trait_df <- traits %>%
           as.data.frame() %>%
           dplyr::mutate(species = temp_tree$tip.label) %>%
           dplyr::select(species, dplyr::everything())
          warning("traits object has no rownames, assuming same order as tree$tip.label")
        }
    }
    tree_df <- tree_dat %>%
      dplyr::left_join(trait_df, by = c("label" = "species")) %>%
      dplyr::rename(species = label)
  }
  
  if(standardise_traits) {
    mean_sd <- tree_df %>%
      dplyr::summarise_at(dplyr::vars(-species, -node, -branch.length, -parent),
                          list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))
    tree_df <- tree_df %>%
      dplyr::mutate_at(dplyr::vars(-species, -node, -branch.length, -parent),
                       ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))
  }

  phurl_ob <- list(phylo = temp_tree, 
                   y_data = tree_df, 
                   x_data_tips = tip_df,
                   x_data_nodes = node_df,
                   nodes_on_path_to_tips = nodes_on_path_to_tips,
                   nodes_on_path_to_nodes = nodes_on_path_to_nodes,
                   model_fit = NULL, 
                   plots = NULL,
                   phylo_unmodified = tree,
                   edge_scaled = scale_edges,
                   edge_centred = centre,
                   y_standardisation = mean_sd)
  class(phurl_ob) <- "phurl"
  phurl_ob
}

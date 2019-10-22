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
#'
#' @return A object of class \code{\link{phurl}}
#' @export
#'
#' @examples
phurl_prepare_data <- function(tree, traits, species_column = NULL) {
  temp_tree <- tree
  max_depth <- max(ape::node.depth.edgelength(temp_tree))
  temp_tree$edge.length <- temp_tree$edge.length / max_depth
  temp_tree$root.edge <- 1
  tree_dat <- tidytree::as_tibble(temp_tree)

  edge_mat <- RRphylo::makeL(temp_tree)
  tip_df <- t(apply(edge_mat, 1, function(x) {x[x != 0][-1] <- rev(cumsum(rev(x[x != 0][-1]))); x})) %>%
    dplyr::as_tibble() %>%
    stats::setNames(paste0("edge_beginning_at_node_", names(.))) %>%
    dplyr::mutate(species = rownames(edge_mat)) %>%
    dplyr::select(species, dplyr::everything())
  
  root_node <- list(root_value = dplyr::sym(paste0("edge_beginning_at_node_", length(temp_tree$tip.label) + 1)))

  first_splits_edges <- temp_tree$edge[temp_tree$edge[ , 1] == (length(temp_tree$tip.label) + 1), 2]
  first_splits_name <- paste0("first_split_", first_splits_edges)
  first_splits = list()
  first_splits[[first_splits_name[1]]] <- dplyr::sym(paste0("edge_beginning_at_node_", first_splits_edges[1]))
  first_splits[[first_splits_name[2]]] <- dplyr::sym(paste0("edge_beginning_at_node_", first_splits_edges[2]))
  
  tip_df <- tip_df %>%
    dplyr::rename(!!!root_node) %>%
    dplyr::rename(!!!first_splits)
  
  
  node_mat <- RRphylo::makeL1(temp_tree)
  node_df <- t(apply(node_mat, 1, function(x) {x[x != 0][-1] <- rev(cumsum(rev(x[x != 0][-1]))); x})) %>%
    dplyr::as_tibble() %>%
    stats::setNames(paste0("edge_beginning_at_node_", names(.))) %>%
    dplyr::mutate(node = paste0("node_", rownames(node_mat))) %>%
    dplyr::select(node, dplyr::everything())
  
  node_df <- node_df %>%
    dplyr::rename(!!!root_node) %>%
    dplyr::rename(!!!first_splits)
  
  if(!is.null(species_column)){
    tree_df <- tree_dat %>%
      dplyr::left_join(traits, by = c("label" = species_column))
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
      dplyr::left_join(trait_df, by = c("label" = "species"))
  }

  phurl_ob <- list(phylo = temp_tree, 
                   y_data = tree_df, 
                   x_data_tips = tip_df,
                   x_data_nodes = node_df,
                   model_fit = NULL, 
                   plots = NULL,
                   phylo_unmodified = tree)
  class(phurl_ob) <- "phurl"
  phurl_ob
}

#' Simulate Continuous Traits on a Phylogeny
#'
#' @param x Object to be simulated from. Can be either a \code{phurl} object generated from
#' \code{\link{phurl_prepare_data}}, or an object of class \code{phylo}
#' @param family What error structure to use in the simulations? Can be a character vector of length \code{n_trait}
#' to use different families for different traits.
#' @param n_trait Number of traits to simulate.
#' @param root_value A vector of root values for each trait. Will be recycled to length \code{n_trait}
#' @param initial_rate Vector of standard deviations for initial rates. Rate at which traits are diverging after first split (will
#' be drawn from a normal distribution with \code{sd = initial_rate}). Will be recycled to length 
#' \code{n_trait}
#' @param tau A vector tau values. This parameter controls the degree to which trait evolution rates change at
#' phylogeny nodes. Rate changes are sampled from a gaussian distribution with mean zero and sd equal to tau,
#' scaled by the edge length of the parent edge. Will be recycled to length \code{n_trait}
#' @param sd Gaussian error of final tip-level observations. Standard deviation of noise added to tip-level
#' trait values. Only used if \code{family="gaussian"}
#'
#' @return
#' @export
#'
#' @examples
prl_sim_traits_continuous <- function(x, family = c("gaussian"), n_trait = 1, root_value = 0, initial_rate = 2, tau = 5, sd = 0.25) {
  
  if(!inherits(x, "phurl")) {
    x <- prl_prepare_data(tree = x, centre = FALSE)
  }
  
  design_tips <- x$x_data_tips %>%
    dplyr::select(-species, -root_value) %>%
    as.matrix
  
  design_nodes <- x$x_data_nodes %>%
    dplyr::select(-node, -root_value) %>%
    as.matrix() %>%
    cbind(matrix(0, ncol = ape::Ntip(x$phylo), nrow = nrow(x$x_data_nodes)))
  
  tree_dat <- tidytree::as_tibble(x$phylo)
  parent_branches <- purrr::map_dfr(tree_dat$node, ~tidytree::parent(tree_dat, .) %>%
                                  dplyr::select(branch.length) %>%
                                  dplyr::mutate(node = .x)) %>%
    dplyr::mutate(node = paste0("node_", node))
  
  first_split_replace <- c(names(x$first_splits), names(x$root_node))
  names(first_split_replace) <- paste0(c(as.character(unlist(x$first_splits)), as.character(x$root_node[[1]])), "$")
  
  parent_branches$node <- stringr::str_replace_all(parent_branches$node,
                                                   first_split_replace)
  
  edge_weights <- parent_branches$branch.length
  names(edge_weights) <- parent_branches$node
  
  generate_trait <- function() {
  
    coefs <- stats::rnorm(ncol(design_tips), 0, tau * edge_weights[colnames(design_tips)])
    
    first_splits <- which(startsWith(colnames(design_tips), "first_split"))
    coefs[first_splits] <- rnorm(2, sd = initial_rate)
    
    mu_tips <- root_value + design_tips %*% coefs
    mu_nodes <- root_value + design_nodes %*% coefs
    
    y_tips <- stats::rnorm(length(mu_tips), mu_tips, sd)
    y_nodes <- mu_nodes
    
    list(y_tips = y_tips, y_nodes = y_nodes)
  }
  
  ys <- replicate(n_trait, generate_trait(), simplify = FALSE)
  names(ys) <- paste0("trait_", seq_len(n_trait))
  
  tips_ys <- lapply(ys, function(z) z$y_tips) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(label = x$x_data_tips$species)
  
  nodes_ys <- lapply(ys, function(z) z$y_nodes) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(node = x$x_data_nodes$node)
  
  trait_df <- tidytree::as_tibble(x$phylo) %>%
    dplyr::filter(node <= ape::Ntip(x$phylo)) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(tips_ys, by = "label") %>%
    dplyr::rename(species = label) %>%
    dplyr::select(-parent, -node, -branch.length)
  
  trait_nodes_df <- dplyr::tibble(node = x$x_data_nodes$node) %>%
    dplyr::left_join(nodes_ys, by = "node") 
  
  res <- list(traits = trait_df, node_traits = trait_nodes_df)
  class(res) <- "phurl-sim"
  return(res)
}

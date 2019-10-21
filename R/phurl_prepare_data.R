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
  tree_dat <- tidytree::as_tibble(tree)

  edge_mat <- RRphylo::makeL(tree)
  new_mat <- t(apply(edge_mat, 1, function(x) {x[x != 0] <- rev(cumsum(rev(x[x != 0]))); x}))

  if(!is.null(species_column)){
    tree_df <- tree_dat %>%
      dplyr::left_join(traits, by = c("label" = species_column))
  } else {
    if(is.vector(traits)) {
      if(is.null(names(traits))) {
        trait_df <- dplyr::tibble(species = tree$tip.label,
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
           dplyr::mutate(species = tree$tip.label) %>%
           dplyr::select(species, dplyr::everything())
          warning("traits object has no rownames, assuming same order as tree$tip.label")
        }
    }
    tree_df <- tree_dat %>%
      dplyr::left_join(trait_df, by = c("label" = "species"))
  }

  phurl_ob <- list(phylo = tree, y_data = tree_df, model_fit = NULL, plots = NULL)
  class(phurl_ob) <- "phurl"
  phurl_ob
}

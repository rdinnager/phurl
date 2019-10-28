#' Title
#'
#' @param x Object to be simulated from. Can be either a \code{phurl} object generated from
#' \code{\link{phurl_prepare_data}}, or an object of class \code{phylo}
#' @param family What error structure to use in the simulations? Can be a character vector of length \code{n_trait}
#' to use different families for different traits.
#' @param n_trait Number of traits to simulate.
#' @param root_value A vector of root values for each trait. Will be recycled to length \code{n_trait}
#' @param initial_rate A vector of initial rates. Rate at which traits are diverging after first split (will
#' be assumed to be c(initial_rate, -initial_rate) along the first splitting edges). Will be recycled to length 
#' \code{n_trait}
#' @param tau A vector tau values. This parameter controls the degree to which trait evolution rates change at
#' phylogeny nodes. Rate changes are sampled from a gaussian distribution with mean zero and sd equal to tau.
#' Will be recycled to length \code{n_trait}
#' @param sd Gaussian error of final tip-level observations. Standard deviation of noise added to tip-level
#' trait values. Only used if \code{family="gaussian"}
#'
#' @return
#' @export
#'
#' @examples
phurl_sim_traits_continuous <- function(x, family = c("gaussian"), n_trait = 1, root_value = 0, initial_rate, tau = 1, sd = 0.25) {
  
  if(inherits(x, "phurl")) {
    design <- x$x_data_tips %>%
      dplyr::select(-species, -root_value) %>%
      as.matrix
  } else {
    design <- phurl_prepare_data(tree = x)$x_data_tips %>%
      dplyr::select(-species, -root_value) %>%
      as.matrix
  }
  
  coefs <- stats::rnorm(ncol(design), 0, tau)
  
  first_splits <- which(startsWith(colnames(design), "first_split"))
  coefs[first_splits] <- c(initial_rate, -initial_rate)
  
  mu <- root_value + design %*% coefs
  
  y <- stats::rnorm(length(mu), mu, sd)
}

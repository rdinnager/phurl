#' Title Fit a continuous model of trait evolution along a phylogeny, using a Bayesian method
#'
#' This function fits a Bayesian version of Ridge regression on edge paths, producing a phylogenetically
#' smoothed estimate of trait evolution rates along a phylogeny's branches. Shrinkage limits how much
#' change in rates can occur at each split of the phylogeny.
#'
#' @param x Object of class \code{phurl}, such as that produced by \code{\link{phurl_prepare_data}},
#' or a matrix or data.frame of traits values to be fit (with rownames or order matching tree$tip.label).
#' @param family a character vector specifying the statistical families to be used for each trait. If named,
#' names will be used to match the family to a column in the traits data, otherwise they will be used in order.
#' If the length of the vector is less than the number of traits, it will be recycled with a warning.
#' @param method Method to use. Currently the only option is "Ridge". This is mainly for compatibility
#' with \code{\link{fit_trait_shift}}, which allows multiple methods.
#' @param tree A phylogeny of class \code{phylo}. Only used if x is a matrix of trait values (if x is an object of
#' class \code{phurl} it should have a tree already associated with it)
#' @param regularise_first_split Should the edges beginning at the first split in the phylogeny be regularised?
#' If not, these are treated like a a set of 'intercepts' for evolutionary rates, and every subsequent change
#' in these rates along the tree are regularised (penalised).
#' @param n_samples Number of MCMC samples to draw during the fitting process.
#' @param thin Thinning rate for MCMC samples.
#' @param warmup Number of iteration to use for warmup of the MCMC sampling process.
#' @param chains Number of replicate MCMC chains to run.
#' @param n_cores Number of CPU cores to use for sampling.
#' @param ... Any other parameters to pass to the underlying fitting function (in this case \code{\link[greta]{mcmc}})
#'
#' @details Test equations:
#' \out{<pre>y<sub>i</sub> ~ f(μ, φ)
#'
#'              ===
#'              \
#'   g(μ) = α + /   β<sub>j</sub>
#'              ===
#'                j
#'
#' β<sub>j</sub> ~ Normal(0, τ)
#'
#' τ ~ InvGamma(1, 1)</pre>}
#'
#'
#'
#'
#' @return Object of class \code{phurl}, containing a fitted model in the \code{model} slot, and a
#' tibble in the \code{results} slot containing estimates of rates along each edge, estimates of
#' rate changes at each node, and ancestral character estimates (ACEs), including samples from the full
#' posterior distribution.
#' @export
#'
#' @examples
phurl_fit_traits_continuous_bayesian <- function(x, family = c("gaussian"), method = c("Ridge"),
                                                 tree = NULL, regularise_first_split = FALSE,
                                                 n_samples = 1000, thin = 1,
                                                 warmup = 1000, chains = 4, n_cores = NULL, ...) {
  if(!inherits(x, "phurl")) {
    x <- phurl_prepare_data(tree, x)
  }

  if(regularise_first_split) {
    design <- x$x_data_tips %>%
      dplyr::select(-species, -root_value, -dplyr::starts_with("first_split")) %>%
      as.matrix %>%
      greta::as_data()
  } else {
    design <- x$x_data_tips %>%
      dplyr::select(-species, -root_value) %>%
      as.matrix %>%
      greta::as_data()
  }

  y <- x$y_data %>%
    dplyr::select(-parent, -node, -branch.length, -species) %>%
    as.list() %>%
    lapply(greta::as_data)

  draws <- list()
  for(i in 1:length(y)) {
    message('Running model on trait "', names(y)[i], '"')
    #draws[[i]] <- run_greta_ridge(y[[i]], design, n_samples = n_samples, thin = thin, warmup = warmup, chains = chains, n_cores = n_cores, ...)
    draws <- run_greta_ridge(y[[i]], design, n_samples = n_samples, thin = thin, warmup = warmup, chains = chains, n_cores = n_cores)
    names(draws) <- names(y)
  }



}

run_greta_ridge <- function(y, design, n_samples, thin, warmup, chains, n_cores, ...) {
  root_value <- greta::normal(0, 10)
  tau <- greta::inverse_gamma(1, 1)
  coefs <- greta::normal(0, tau, dim = ncol(design))

  sd <- greta::cauchy(0, 3, truncation = c(0, Inf))

  mu <- root_value + greta::`%*%`(design, coefs)

  greta::distribution(y) <- greta::normal(mu, sd)

  m <- greta::model(root_value, coefs, sd, tau)

  draws <- greta::mcmc(m, n_samples = n_samples, thin = thin, warmup = warmup, chains = chains, n_cores = n_cores, ...)

  draws
}

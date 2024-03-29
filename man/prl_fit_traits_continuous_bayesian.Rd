% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_models.R
\name{prl_fit_traits_continuous_bayesian}
\alias{prl_fit_traits_continuous_bayesian}
\title{Title Fit a continuous model of trait evolution along a phylogeny, using a Bayesian method}
\usage{
prl_fit_traits_continuous_bayesian(x, family = c("gaussian"),
  method = c("Ridge"), tree = NULL, regularise_first_split = FALSE,
  n_samples = 1000, thin = 1, warmup = 1000, chains = 4,
  n_cores = NULL, ...)
}
\arguments{
\item{x}{Object of class \code{phurl}, such as that produced by \code{\link{phurl_prepare_data}},
or a matrix or data.frame of traits values to be fit (with rownames or order matching tree$tip.label).}

\item{family}{a character vector specifying the statistical families to be used for each trait. If named,
names will be used to match the family to a column in the traits data, otherwise they will be used in order.
If the length of the vector is less than the number of traits, it will be recycled with a warning.}

\item{method}{Method to use. Currently the only option is "Ridge". This is mainly for compatibility
with \code{\link{fit_trait_shift}}, which allows multiple methods.}

\item{tree}{A phylogeny of class \code{phylo}. Only used if x is a matrix of trait values (if x is an object of
class \code{phurl} it should have a tree already associated with it)}

\item{regularise_first_split}{Should the edges beginning at the first split in the phylogeny be regularised?
If not, these are treated like a a set of 'intercepts' for evolutionary rates, and every subsequent change
in these rates along the tree are regularised (penalised).}

\item{n_samples}{Number of MCMC samples to draw during the fitting process.}

\item{thin}{Thinning rate for MCMC samples.}

\item{warmup}{Number of iteration to use for warmup of the MCMC sampling process.}

\item{chains}{Number of replicate MCMC chains to run.}

\item{n_cores}{Number of CPU cores to use for sampling.}

\item{...}{Any other parameters to pass to the underlying fitting function (in this case \code{\link[greta]{mcmc}})}
}
\value{
Object of class \code{phurl}, containing a fitted model in the \code{model} slot, and a
tibble in the \code{results} slot containing estimates of rates along each edge, estimates of
rate changes at each node, and ancestral character estimates (ACEs), including samples from the full
posterior distribution.
}
\description{
This function fits a Bayesian version of Ridge regression on edge paths, producing a phylogenetically
smoothed estimate of trait evolution rates along a phylogeny's branches. Shrinkage limits how much
change in rates can occur at each split of the phylogeny.
}
\details{
Test equations:
\ifelse{html}{\out{<pre>y<sub>i</sub> ~ f(μ, φ)

             ===
             \
  g(μ) = α + /   β<sub>j</sub>
             ===
               j

β<sub>j</sub> ~ Normal(0, τ)

τ ~ InvGamma(1, 1)</pre>}}{\out{<pre>y<sub>i</sub> ~ f(mu, phi)

                  ===
                  \
  g(mu) = alpha + /   beta<sub>j</sub>
                  ===
                    j

beta<sub>j</sub> ~ Normal(0, tau)

tau ~ InvGamma(1, 1)</pre>}}
}

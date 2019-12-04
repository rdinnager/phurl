test_that("multiplication works", {
  traits <- prl_sim_traits_continuous(RRphylo::DataCetaceans$treecet, n_trait = 5, tau = 0.2, initial_rate = 0.5, sd = 0)
  phurl_ob <- prl_prepare_data(RRphylo::DataCetaceans$treecet, traits = traits$traits, species_column = "species",
                               standardise_traits = FALSE)
  with_plot <- prl_plot_traits(phurl_ob)
})

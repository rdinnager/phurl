test_that("different classes of trait data can be loaded and are the same", {
  tree <- RRphylo::DataOrnithodirans$treedino
  traits_vec <- RRphylo::DataOrnithodirans$massdino
  traits_vec_no_names <- unname(traits_vec)
  traits_mat_no_names <- matrix(traits_vec, ncol = 1) 
  colnames(traits_mat_no_names) <- "trait"
  traits_mat_names <- traits_mat_no_names
  rownames(traits_mat_names) <- names(traits_vec)
  traits_dataframe_names <- as.data.frame(traits_mat_names)
  traits_tibble_no_names_col <- traits_dataframe_names %>%
    dplyr::mutate(species = names(traits_vec)) 
  traits_dataframe_names_col <- traits_tibble_no_names_col %>%
    as.data.frame()
  
  expect_warning(ob1 <- phurl_prepare_data(tree, traits_vec_no_names))
  expect_warning(ob2 <- phurl_prepare_data(tree, traits_mat_no_names))
  
  expect_identical(ob1, ob2)
  
  ob_vec <- phurl_prepare_data(tree, traits_vec)
  ob_mat <- phurl_prepare_data(tree, traits_mat_names)
  ob_df_no_col <- phurl_prepare_data(tree, traits_dataframe_names)
  ob_df_col <- phurl_prepare_data(tree, traits_dataframe_names_col, species_column = "species")
  ob_tibble_col <- phurl_prepare_data(tree, traits_tibble_no_names_col, species_column = "species")
  
  expect_identical(ob_vec, ob_mat)
  expect_identical(ob_vec, ob_df_no_col)
  expect_identical(ob_vec, ob_df_col)
  expect_identical(ob_vec, ob_tibble_col)
  expect_identical(ob_mat, ob_df_no_col)
  expect_identical(ob_mat, ob_df_col)
  expect_identical(ob_mat, ob_tibble_col)
  expect_identical(ob_df_no_col, ob_df_col)
  expect_identical(ob_df_no_col, ob_tibble_col)
  
  expect_s3_class(ob_vec, "phurl")
  
  expect_null(ob_vec$model_fit)
  expect_null(ob_vec$plots)
  
  expect_equal(length(ob_vec$nodes_on_path_to_tips), ape::Ntip(ob_vec$phylo))
  expect_equal(length(ob_vec$nodes_on_path_to_nodes), ape::Nnode(ob_vec$phylo))
  
  expect_equal(ncol(ob_vec$x_data_tips), ape::Ntip(ob_vec$phylo) + ape::Nnode(ob_vec$phylo) + 1)
  expect_equal(ncol(ob_vec$x_data_nodes), ape::Nnode(ob_vec$phylo) + 1)
  
})

prl_plot_traits <- function(x, palette = "berlin") {
  if(is.null(x$y_data)) {
    stop("phurl object does not contain trait data")
  }
  
  tree_dat <- tidytree::as_tibble(x$phylo_unmodified) %>%
    dplyr::left_join(x$y_data %>%
                       dplyr::select(-parent, -node, - branch.length), 
                     by = c("label" = "species")) %>%
    tidytree::as.treedata()
  
  p <- ggtree::ggtree(tree_dat) +
    ggtree::theme_tree2()
  
  trait_df <- tidytree::as_tibble(tree_dat) %>%
    dplyr::filter(node <= ape::Ntip(tree_dat@phylo))
    
  tip_names <- trait_df$label
  
  trait_df <- trait_df %>%
    dplyr::select(-parent, -node, -branch.length, -label) %>%
    as.data.frame()
  
  rownames(trait_df) <- tip_names
  
  
  p <- ggtree::gheatmap(p, trait_df) +
    ggtree::scale_x_ggtree() + 
    ggplot2::scale_y_continuous(expand=c(0, 0)) +
    scico::scale_fill_scico(palette = palette)
  
  plot(p)
  
  if(is.null(x$plots)) {
    x$plots <- list(trait_plot = p)
  } else {
    x$plots <- c(x_plots, trait_plot = p)
  }
  
  return(x)
  
}
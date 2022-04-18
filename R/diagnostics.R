#' @name get_avg_nodes
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Gets node averages
#' @description Get node averages per tree
#' @param model The HEBART model object
#' @return The correspondent averages

get_avg_nodes <- function(model){
  all_trees <- model$final_trees$tree_data

  get_parent <- function(data){
    dplyr::pull(data, parent)
  }

  avg_nodes <- function(trees){
    trees %>%
      dplyr::slice(-1) %>%
      dplyr::mutate(nodes   = purrr::map(est_tree, get_parent),
                    n_nodes = purrr::map_dbl(nodes, dplyr::n_distinct)) %>%
      dplyr::pull(n_nodes) %>%
      mean()
  }

  all_avg <- all_trees %>%
    purrr::map_dbl(avg_nodes)
  return(all_avg)
}


#' @name plot_avg_nodes
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Plots node averages
#' @description Gets node averages per tree and plots it
#' @param model The HEBART model object
#' @return The correspondent averages

plot_avg_nodes <- function(model){
  avg_nodes <- dplyr::tibble(
    `Node Average` = get_avg_nodes(model)
  ) %>%
    dplyr::mutate(`Tree Index` = 1:n())

  ggplot2::ggplot(data = avg_nodes, ggplot2::aes(y = `Node Average`, x = `Tree Index`)) +
    ggplot2::geom_hline(yintercept = mean(avg_nodes$`Node Average`),
                        colour = '#c95a49', size = 0.5, linetype = 'dotted') +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::theme_bw(12)
}


#' @name diagnostics_density_plot
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Density plots for tau, k1 or sqrt(k1/tau)
#' @description Plots the sampled of values as a density
#' @param model The HEBART model object
#' @param type Type of plot: tau, k1 or sqrt(k1/tau)
#' @param sqrt Logical to decide whether to plot 1/sqrt(tau) instead
#' (residual precision)
#' @return The correspondent plot

diagnostics_density_plot <- function(model, type = 'tau', sqrt = FALSE){
  if(type == 'tau'){
    df_tau <- data.frame(tau = model$tau_post)
    if(sqrt){
      df_tau$tau <- 1/sqrt(df_tau$tau)
      label_x <- expression('Samples of 1/sqrt('~tau~')')
    } else {
      label_x <- expression('Samples of '~tau)
    }

    ggplot2::ggplot(df_tau, ggplot2::aes(x = tau)) +
      ggplot2::geom_vline(xintercept = mean(df_tau$tau),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::labs(x = label_x, y = "Density") +
      ggplot2::theme_bw(12)
  } else if( type == 'k1'){

    df_k1   <- data.frame(k1 = model$sampled_k1)
    label_x <- expression('Samples of '~k[1])

    ggplot2::ggplot(df_k1, ggplot2::aes(x = k1)) +
      ggplot2::geom_vline(xintercept = mean(df_k1$k1),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(x = label_x, y = "Density") +
      ggplot2::theme_bw(12)
  } else if( type == 'k1_tau'){

    df   <- data.frame(k1 = model$sampled_k1,
                       tau = model$tau_post) %>%
      dplyr::mutate(value = sqrt(k1/tau))
    label_x <-  expression('Samples of sqrt('~k[1]~'/'~tau~')')

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_vline(xintercept = mean(df$value),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(x = label_x, y = "Density") +
      ggplot2::theme_bw(12)
  } else{
    stop("Type of plot not available")
  }
}


#' @name diagnostics_traceplot
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Traceplot for tau, k1 or sqrt(k1/tau)
#' @description Plots the sampled of values of tau, k1 or sqrt(k1/tau) per iteration
#' @param model The HEBART model object
#' @param type Type of plot: tau, k1 or sqrt(k1/tau)
#' @param sqrt Logical to decide whether to plot 1/sqrt(tau) instead
#' (residual precision)
#' @return The correspondent plot

diagnostics_traceplot <- function(model, type = 'tau', sqrt = FALSE){
  if(type == 'tau'){
    df_tau <- data.frame(tau = model$tau_post) %>%
      dplyr::mutate(iter = 1:n())
    if(sqrt){
      df_tau$tau <- 1/sqrt(df_tau$tau)
      label_y <- expression('Samples of 1/sqrt('~tau~')')
    } else {
      label_y <- expression('Samples of '~tau)
    }

    ggplot2::ggplot(df_tau, ggplot2::aes(y = tau, x = iter)) +
      ggplot2::geom_hline(yintercept = mean(df_tau$tau),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(y = label_y, x = "Iteration") +
      ggplot2::theme_bw(12)

  } else if( type == 'k1'){

    df_k1   <- data.frame(k1 = model$sampled_k1) %>%
      dplyr::mutate(iter = 1:n())
    label_y <- expression('Samples of '~k[1])

    ggplot2::ggplot(df_k1, ggplot2::aes(y = k1, x = iter)) +
      ggplot2::geom_hline(yintercept = mean(df_k1$k1),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(y = label_y, x = "Iteration") +
      ggplot2::theme_bw(12)

  } else if( type == 'k1_tau'){

    df   <- data.frame(k1 = model$sampled_k1,
                       tau = model$tau_post) %>%
      dplyr::mutate(value = sqrt(k1/tau), iter = 1:n())
    label_y <-  expression('Samples of sqrt('~k[1]~'/'~tau~')')

    ggplot2::ggplot(df, ggplot2::aes(y = value, x = iter)) +
      ggplot2::geom_hline(yintercept = mean(df$value),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(y = label_y, x = "Iteration") +
      ggplot2::theme_bw(12)

  } else{
    stop("Type of plot not available")
  }
}

#' @name plot_mse_iter
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title MSE per iteration
#' @description Plots training MSE values per iteration
#' @param model The HEBART model object
#' @return The correspondent plot

plot_mse_iter <- function(model){

  get_mse_tree <-  function(tree){
    tree %>%
      dplyr::summarise(mse = mean((y - sampled_mu_j)^2)) %>%
      dplyr::pull(mse)
  }

  get_mse <- function(tree_data){
    all_trees <- tree_data$est_tree
    all_trees <- all_trees[-1]
    mses <- all_trees %>% purrr::map_dbl(get_mse_tree)
    mses
  }

  all_mse <- dplyr::tibble(
    trees = model$final_trees$tree_data
  ) %>%
    dplyr::mutate(mses = purrr::map(trees, get_mse))

  all_mse <- all_mse %>%
    dplyr::mutate(n = 1:n()) %>%
    tidyr::unnest(mses) %>%
    dplyr::group_by(n) %>%
    dplyr::mutate(iter = 1:n()) %>%
    dplyr::ungroup()

  df_avg <- all_mse %>%
    dplyr::select(iter, mses) %>%
    dplyr::group_by(iter) %>%
    dplyr::summarise(avg_mse = mean(mses))

  # Plotting -----
  label_y <-  expression('Average MSE')
  ggplot2::ggplot(df_avg, ggplot2::aes(y = avg_mse, x = iter)) +
    ggplot2::geom_hline(yintercept = mean(df_avg$avg_mse),
                        colour = '#c95a49', size = 0.5, linetype = 'dotted') +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::labs(y = label_y, x = "Iteration") +
    ggplot2::theme_bw(12)

}


#' @name diagnostics
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Plots all diagnostics
#' @description Plots all diagnostics plots
#' @param model The HEBART model object
#' @return The correspondent plot

diagnostics <- function(model){
  # Traceplots --------------------
  p1t <- diagnostics_traceplot(model, type = "tau", sqrt = TRUE)
  p2t <- diagnostics_traceplot(model, type = "k1")
  p3t <- diagnostics_traceplot(model, type = "k1_tau")
  p1t <- p1t +
    ggplot2::ggtitle("Traceplots") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "italic"))

  # Density plots --------------------
  p1d <- diagnostics_density_plot(model, type = "tau", sqrt = TRUE)
  p2d <- diagnostics_density_plot(model, type = "k1")
  p3d <- diagnostics_density_plot(model, type = "k1_tau")
  p1d <- p1d +
    ggplot2::ggtitle("Density plots") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "italic"))

  # Number of nodes ------------------
  p_node <- plot_avg_nodes(model) +
    ggplot2::ggtitle("Nodes per tree (average)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "italic"))

  # MSE per iteration ------------------
  p_mse <- plot_mse_iter(model) +
    ggplot2::ggtitle("MSE per iteration (average)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "italic"))


  (p1t + p2t + p3t) /
    (p1d + p2d + p3d)  / (p_node + p_mse)


}

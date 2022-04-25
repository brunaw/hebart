#' @name faster_det
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title A function to calculate the determinant of M (faster than solve)
#' @param k_1_d The current value of k1
#' @param k_2_d The current value of k2
#' @param M_d The M matrix
#' @return The determinant

faster_det <- function(k_1_d, k_2_d, M_d) {
  n   <- nrow(M_d)
  n_j <- colSums(M_d)
  tMM <- crossprod(x = M_d)
  diag_k1_n <- k_1_d/(1 + k_1_d*n_j)
  #if(diag_k1_n < 1){ diag_k1_n <- matrix(-1) }
  Psi_tilde_inv <- diag(n) - M_d%*%diag(diag_k1_n)%*%t(M_d)


  return(log(k_2_d) + log(1/k_2_d + sum(Psi_tilde_inv)) + sum(log(1 + k_1_d*n_j)))
}

#' @name marginal_dist
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title A function to calculate the marginal distribution of y
#' @param y The y vector
#' @param k_1 The current value of k1
#' @param k_2 The current value of k2
#' @param M The M matrix
#' @param mu_mu The value for mu_mu (usually 0)
#' @param alpha The alpha value from the prior of tau
#' @param beta The beta value from the prior of tau
#' @return The resulting value

marginal_dist <- function(y, k_1, k_2, M, mu_mu, alpha, beta) {
  n <- nrow(M)
  W_0 <- rep(mu_mu, n)
  ymW_0 <- y - W_0

  term_1 <- -(n/2)*log(2*pi)
  term_2 <- - 0.5 * faster_det(k_1_d = k_1, k_2_d = k_2, M_d = M)
  term_3 <- lgamma(n/2 + alpha)

  term_4 <- - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta)
  all <- term_1 + term_3 + term_2 + term_4

  return(all)
}

#' @name calculate_all
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title A function to calculate the likelihood of the fed parameters
#' @param current_tree_mh The current trees set
#' @param new_k1 The current value of k2
#' @param pars The hyperparameters list
#' @return The resulting value
#'
calculate_all <- function(current_tree_mh, new_k1, pars){

  res <- current_tree_mh %>%
    dplyr::select(tree_index, y, sampled_mu_j, group) %>%
    tidyr::pivot_wider(names_from = tree_index, values_from = sampled_mu_j,
                       values_fn = mean) %>%
    dplyr::mutate(pred = rowSums(.[3:ncol(.)]))

    M       <- stats::model.matrix(~ factor(res$group) - 1)
    tot_lik <- marginal_dist(
      y = res$y, k_1 = new_k1, k_2 = pars$k2,
      M = M, mu_mu = pars$mu_mu, alpha = pars$alpha,
      beta = pars$beta)

  return(tot_lik)

}


#' @name lk_ratio_k
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for the sampled k1
#' @description Likelihood ratio of the candidate tree and the previous
#' tree, given that the action step was a prune.
#' @param new_k1 The sampled k1
#' @param k1 The current value of k1
#' @param current_tree The current tree
#' @param pars The full list of hyperparameters
#' @param group The grouping variable
#' @param tau_post The current value of tau
#' @param M The M matrix
#' @return The likelihood ratio


lk_ratio_k <- function(new_k1, k1, current_tree, pars, group,
                       tau_post = tau_post, M){
  beta    <- pars$beta
  alpha   <- pars$alpha
  mu_mu   <- pars$mu_mu
  k2      <- pars$k2
  # Sample new mu ------
  n_nodes <- dplyr::n_distinct(current_tree$node)
  mu_post <- unique(current_tree$mu_sampled)

  if(length(mu_post) == 1){ mu_post <- rep(mu_post, n_nodes) }
  njs <-  current_tree %>%
    dplyr::group_by(node, group) %>%
    dplyr::count()
  split_nodes <- njs %>% split(.$node)

  n_j_means <- current_tree %>%
    dplyr::group_by(node, group) %>%
    dplyr::summarise(mean_y_j = mean(y)) %>%
    dplyr::arrange(group) %>%
    split(.$node)

  mu_js_post <- list()
  nodes_unique <- unique(njs$node)

  for(m in 1:n_nodes){
    mean_mu      <- c()
    var_mu       <- c()

    nj_node      <- split_nodes[[nodes_unique[m]]] %>% dplyr::pull(n)
    y_bar_node   <-  n_j_means[[nodes_unique[m]]] %>% dplyr::pull(mean_y_j)

    for(j in sort(unique(group))){
      y_bar_j    <- y_bar_node[j]
      mean_mu[j] <- ((mu_post[m]/new_k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/new_k1)
      var_mu[j]  <- (tau_post*(nj_node[j] + 1/new_k1))^(-1)
    }
    mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
  }

  new_error <- tidyr::expand_grid(
    group = 1:25,
    node = unique(current_tree$node)) %>%
    dplyr::arrange(node, group) %>%
    dplyr::mutate(muj = c(mu_js_post[[1]], mu_js_post[[2]])) %>%
    dplyr::left_join(dplyr::select(current_tree, y, node, group),
              by = c("node", "group")) %>%
    dplyr::mutate(err_y = (y - muj)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y)) %>%
    dplyr::pull(sum_errors_y)

  # --------------------------------------------------
  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------
  error_y <- current_tree %>%
    dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y)) %>%
    dplyr::pull(sum_errors_y)

  y      <- current_tree$y
  inner  <- error_y/2 +  beta
  M_mat  <- stats::model.matrix(y ~ factor(group) - 1, current_tree)
  N      <- nrow(current_tree)
  W_0    <- rep(mu_mu, N)
  k2_mat <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
               t(diag(x = 1, nrow = N, ncol = 1)))
  #------------------------------------------------------------------
  psi    <- (k1 * M_mat %*% t(M_mat)) + diag(N)
  W_1    <- k2_mat + psi

  p_cond_current <- (-1/2) * log(det(W_1)) + (
    log(inner)*(-(N/2 + alpha)))

  # Candidates ------------------------------
  psi_candidate   <- (new_k1 * M_mat %*% t(M_mat)) + diag(N)
  W_1_candidate   <- k2_mat + psi_candidate
  inner_new       <-  new_error/2 + beta
  p_cond_new      <- (-1/2)*log(det(W_1_candidate)) +(
    log(inner_new)*(-(N/2 + alpha)))

  lk_ratio        <- p_cond_new - p_cond_current
  r               <- min(1, exp(lk_ratio))
  return(r)
}

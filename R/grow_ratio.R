#' @name grow_tree
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Grows the current tree
#' @param current_tree The current tree
#' @param selec_var The node to split on node
#' @param drawn_node The node to split on
#' @param rule The variable value to split on
#' @return The new tree
grow_tree <- function(current_tree, selec_var, drawn_node, rule){

  # Is this a binary split?
  available_splits <- current_tree %>% dplyr::pull(selec_var) %>% unique()

  if(length(available_splits) == 2 & rule == max(available_splits)){
    rule <- min(available_splits)
  }

  # is there data in both sides of the split?
  dd <- current_tree %>%
    dplyr::filter(node == drawn_node) %>%
    dplyr::mutate(cond = !!rlang::sym(selec_var) > rule) %>%
    dplyr::count(cond)

  if(nrow(dd) == 1){
    return(current_tree)
  } else{
    current_tree %>%
      dplyr::mutate(
        # Increasing the depth of the node
        d =  ifelse(node == drawn_node, d + 1, d),
        # Updating the parent of the split node
        parent = ifelse(node == drawn_node, drawn_node, parent),
        # Changing the node "side" of each observation: left and right
        criteria = ifelse(
          node == drawn_node,
          ifelse(!!rlang::sym(selec_var) > rule, "left", "right"), "no split"),

        # Updating the node with the new split
        node = ifelse(node == drawn_node,
                      ifelse(!!rlang::sym(selec_var) > rule,
                             paste(node, selec_var, "left"),
                             paste(node, selec_var, "right")), node),
        # Updating the node index
        node_index =  as.numeric(as.factor(node)))
  }
}

#' @name transition_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Transition ratio for grow.
#' @description Transition probability of going to the candidate tree,
#' given that the action step was a grow
#' @param tree The current tree
#' @param current_node The current node
#' @param p_vars The number of predictors still available
#' @param current_selec_var The variable selected for the split
#' @param p_grow The probability of growing the tree
#' @param i The current iteration
#' @return The transition ratio
#' @details When transitioning to a grow, we need the probabilities of:
#'  1. Growing the tree
#'  2. Growing from the specific node, that has to consider the
#'  available predictors and available values for that grow.
#' When transitioning from a grow back to the split,
#' we need the probabilities of:
#' 1. Pruning the tree
#' 2. Selecting  node to prune
transition_ratio_grow <- function(tree, current_node, p_vars,
                                  current_selec_var,
                                  #results_f,
                                  p_grow, i){

  # #p_grow = 0.5
  # p_grow = 0.95
  p_prune = 1 - p_grow
  # Number of available final nodes to split on -------
  b <-  length(unique(tree$node_index))

  # Probability of splitting (n variables) ------------
  p_adj <-  1/p_vars

  # Available values to split given the variable selected for the
  # split and the node selected ----------------------
  n_j_adj <-  tree %>%
    dplyr::filter(parent == current_node) %>%
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()


  # Probability of the transition -------------------
  p_tstar_to_t <- p_grow * (1/b) * (p_adj) * (1/n_j_adj)


  #  Probability of transitioning from the new tree
  # back to the original -------------------------
  p_prune_the_node <- 1/(b+2)
  p_t_to_tstar <- p_prune * p_prune_the_node

  trans_ratio <- log(p_tstar_to_t/p_t_to_tstar)
  return(trans_ratio)
}

#' @name inv2
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title A function to calculate the inversion of the M matrix
#' @param k_1 The current value of k1
#' @param k_2 The current value of k2
#' @param M The M matrix
#' @return The inverted matrix
inv2 <- function(k_1, k_2, M) {
  n <- nrow(M)
  n_j <- colSums(M)
  Psi_tilde_inv <- diag(n) - M%*%diag(k_1/(1 + k_1*n_j))%*%t(M)
  k_3 <- 1/k_2 + sum(Psi_tilde_inv)
  Psi_row_sums <- rowSums(Psi_tilde_inv)
  Psi_col_sums <- colSums(Psi_tilde_inv)
  inv_M <- Psi_tilde_inv - tcrossprod(Psi_row_sums, Psi_col_sums) / k_3
  return(inv_M)
}

#' @name inv2
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Full parent conditional
#' @param data_cond The data to calculate the conditional distribution on
#' @param pars The hyperparameters list
#' @return The conditional probability
cond_calculation <- function(data_cond, pars){

  n_groups <- length(unique(factor(data_cond$group)))
  if(n_groups > 1){
    M     <-  stats::model.matrix(~ factor(data_cond$group) - 1)
  } else{
    M <- matrix(rep(1, nrow(data_cond), ncol = 1))

  }
  y     <-  data_cond$y
  k_1   <-  pars$k1
  k_2   <-  pars$k2
  mu_mu <-  pars$mu_mu
  alpha <-  pars$alpha
  beta  <-  pars$beta
  n     <- nrow(M)
  W_0   <- rep(mu_mu, n)
  ymW_0 <- y - W_0

  term_1 <- -(n/2)*log(2*pi)
  if(n_groups > 1){
    term_2 <- - 0.5 * log(faster_det(k_1_d = k_1, k_2_d = k_2, M_d = M))
  } else {
    W_1 <- k_1 * M %*% t(M) + diag(n) + k_2
    term_2 <- - 0.5 * log(det(W_1))
  }
  term_3 <- lgamma(n/2 + alpha)
  if(n_groups > 1){
    inv_2 <- inv2(k_1, k_2, M)
  } else {
    inv_2 <- solve(W_1)
  }
  term_4 <- - (n/2 + alpha)*log((0.5 * t(ymW_0)%*%inv_2%*%ymW_0) + beta)
  p_cond_y <- term_1 + term_4 + term_2 + term_3
  return(p_cond_y)
}

#' @name lk_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for grow
#' @description Likelihood ratio of the candidate tree and the previous
#' tree, given that the action step was a grow
#' @param tree The current tree
#' @param current_node The current pruned node
#' @param pars The full list of parameters
#' @return The likelihood ratio
#' @details The likelihood ratio for growth needs the
#' joint distribution of each region (of each node), given the
#' tau parameter
lk_ratio_grow <- function(tree, current_node, pars){

  filtered_tree <- tree %>% dplyr::filter(parent == current_node)
  cond_parent <- cond_calculation(data_cond = filtered_tree, pars = pars)

  nam <- unique(filtered_tree$node)
  cond_node <- 0

  for(name in  nam){
    data_set <- tree %>% dplyr::filter(node == name)
    marg <- cond_calculation(data_cond = data_set, pars = pars)
    cond_node <- marg + cond_node
  }

  lk_ratio <-  cond_node - cond_parent

  return(lk_ratio)
}

#' @name structure_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree structure ratio for grow
#' @description Tree structure ratio of the candidate tree and
#' the previous tree, given that the action step was a grow
#' @param tree The current tree
#' @param current_node The current grown node
#' @param current_selec_var The variable selected for the split
#' @param p_vars The number of available predictors.
#' @param alpha_grow The alpha of the growing probability
#' @param beta_grow The beta of the growing probability
#' @return The tree structure ratio
#' @details For the tree structure ratio of the new tree, we need
#' the probabilities of:
# 1. Splitting at node n
# 2. Splitting at the node of the left
# 3. Splitting at the node of the right
# 4. Using each rule at node n

structure_ratio_grow <- function(tree, current_node,
                                 current_selec_var, p_vars,
                                 alpha_grow, beta_grow){

  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p_adj <- 1/p_vars

  # Counting the distinct rule options from
  # this available predictor -------------------------------
  filter_tree <- tree %>%
    dplyr::filter(parent == current_node)
  #depth <- unique(filter_tree$d)

  n_j_adj <-  filter_tree %>%
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()

  # Calculating the probability of the chosen rule --------
  p_rule <- p_adj * (1/n_j_adj)

  # Calculating the probability of split
  #terminal_nodes <- tree %>% dplyr::distinct(node_index) %>% nrow()

  # Counting terminal & internal nodes
  internal_nodes <- unique(tree$parent)
  n_int          <- length(internal_nodes)
  log_p_inter    <- 0

  for(i in 1:n_int){
    depth_1 <- tree |>
      dplyr::filter(parent == internal_nodes[i]) |>
      dplyr::pull(d) |>
      unique()
    log_p_inter <- log_p_inter + log(alpha_grow) - beta_grow * log(depth_1)
  }

  terminal_nodes <- unique(tree$node)
  n_term         <- length(terminal_nodes)
  log_p_term     <- 0

  for(i in 1:n_term){
    depth <- tree |>
      dplyr::filter(node == terminal_nodes[i]) |>
      dplyr::pull(d) |>
      unique()
    log_p_inter <- log_p_inter + log(alpha_grow) - beta_grow * log(1 + depth)
  }

  p_tree <- exp(log_p_term + log_p_inter)
  #p_right <- p_left <- alpha_grow*(1 + depth + 1)^(-beta_grow)

  # !!! this is a hack, that should be fixed to no split
  # if(internal_nodes == 1){
  #   p_split <-  0.99
  # } else {
  #   p_split <- 1/terminal_nodes
  # }

  p_t_star <-  p_tree*p_rule
  p_t      <-  (1 - exp(log_p_inter))

  st_ratio <- log(p_t_star/p_t)

  return(st_ratio)
}

#' @name tree_prior
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree prior
#' @description Calculates the prior for each tree structure
#' @param tree The current tree
#' @param alpha_grow The alpha of the growing probability
#' @param beta_grow The beta of the growing probability
#' @return The tree prior
tree_prior <- function(tree, alpha_grow, beta_grow){

  # Counting terminal & internal nodes
  internal_nodes <- unique(tree$parent)
  n_int          <- length(internal_nodes)
  log_prior      <- 0

  for(i in 1:n_int){
    depth_1 <- tree |>
      dplyr::filter(parent == internal_nodes[i]) |>
      dplyr::pull(d) |>
      unique()
    log_prior <- log_prior + log(alpha_grow) - beta_grow * log(depth_1 + 1)
  }

  terminal_nodes <- unique(tree$node)
  n_term         <- length(terminal_nodes)

  for(i in 1:n_term){
    depth <- tree |>
      dplyr::filter(node == terminal_nodes[i]) |>
      dplyr::pull(d) |>
      unique()
    log_prior <- log_prior + log(alpha_grow) - beta_grow * log(1 + depth + 1)
  }

  return(log_prior)
}


#' @name ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Final ratio for a growth step.
#' @description The final ratio is to be used as the acceptance
#' criteria in the MCMC of the hebart model
#' @param tree The current tree
#' @param old_tree The old tree
#' @param current_node The current grown node
#' @param pars The full list of parameters
#' @param alpha_grow The alpha of the growing probability
#' @param beta_grow The beta of the growing probability
#' @return The final ratio for the candidate tree

ratio_grow <- function(tree, old_tree,
                       current_node,
                       pars,
                       alpha_grow, beta_grow){

  # All ratios -- old values
  # trans <- transition_ratio_grow(tree, current_node,
  #                                current_selec_var = current_selec_var,
  #                                p_vars = p_vars, i = i,
  #                                p_grow = p_grow)
  #
  # struct <- structure_ratio_grow(tree, current_node,
  #                                current_selec_var = current_selec_var,
  #                                p_vars = p_vars,
  #                                alpha_grow, beta_grow)

  lk             <- lk_ratio_grow(tree, current_node, pars)
  new_tree_prior <- tree_prior(tree, alpha_grow, beta_grow)
  old_tree_prior <- tree_prior(old_tree, alpha_grow, beta_grow)
  pr_ratio       <- new_tree_prior - old_tree_prior

  r <- min(1, exp(lk + pr_ratio))
  return(r)
}

#' @name hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Hierarchical Embedded Bayesian Additive Regression Trees
#' @description This function runs a BCART model and returns the tree and
#' other results obtained in the last iteration of the MCMC
#' @param formula The model formula
#' @param dataset The data to be used in the modeling
#' @param iter The number of MCMC iterations
#' @param group_variable The grouping variable
#' @param pars The hyperparameters set/list
#' @param min_u Integer representing the lower interval of the
#' Uniform distribution used to sample k1
#' @param max_u Integer representing the upper interval of the
#' @param prior_k1 Logical to decide whether or not use a prior for k1
#' Uniform distribution used to sample k1
#' @param P The number of trees
#' @param sample_k1 Logical to decide whether to sample_k1 or not
#' @param burn_in The number of burn-in iterations
#' @param alpha_grow Number between 0 and 1 used in the growing probability
#' calculation
#' @param beta_grow Number between 0 and 1 used in the growing probability
#' calculation
#' @param ... Other parameters
#' @return A list containing:
#'  the sampled values of tau and k1, the final trees
#' @details
#' Priors used ----------------------------------------------------------
#' y_{ij} ~ Normal(m_j, tau^-{1})
#' tau    ~ Gamma(alpha, beta)
#' mu     ~ Normal(0, tau_mu = k2*tau^-{1})
#' mu_j   ~ Normal(mu, k1*tau^-{1})
#' ----------------------------------------------------------------------

hebart <- function(formula, dataset, iter = 500,
                   group_variable = 'group',
                   pars,
                   # HEBART parameters
                   min_u = 0, max_u = 20, prior_k1 = TRUE,
                   P = 5, sample_k1 = TRUE, burn_in = 50,
                   alpha_grow = 0.95,
                   beta_grow  = 0.15,
                   ...
){

  options(dplyr.summarise.inform = FALSE)

  # HEBART Function
  #---------------------------------------------------------------------
  # Handling initial dataset
  #---------------------------------------------------------------------
  results_data <- data_handler(formula, data = dataset, group_variable)
  data         <- results_data$data
  depara_names <- results_data$names
  names(data)[names(data) == group_variable] <- "group"
  N <- n       <- nrow(data)
  group        <- results_data$group
  mf           <- stats::model.frame(formula, dataset)
  y            <- stats::model.extract(mf, "response")

  name_y       <- names(mf)[1]
  names(data)[names(data) == name_y] <- "y"

  #---------------------------------------------------------------------
  # Defining current distribution parameters
  #---------------------------------------------------------------------
  # Prior hyperparameters -------------
  J <- dplyr::n_distinct(group)
  beta <- pars$beta
  alpha <- pars$alpha
  mu_mu <- pars$mu_mu
  k1 <- pars$k1
  k2 <- pars$k2
  # Minimum batch for each node
  keep_node <- 0.05 * nrow(data)
  x_vars <-  all.vars(formula[[3]])
  p_vars <- length(x_vars)  # number of vars
  to_do <- vector()         # Actions that can be taken
  #---------------------------------------------------------------------
  # Initializing useful vectors
  #---------------------------------------------------------------------
  # For grow or prune
  action_taken = vector()       # To save the actions taken in the algorithm
  selec_var = vector()          # To save the selected variable when growing
  rule = vector()               # To save the selected splitting rule when growing
  drawn_node = vector()         # To save the selected node to grow or prune

  r <- r_k <- vector()          # To save the ratios of grow or prune
  u <- u_k <-  vector()         # To save the sampled uniform values
  sampled_k1 <- vector()        # To save the sampled values for k1
  sampled_k1[1] <- pars$k1

  # For the trees  ------------
  my_trees_l <-  list()          # To save each new tree
  my_trees_l[[1]] <- data        # Initializing the first tree as the
                                 # data will be the  'root' tree, with no nodes

  # For the sampling of posterior values  -------------
  tau_post <- vector()            # To save posterior values of
  tau_post[1] <- stats::rgamma(n = 1, 1/alpha, beta)
  parent_action <- vector()
  results <- data.frame(node = NA, var = NA, rule = NA, action = NA)

  my_trees <- dplyr::tibble(est_tree = my_trees_l, parent_action = NA)

  # One results and one tree_data per tree
  all_trees <- dplyr::tibble(tree_index = 1:P,
                             tree_data = list(my_trees),
                             results = list(results))

  #---------------------------------------------------------------------
  # A simple progress bar
  pb <- progress::progress_bar$new(
    format = "  Iterations of the HEBART model [:bar]  :current/:total (:percent)",
    clear = FALSE, width = 60, total = iter)

  # Loop to perform the HEBART model
  for(i in 1:iter){

  for(p in 1:P){

    # Growing, pruning or staying on the same tree ---------------------
    # Unnesting trees data by i and p
    current_tree <- tidyr::unnest(all_trees[p, ], tree_data) %>%
      dplyr::slice(i)

    my_tree <- tidyr::unnest(current_tree, est_tree) %>%
      dplyr::select(dplyr::starts_with("X"), y, node, d, group, parent, node_index)

    results_current <- tidyr::unnest(
      dplyr::select(current_tree, results), results)
    #----------------------------------------------------------------------
    # Sampling details ----------------------------------------------------
    action_taken = "grow"

    dn     <- dplyr::n_distinct(my_tree$d)
    p_grow <- alpha_grow*(1 + dn)^(-beta_grow)
    u_grow <- stats::runif(1)
    depth  <- dplyr::n_distinct(my_tree$node)

    # deciding on growing or pruning
    if(u_grow > p_grow & depth > 1){
      action_taken = 'prune'
    } else{
      action_taken = 'grow'
    }

    if(action_taken == 'grow'){
      # Selecting the node to grow, uniformly
      drawn_node <- sample(unique(my_tree$node), size = 1)
      # Selecting the variable and splitting rule, uniformly
      selec_var <- depara_names$new[sample(1:p_vars, size = 1)]
      rule  <- p_rule(variable_index = selec_var,
                      data = my_tree, sel_node = drawn_node)


        # all_actions <- paste0(
        #   results_current$var, " ",
        #   results_current$rule
        # )
        #
        # current_action <- paste0(selec_var, " ", rule)
        # if(current_action %in% all_actions){
        #   sample_tree <- my_tree
        #   parent_action <- NA
        # } else{

          if(is.na(rule)){
            sample_tree <- my_tree
            parent_action <- NA
          } else {
            # Grow the tree
            sample_tree <- grow_tree(
              current_tree = my_tree, selec_var = selec_var,
              drawn_node = drawn_node, rule = rule
            )
          }
        #}

      # Checking whether all of the nodes have the minimum
      # of observations required -- per group as well
      temps <-  sample_tree  %>%
        dplyr::count(node) %>%
        dplyr::pull(n)

      # temps_g <-  sample_tree  %>%
      #   dplyr::group_by(node) %>%
      #   dplyr::summarise(n_group = n_distinct(group),
      #                    n_group_b = n_group < J ,
      #                    sum_b = sum(n_group_b))


      #if(sum(temps <= keep_node) > 0 | sum(temps_g$sum_b) > 0){
      if(sum(temps <= keep_node) > 0){
        # If all nodes don't have the minimum of observations,
        # just keep the previous tree
        sample_tree <- my_tree
        parent_action <- NA
        #current_results <- results[[1]]
      } else {

        # Saving the parent of the new node to use in the
        # calculation of the transition ratio
        # if(is.na(rule)){
        #   sample_tree <- my_tree
        #   parent_action <- NA
        #   #r <- 0
        # } else{
        if(!identical(my_tree, sample_tree)){
        parent_action <- sample_tree %>%
          dplyr::filter(parent == drawn_node) %>%
          dplyr::pull(parent) %>%
          unique()

          results_new <- suppressWarnings(
            dplyr::bind_rows(results_current,
                             data.frame(
                               node = drawn_node,
                               var = selec_var,
                               rule = rule,
                               action = action_taken)))



        # Calculating the acceptance ratio for the growing of this tree,
        # which uses the ratios of:
        # 1. The transition probabilities of the trees
        # 2. The likelihoods of the trees
        # 3. The probability of the tree structures

        r <- ratio_grow(tree = sample_tree,
                        current_node =  parent_action,
                        pars = pars,
                        p = p_vars,
                        current_selec_var = selec_var,
                        i = i,
                        p_grow = p_grow)

        } else{ r <- 0 }
      }
      # Should we prune the tree?
    } else {
      #(to_do[i] >  p_grow && depths > 0){
      #action_taken[i] = "prune"

      # Selecting node to prune, uniformly
       drawn_node <- sample(unique(my_tree$node), size = 1)

      parent_action <- my_tree %>%
        dplyr::filter(node == drawn_node) %>%
        dplyr::distinct(parent) %>%
        dplyr::pull(parent)

      #parent_prune <- stringr::str_remove(parent_action, '( right| left)$')

      # results_current <- results_current %>%
      #   dplyr::filter(!stringr::str_detect(node, parent_prune))

      #results_current$action <- action_taken

      selec_var <- stringr::str_extract(
        drawn_node,'X[0-9][^X[0-9]]*$') %>% stringr::str_remove(" left| right")

      rule <- NA

      # Prune the tree
      sample_tree <- prune_tree(
        current_tree = my_tree, drawn_node, selec_var
      )
      # table(sample_tree$node)
      # table(my_tree$node)

      results_new <- results_current %>%
      dplyr::filter(!stringr::str_detect(node, parent_action))

      # Calculating the acceptance ratio for the prune,
      # which uses the ratios of:
      # 1. The transition probabilities of the trees
      # 2. The likelihoods of the trees
      # 3. The probability of the tree structures

      r <- ratio_prune(old_tree = my_tree,
                       tree = sample_tree,
                       pars = pars,
                       current_node = parent_action,
                       var_in_prune = selec_var,
                       p_split = 0.5,
                       nodes_to_prune = nodes_to_prune,
                       p_grow = p_grow,
                       i = i)

    }


    #
    #   # Should we stay in the same tree?
    # } else {
    #   action_taken[i] = "same"
    #   my_trees[[i + 1]] <- my_trees[[i]]
    #   results[[i+1]] <- results[[i]]
    #   r[i] <- 0
    # }


    # Checking if the tree will be accepted or not ---------------------
    # TO REVIEW
    if(!identical(
      dplyr::select(my_tree, node, parent),
      dplyr::select(sample_tree, node, parent))){

      # Should we accept the new tree?
      u <- stats::runif(1)

      # Checking if the tree should be accepted or not, based
      # on the acceptance ratio calculated and a value sampled
      # from a uniform distribution
      if(u >= r){
        # If that, do not accept tree
        sample_tree <- my_tree
      } else {
        results_current <- results_new

        #Need to find a better way to do this
        # # Is this the same action as the latest iteration
        # if(nrow(results_current) > 1 && action_taken == 'grow'){
        #
        #   all_actions <- paste0(
        #     results_current$var, " ",
        #     results_current$rule
        #   )
        #   current_action <- paste0(selec_var, " ", rule)
        #
        #   if(current_action %in% all_actions){
        #     sample_tree <- my_tree
        #   } else {
        #     results_current <- results_new
        #   }
        # } else {
        #   results_current <- results_new
        # }
    }

    }

    # Updating posteriors -----------
    if(p == 1){
      sample_tree$res <- sample_tree$y
    } else{
      sample_tree$res <- res_previous
      rm(res_previous)
    }

    mu_post <- sample_parameters(
      type = "mu",
      current_tree_post = sample_tree,  P = P,
      k1 = pars$k1, k2 = pars$k2, tau_post = tau_post[i])

    mu_js_post <- sample_parameters(
      type = "muj", current_tree_post = sample_tree,
      k1 = pars$k1, k2 = pars$k2, J = J, P = P, mu_post = mu_post,
      tau_post = tau_post[i])


    sample_tree <- sample_tree %>%
      dplyr::left_join(mu_post, by = "node") %>%
      dplyr::left_join(mu_js_post, by = c("node", "group"))

    # Res = y - sum of all mu_j so far
    sample_tree$res <- sample_tree$res - sample_tree$sampled_mu_j

    # Getting residuals from previous tree
    res_previous <- sample_tree$res

    all_trees[p, ]$tree_data <-
      all_trees[p, ]$tree_data %>%
      purrr::map(~{
        .x %>% add_row(
          est_tree = list(sample_tree),
          parent_action = parent_action
        )})


    all_trees[p, ]$results[[1]] <- results_current

  }


    samp_aux <- all_trees %>%
      dplyr::mutate(current_tree = purrr::map(tree_data, ~{ tail(.x, 1)})) %>%
      dplyr::select(tree_index, current_tree) %>%
      tidyr::unnest(current_tree) %>%
      dplyr::select(tree_index, est_tree) %>%
      tidyr::unnest(est_tree) %>%
      dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>%
      dplyr::ungroup()


    # -----------------------------------------------------------
    # Sampling from the posterior distribution of tau -----------
    tau_post[i + 1] <- sample_parameters(
      type = "tau", i = i,
      k1 = pars$k1,
      k2 = pars$k2,
      P = P,
      current_tree_post = samp_aux,
      alpha = alpha, beta = beta, N = N)

    # min_u = 0
    # max_u = 20
    # sample_k1 = T
    # ------------------------------------------------------------
    # Sampling K1 -----------------------------------------------
    if(sample_k1){
    # We can set these parameters more smartly
    samp_k1 <- MH_update(
      current_tree_mh = samp_aux,
      i = i,
      prior = prior_k1,
      min_u = min_u, max_u = max_u, pars = pars)

    if(samp_k1 == pars$k1){ sampled_k1[i] <- pars$k1
    } else {
      pars$k1 <- samp_k1
      sampled_k1[i] <- samp_k1
    }
    } else{
      sampled_k1[i] <- pars$k1
    }

    pb$tick()
  }



  return(list(
    tau_post = tau_post[-c(1:burn_in)],
    final_trees = all_trees,
    sampled_k1 = sampled_k1[-c(1:burn_in)]
  ))
}


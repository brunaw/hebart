#' @name sample_pred
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Prediction sampling
#' @description This function predicts for the final trees of a
#' Hierachical Embedded BART model
#' @param mu The mu value to use as the mean
#' @param k1 The latest sampled k1 value
#' @param tau The latest sampled tau value
#' @return The sampled values

sample_pred <- function(mu, k1, tau){
  tau_mu <- 1/(tau/k1)
  samp <- stats::rnorm(1, mean = mu, sd = sqrt(tau_mu))
  return(samp)
}

#' @name calculate_preds
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Calculates predictions
#' @description Calculates predictions for the HEBART model
#' @param results The results (description to be completed)
#' @param tree The latest tree results
#' @param newdata The newdata to predict for
#' @param k1 The latest sampled k1 value
#' @param tau The latest sampled tau value
#' @return The sampled values

calculate_preds <- function(results, tree, newdata, k1, tau){

  tree <- tree %>%
    dplyr::select(est_tree) %>%
    tidyr::unnest(cols = est_tree)

  nodes <- unique(tree$node)

  # names(X) <-  paste0("X", 1:ncol(X))
  # data <- data.frame(y, X, group)

  model_muj <- tree %>%
    dplyr::group_by(node, group) %>%
    dplyr::summarise(sampled_mu_j = unique(sampled_mu_j), .groups = 'drop')

  model_mu <- tree %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(sampled_mu = unique(sampled_mu), .groups = 'drop')

  newdata$node <- "root"
  pred <- newdata


  if(length(nodes) == 1 && nodes == "root"){
    pred_final <- pred %>%
      dplyr::left_join(model_mu, by = c("node")) %>%
      dplyr::left_join(model_muj, by = c("node", "group")) %>%
      dplyr::mutate(
        sampled_mu_j =
          ifelse(is.na(sampled_mu_j),
                 #sample_pred(sampled_mu, k1, tau),
                 sampled_mu,
                 sampled_mu_j))

    return(pred_final)

  } else {


  results <- dplyr::filter_all(results, dplyr::any_vars(!is.na(.)))
  results <- results %>% dplyr::mutate(id = 1:n())

  # for each grow followed by prune, remove the parent
  # to_prune <- dplyr::filter(results, action == 'prune')
  #
  #
  # if(nrow(to_prune) > 0){
  #   for(i in 1:nrow(to_prune)){
  #
  #     # results_change <- results %>%
  #     #   dplyr::group_by(node) %>%
  #     #   dplyr::mutate(n_row = 1:n()) %>%
  #     #   dplyr::ungroup()
  #
  #     results_change <- results %>%
  #       filter(id <= to_prune$id[i])
  #
  #     results_rest <- results %>%
  #       filter(id > to_prune$id[i])
  #
  #     results_f <- results_change %>%
  #       filter(node %in% to_prune$node[i]) %>%
  #       # getting the closest node
  #       mutate(diff_id = to_prune$id[i] - id) %>%
  #       dplyr::mutate(to_filter = diff_id %in% c(0, sort(diff_id)[2])) %>%
  #       dplyr::filter(to_filter)
  #
  #     also_to_filter <- results_change %>%
  #       rowwise() %>%
  #       mutate(to_filter = str_detect(node, pattern = to_prune$node[i])) %>%
  #       dplyr::filter(to_filter)
  #
  #     results_change <- results_change %>%
  #       filter(!id %in% c(results_f$id, also_to_filter$id))
  #
  #     results <- bind_rows(results_change, results_rest)
  #
  #   }
  #
  #   #results <- dplyr::select(results, -to_filter, -n_row)
  # }

  # create nodes in new data then left join with main trees
  for(i in 1:nrow(results)){
    pred <- pred %>%
      dplyr::mutate(
        node =
          ifelse(
            node == results$node[i],
            ifelse(!!rlang::sym(results$var[i]) > results$rule[i],
                   paste(node, results$var[i], "left"),
                   paste(node, results$var[i], "right")), node))
  }

  pred_final <-   pred %>%
    dplyr::left_join(model_mu, by = c("node")) %>%
    dplyr::left_join(model_muj, by = c("node", "group")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sampled_mu_j =
        ifelse(is.na(sampled_mu_j),
               #sample_pred(sampled_mu, k1, tau),
               sampled_mu,
               sampled_mu_j)) %>%
    dplyr::ungroup()


  return(pred_final)

  }


}

#' @name predict_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Predictions for the HEBART model
#' @description This function predicts for the final trees of a
#' Hierachical Embedded BART model
#' @param model The model object
#' @param newdata The new data to predict
#' @param formula The model formula
#' @param group_variable The grouping variable
#' @return A dataframe with the prediction column

predict_hebart <- function(model, newdata, formula, group_variable){

  tau <- model$tau_post[length(model$tau_post)]
  k1  <- model$sampled_k1[length(model$sampled_k1)]

  model  <- model$final_trees

  var_mus <-  model$tree_data[[1]] %>%
    utils::tail(1) %>%
    tidyr::unnest(est_tree) %>%
    dplyr::count(group) %>%
    dplyr::mutate(var_mu = 1/((n + 1/k1)*tau)) %>%
    dplyr::select(-n)

  # Removing the intercept
  formula <- stats::as.formula(paste(c(formula), "- 1"))
  #formula <- stats::as.formula(paste(c(formula), "- 1"))

  response_name <- all.vars(formula)[1]
  group         <- dplyr::pull(newdata, !!group_variable)
  newdata$y     <- 0
  m             <- stats::model.frame(formula, data = newdata)
  X             <- stats::model.matrix(formula, m) %>%
    as.data.frame()

  names(X) <-  paste0("X", 1:ncol(X))
  #y <- newdata[ , response_name]
  group <- dplyr::pull(newdata, !!group_variable)

  newdata <- data.frame(X, group)

  all_results  <- model$results
  P <- length(all_results)

  newdata <- newdata %>%
    dplyr::mutate(id = 1:n())

  all_preds <- model %>%
    tibble::add_column(newdata = list(newdata)) %>%
    dplyr::mutate(final_tree = purrr::map(tree_data, ~utils::tail(.x, 1))) %>%
    #dplyr::rowwise() %>%
    dplyr::mutate(preds = purrr::pmap(list(results, final_tree, newdata, k1, tau), calculate_preds))  %>%
    dplyr::select(tree_index, preds)

  all_preds_wrangle <- all_preds %>%
    tidyr::unnest(cols = preds) %>%
    dplyr::select(tree_index, sampled_mu_j, group, id) %>%
    tidyr::pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>%
    tidyr::unnest(cols = dplyr::everything()) %>%
    dplyr::mutate(pred = rowSums(.[3:ncol(.)])) %>%
    dplyr::select(id, group, pred)

  all_preds_wrangle_final <- all_preds_wrangle %>%
    dplyr::left_join(newdata, by = c('id', 'group')) %>%
    dplyr::select(-id) %>%
    dplyr::left_join(var_mus, by = "group") %>%
    dplyr::mutate(var_mu = ifelse(is.na(var_mu), (1/tau), var_mu))

  return(all_preds_wrangle_final)
}

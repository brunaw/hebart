# Get rid of NOTES
globalVariables(
  c(".",            "%>%",       "d",            "depara_names",   "d_reduction",
    "cond",         "err",       "err_y",        "est_tree",       "final_tree",
    "group",        "id",        "k1",           "k2",             "mean_res",
    "mean_y_j",     "mu",        "mu_avg",       "mu_js_sampled",  "mu_samp",
    "muj",          "n",         "node",         "node_index",     "nodes_to_prune",
    "P",            "p_prune",   "p1_mu",        "p1_muj",         "p2_mu",
    "p2_muj",       "parent",    "pred",         "preds",          "prior_k1",
    "results",      "rn",        "sampled_mu_j", "sampled_mu",     "sum_errors",
    "sum_errors_y", "temp_node", "tree_data",    "tree_index",     "var",
    "var_mu",       "y",

    # New variables
    "Node Average",  "Tree Index",  "avg_mse",  "get_parent",
    "iter",  "mse",  "mses",  "n_nodes", "nodes",
    "tau",  "trees", "value", "id_iter", "id_obs", "id_tree", "se"

    ))

requireNamespace("ggplot2", quietly = TRUE)
requireNamespace("patchwork", quietly = TRUE)

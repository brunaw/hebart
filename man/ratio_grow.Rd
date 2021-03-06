% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grow_ratio.R
\name{ratio_grow}
\alias{ratio_grow}
\title{Final ratio for a growth step.}
\usage{
ratio_grow(tree, old_tree, current_node, pars, alpha_grow, beta_grow)
}
\arguments{
\item{tree}{The current tree}

\item{old_tree}{The old tree}

\item{current_node}{The current grown node}

\item{pars}{The full list of parameters}

\item{alpha_grow}{The alpha of the growing probability}

\item{beta_grow}{The beta of the growing probability}
}
\value{
The final ratio for the candidate tree
}
\description{
The final ratio is to be used as the acceptance
criteria in the MCMC of the hebart model
}
\author{
Bruna Wundervald, \email{brunadaviesw@gmail.com}.
}

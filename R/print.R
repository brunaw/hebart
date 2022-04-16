##' @title Print hebart
##' @param x Object of class 'hebart'
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{hebart}}
##' @author Bruna Wundervald
##' @export
print.hebart <- function(x, ...) {
  cat("Hebart result\n")
  cat("-----------------------------------\n")
  cat("Formula:\n", deparse(x$formula), "\n\n")
  cat("Number of trees:        ", x$P, "\n")
  cat("Number of covariates:   ", x$num.variables, "\n")
  cat("Prediction error (MSE): "   , x$mse, "\n")
  cat("R squared:              ",      x$r.squared, "\n")
}

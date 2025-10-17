#' Compute Residuals for HZIP Models
#'
#' This function calculates residuals for objects of class \code{HZIP}
#' using either Pearson residuals, randomized quantile residuals, or
#' adjusted quantile residuals. The computation is performed efficiently
#' using C++ functions for predicting random effects and calculating
#' residuals.
#'
#' @param object An object of class \code{HZIP} returned by the fitting function.
#' @param type A character string specifying the type of residuals to compute.
#'   Options are \code{"Pearson"}, \code{"quantile"}, or \code{"Adj.quantile"}.
#'   Default is \code{"quantile"}.
#' @param nodes Numeric vector of quadrature nodes for approximating integrals.
#' @param weights Numeric vector of quadrature weights corresponding to \code{nodes}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector of residuals with length equal to the total number
#'   of observations in the dataset.
#'
#' @details
#' The function internally groups the data by individual (\code{Ind}), constructs
#' model matrices for both zero-inflation and count parts of the model, and then
#' calls the C++ functions \code{predict_HZIP_cpp_vec} and \code{r_ij_cpp_vec}
#' to efficiently compute the residuals. Random effects are integrated using
#' adaptive quadrature based on the supplied \code{nodes} and \code{weights}.
#'
#' @examples
#' \dontrun{
#' # Assuming fit is an object of class HZIP
#' residuals(fit, type="Pearson", nodes=nodes_vec, weights=weights_vec)
#' residuals(fit, type="quantile", nodes=nodes_vec, weights=weights_vec)
#' }
#'
#' @importFrom stats model.frame model.matrix model.response qnorm
#' @importFrom dplyr group_split group_by
#' @export
residuals.HZIP <- function(object, type=c("Pearson","quantile","Adj.quantile"),
                           Q=15, ...) {
  type <- match.arg(type)

  formula <- object$formula
  data <- object$data
  theta1 <- object$coefficients_zero
  theta2 <- object$coefficients_count

  data_list <- dplyr::group_split(dplyr::group_by(data, Ind))
  xlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs=1))
  wlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs=2))
  ylist <- lapply(data_list, function(df) model.response(model.frame(Formula(formula), df)))

  QGauss <- statmod::gauss.quad(Q, kind = "hermite")
  nodes <- QGauss$nodes
  weights <- QGauss$weights

  vB <- predict_HZIP_cpp_vec(ylist, xlist, wlist, theta1, theta2, nodes, weights)

  res <- vector("list", length(ylist))
  for(i in seq_along(ylist))
    res[[i]] <- r_ij_cpp_vec(theta1, theta2, vB[i,,drop=FALSE], ylist[[i]],
                             xlist[[i]], wlist[[i]], type)

  do.call(c, res)
}

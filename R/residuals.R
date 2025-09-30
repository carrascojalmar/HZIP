#' @importFrom stats model.frame model.matrix model.response qnorm
#' @importFrom dplyr group_split group_by
#' @export
residuals.HZIP <- function(object,...) {

  formula <- object$formula
  data <-object$data
  coefficients_zero <- object$coefficients_zero
  coefficients_count <- object$coefficients_zero

  data_list <- dplyr::group_split(dplyr::group_by(data, Ind))

  xlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs = 1))
  wlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs = 2))
  ylist <- lapply(data_list, function(df) model.response(model.frame(Formula(formula), df)))

  Fq <- CDF(coefficients_zero,coefficients_count,xlist,wlist,ylist)

  #+runif(n)*PMF(coefficients_zero,coefficients_count,xlist,wlist,ylist)

  rq <- qnorm(Fq)

  return(rq)

}

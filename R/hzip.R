#' Fit a Hierarchical Zero-Inflated Poisson (HZIP) Model via GHQ
#'
#' \code{hzip()} fits a longitudinal/clustered zero-inflated Poisson model with
#' subject-level random effects by maximizing a (marginal) likelihood approximated
#' with Gauss–Hermite quadrature. The model uses a two-part \link[Formula]{Formula}:
#' \eqn{y ~ \text{count part} \mid \text{zero part}}, where the count intensity
#' (Poisson mean) and the zero-inflation probability are linked to (possibly
#' different) sets of covariates. Initial values are obtained from
#' \code{pscl::zeroinfl(..., dist = "poisson", link = "cloglog")}.
#'
#' @param formula A two-part \link[Formula]{Formula} of the form
#'   \code{y ~ x_count + ... | w_zero + ...}, where the right-hand side before
#'   the bar specifies covariates for the Poisson mean and the right-hand side
#'   after the bar specifies covariates for the zero-inflation component.
#' @param data A \code{data.frame} containing all variables used in \code{formula}
#'   and a subject identifier named \code{Ind} (one row per observation).
#' @param hessian Logical; if \code{TRUE} (default) the observed Hessian at the
#'   optimum is returned and used to compute standard-error estimates.
#' @param method Character string passed to \code{\link[stats]{optim}}
#'   (default \code{"BFGS"}).
#' @param Q Integer; number of Gauss–Hermite nodes for quadrature (default \code{15}).
#'   Larger values improve accuracy at higher computational cost.
#' @param control Optional \code{list} passed to \code{\link[stats]{optim}}'s
#'   \code{control=} argument (e.g., \code{list(maxit = 500)}).
#' @param ... Further arguments passed to \code{\link[stats]{optim}}.
#'
#' @details
#' Let \eqn{y_{ij}} denote the count response for subject \eqn{i} at occasion \eqn{j}.
#' The HZIP model assumes
#' \deqn{P(y_{ij}=0 \mid u_i) = \pi_{ij}(u_i) + \{1-\pi_{ij}(u_i)\}\exp\{-\mu_{ij}(u_i)\},}
#' \deqn{P(y_{ij}=k \mid u_i) = \{1-\pi_{ij}(u_i)\}\frac{\mu_{ij}(u_i)^k e^{-\mu_{ij}(u_i)}}{k!},\quad k\ge 1,}
#' with linear predictors for the count and zero parts (links typically \code{log}
#' for the Poisson mean and \code{cloglog} for the zero-inflation). Subject-specific
#' random effects \eqn{u_i} induce within-subject dependence; the marginal likelihood
#' is approximated by Gauss–Hermite quadrature with \code{Q} nodes.
#'
#' @return An object of class \code{"HZIP"}, a \code{list} with elements:
#' \item{call}{The matched call.}
#' \item{formula}{The model \code{Formula}.}
#' \item{coefficients_zero}{Estimated coefficients for the zero-inflation part.}
#' \item{coefficients_count}{Estimated coefficients for the count part.}
#' \item{scale_zero}{Estimated scale (zero part).}
#' \item{scale_count}{Estimated scale (count part).}
#' \item{loglik}{Optimized objective value returned by \code{optim}.
#'   (Note: depending on \code{lvero}, this may be the negative log-likelihood.)}
#' \item{convergence}{\code{optim} convergence code.}
#' \item{n}{Number of observations or subjects (see \strong{Note}).}
#' \item{m}{Cluster sizes per subject (vector ordered by \code{Ind}).}
#' \item{ep}{Approximate standard errors (square roots of the diagonal of the inverse Hessian).}
#' \item{iter}{Number of \code{optim} iterations.}
#' \item{method}{Optimization method.}
#' \item{optim}{Raw \code{optim} output.}
#' \item{data}{The input \code{data}.}
#'
#' @note
#' The subject identifier must be named \code{Ind}. The sign convention for the
#' zero-part coefficients in the initial values follows \code{pscl::zeroinfl};
#' the internal parameter vector is \code{c(scale_zero, -beta_zero, scale_count, beta_count)}.
#' Also verify whether \code{loglik} is the (negative) log-likelihood as returned by
#' your objective \code{lvero}; if it is the negative log-likelihood, you may want
#' to store \code{logLik = -op$value} for user convenience.
#'
#' @references
#' Min, Y., & Agresti, A. (2005). Random effect models for repeated measures of
#' zero-inflated count data. \emph{Statistical Modelling}, 5(1), 1–19.
#'
#' Jackman, S. (2020). \emph{pscl}: Classes and Methods for R Developed in the
#' Political Science Computational Laboratory. R package version 1.5.5.
#'
#' Zeileis, A., & Croissant, Y. (2010). Extended model formulas in R:
#' \emph{Journal of Statistical Software}, 34(1), 1–13. (\pkg{Formula})
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # toy panel structure
#' n_id <- 50; m_i <- 4
#' Ind  <- rep(1:n_id, each = m_i)
#' x1   <- rnorm(n_id * m_i)
#' w1   <- rbinom(n_id * m_i, 1, 0.5)
#' y    <- rpois(n_id * m_i, lambda = exp(0.3 + 0.5*x1))  # placeholder
#' dat  <- data.frame(Ind, y, x1, w1)
#'
#' # count part uses x1; zero part uses w1
#' fit <- hzip(y ~ x1 | w1, data = dat, Q = 15)
#' str(fit)
#' }
#'
#' @importFrom stats glm model.frame model.matrix model.response optim pnorm
#' @importFrom dplyr group_by group_split
#' @import Formula
#' @import pscl
#' @import statmod
#' @export
hzip <- function(formula, data, hessian = TRUE, method = "BFGS", Q = 15,
                 control=NULL,...) {


  data_list <- dplyr::group_split(dplyr::group_by(data, Ind))

  xlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs = 1))
  wlist <- lapply(data_list, function(df) model.matrix(Formula(formula), df, rhs = 2))
  ylist <- lapply(data_list, function(df) model.response(model.frame(formula, df)))

  lhs <- formula(Formula(formula), lhs = 1, rhs = 0)
  rhs1 <- formula(Formula(formula), lhs = 0, rhs = 1)
  rhs2 <- formula(Formula(formula), lhs = 0, rhs = 2)

  fAux <- paste(deparse(lhs[[2]]), "~",
                      deparse(rhs2[[2]]), "|",
                      deparse(rhs1[[2]]))

  fit.aux <- zeroinfl(Formula(as.formula(fAux)),data=data, dist = "poisson",
                      link="cloglog")


  initial <- c(1, -as.numeric(fit.aux$coefficients$zero), 1,
               as.numeric(fit.aux$coefficients$count))

  QGauss <- statmod::gauss.quad(Q, kind = "hermite")
  Qnodes <- QGauss$nodes
  Qweights <- QGauss$weights

  op <- optim(
    par = initial,
    fn = lvero,
    ylist = ylist,
    xlist = xlist,
    wlist = wlist,
    Qnodes = Qnodes,
    Qweights = Qweights,
    method = method,
    hessian = hessian,
    control = list(),
    ...
  )

  p1 <- dim(xlist[[1]])[2]
  p2 <- dim(wlist[[1]])[2]


  fit.hzip <- list(
    call = match.call(),
    formula = formula,
    coefficients_zero = op$par[2:(p1+1)],
    coefficients_count = op$par[(p1+3):(p1+p2+2)],
    scale_zero = op$par[1],
    scale_count= op$par[(p1+2)],
    loglik = op$value,
    convergence =op$convergence,
    n = nrow(data),
    n = max(data$Ind),
    m = as.numeric(table(data$Ind)),
    ep = sqrt(diag(solve(op$hessian))),
    iter = op$counts[1],
    method = method,
    optim = op,
    data = data
  )

  class(fit.hzip) <- "HZIP"
  return(fit.hzip)
}

#' @export
print.MCR <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients (Zero part):\n")
  print(round(x$coefficients_zero, 4))

  cat("\nCoefficients (Count part):\n")
  print(round(x$coefficients_count, 4))

  cat("\nZero scale:\n")
  print(round(x$scale_zero, 4))

  cat("\nCount scale:\n")
  print(round(x$scale_count, 4))

  cat("\nLog-likelihood:", round(x$loglik, 4), "\n")
}

#' @export
summary.HZIP <- function(object, ...) {

  y <- model.response(model.frame(Formula(object$formula), data = eval(object$call$data)))
  x <- model.matrix(Formula(object$formula), data = eval(object$call$data), rhs = 1)
  w <- model.matrix(Formula(object$formula), data = eval(object$call$data), rhs = 2)

  coef_zero <- object$coefficients_zero
  coef_count <- object$coefficients_count
  scale_zero <- object$scale_zero
  scale_count <- object$scale_count
  ep <- object$ep
  iter <- object$iter
  value <- object$loglik
  n <- object$n
  convergence <- object$convergence

  p1 <- ncol(x)
  p2 <- ncol(w)

  std_scale_zero <- ep[1]
  std_coeff_zero <- ep[2:(p1+1)]

  std_scale_count <- ep[(p1+2)]
  std_coeff_count <- ep[(p1+3):(p1+p2+2)]


  z_beta_zero <- coef_zero / std_coeff_zero
  p_beta_zero <- 2 * (1 - pnorm(abs(z_beta_zero)))

  z_eta_count <- coef_count / std_coeff_count
  p_eta_count <- 2 * (1 - pnorm(abs(z_eta_count)))

  df_coef_zero <- data.frame(
    Value = coef_zero,
    `Std. Error` = std_coeff_zero,
    z = z_beta_zero,
    p = format.pval(p_beta_zero, digits = 2, eps = 2e-16),
    row.names = colnames(x),
    check.names = FALSE
  )

  df_coef_count <- data.frame(
    Value = coef_count,
    `Std. Error` = std_coeff_count,
    z = z_eta_count,
    p = format.pval(p_eta_count, digits = 2 , eps = 2e-16),
    row.names = colnames(w),
    check.names = FALSE
  )

  df_coef_scale_zero <- data.frame(
      Estimate = scale_zero,
      Std.Error = std_scale_zero,
      row.names = "scale"
    )

  df_coef_scale_count <- data.frame(
    Estimate = scale_count,
    Std.Error = std_scale_count,
    row.names = "scale"
  )

  out <- list(
    call = object$call,
    loglik = value,
    AIC = -2 * value + 2 * length(ep),
    BIC = -2 * value + log(n) * length(ep),
    iter = iter,
    coef_zero = df_coef_zero,
    coef_count = df_coef_count,
    scale_zero = df_coef_scale_zero,
    scale_count = df_coef_scale_count,
    convergence = convergence
  )

  class(out) <- "summary.HZIP"
  return(out)
}

#' @export
print.summary.HZIP <- function(x, digits = 5, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Log-Likelihood:", formatC(x$loglik, digits = digits, format = "f"), "\n")
  cat("AIC:", formatC(x$AIC, digits = digits, format = "f"), "\n")
  cat("BIC:", formatC(x$BIC, digits = digits, format = "f"), "\n")
  cat("Number of Iterations:", x$iter, "\n")
  cat("Convergence:", x$convergence, "\n\n")

  cat("Coefficients (Zero part):\n")
  print(format(x$coef_zero, digits = digits, nsmall = digits), quote = FALSE)

  cat("\nCoefficients (Count part):\n")
  print(format(x$coef_count, digits = digits, nsmall = digits), quote = FALSE)

  cat("\nZero Scale:\n")
  print(format(x$scale_zero, digits = digits, nsmall = digits), quote = FALSE)

  cat("\nCount Scale:\n")
  print(format(x$scale_count, digits = digits, nsmall = digits), quote = FALSE)

}

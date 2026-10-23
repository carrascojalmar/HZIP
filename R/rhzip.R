# rgengamma (interno)
#' @keywords internal
#' @noRd
#' @importFrom stats rgamma rnorm
rgengamma <- function(n, mu=0, sigma=1, lambda, zero=0.0001) {
  # generates n random generalized loggamma variates
  if (abs(lambda) > zero) {
    alpha <- 1/lambda^2
    ei    <- rgamma(n,shape=alpha,rate=1)
    ui    <- (log(ei)-log(alpha))/lambda
  } else {
    ui <- rnorm(n,mean=0,sd=1)
  }
  ui <- ui*sigma + mu
  return(ui)
}

#' Simulate Data from a Hierarchical Zero-Inflated Poisson (HZIP) Model
#'
#' \code{rhzip()} generates panel/longitudinal data from a two–part
#' hierarchical zero-inflated Poisson model with subject-specific random
#' effects for both the zero-inflation and the count components. Random
#' effects are drawn from a generalized log-gamma (GLG) distribution via
#' \code{rgengamma()} (user must provide/attach a function with this name; see
#' \strong{Dependencies}).
#'
#' @param n Integer. Number of subjects.
#' @param m Integer vector of length \code{n} (or a scalar recycled to length
#'   \code{n}). Numbers of repeated measurements per subject.
#' @param para Numeric vector of parameters in the order
#'   \code{c(lambda1, beta1, lambda2, beta2)}, where:
#'   \itemize{
#'     \item \code{lambda1}: GLG scale/shape parameter for the zero part.
#'     \item \code{beta1}:  length-\code{p1} vector of coefficients for the zero part
#'           (matches \code{ncol(x1)}).
#'     \item \code{lambda2}: GLG scale/shape parameter for the count part.
#'     \item \code{beta2}:  length-\code{p2} vector of coefficients for the count part
#'           (matches \code{ncol(x2)}).
#'   }
#'   Internally, the linear predictors are
#'   \eqn{\eta^{(0)}_{ij} = x^{\top}_{1,ij}\beta_1 + b^{(0)}_i} and
#'   \eqn{\eta^{(1)}_{ij} = x^{\top}_{2,ij}\beta_2 + b^{(1)}_i}, with
#'   \eqn{p_{ij} = 1 - \exp\{-\exp(\eta^{(0)}_{ij})\}} (cloglog link for zero-inflation)
#'   and \eqn{\mu_{ij} = \exp(\eta^{(1)}_{ij})} (log link for counts).
#' @param x1 Numeric matrix of covariates for the zero-inflation part
#'   (dimension \code{sum(m)} by \code{p1}). Include an intercept column if
#'   desired; every column supplied is kept in the output (see
#'   \strong{Column naming}).
#' @param x2 Numeric matrix of covariates for the count part
#'   (dimension \code{sum(m)} by \code{p2}). Include an intercept column if
#'   desired; every column supplied is kept in the output (see
#'   \strong{Column naming}).
#'
#' @details
#' For subject \eqn{i=1,\dots,n} with \eqn{m_i} observations, the model draws
#' two subject-level random effects \eqn{b^{(0)}_i} and \eqn{b^{(1)}_i} independently
#' from GLG distributions parameterized by \code{lambda1} and \code{lambda2}.
#' Conditional on these effects, outcomes are generated as
#' \deqn{Y_{ij} = Z_{ij}\times C_{ij},}
#' where \eqn{Z_{ij}\sim \mathrm{Bernoulli}(p_{ij})} controls structural zeros and
#' \eqn{C_{ij}\sim \mathrm{Poisson}(\mu_{ij})} controls the count size.
#'
#' The returned data frame contains subject IDs, the response \code{y}, and the
#' supplied covariates. An attribute \code{"propZeros"} stores a small summary
#' with the percentage of structural zeros, additional Poisson zeros, and total
#' zeros.
#'
#' @return A \code{tibble} with columns:
#' \item{Ind}{Subject identifier (1..n), repeated according to \code{m}.}
#' \item{y}{Simulated response.}
#' \item{x*, w*}{The covariates from \code{x1} and \code{x2} (renamed as described below).}
#'
#' The object has an attribute \code{"propZeros"}: a 3×1 \code{data.frame} with
#' rows \code{"Zeros"} (structural zeros, \%),
#' \code{"Count"} (extra zeros from the Poisson part, \%), and
#' \code{"Total"} (overall zero percentage).
#'
#' @section Column naming:
#' All columns of \code{x1} are kept and renamed \code{x1, x2, ..., xp1}
#' (zero-inflation part). All columns of \code{x2} are kept and renamed
#' \code{w1, w2, ..., wp2} (count part). No column is dropped: if \code{x1}
#' and/or \code{x2} include an intercept column (a column of 1's), it is
#' preserved in the output (as \code{x1} and/or \code{w1}, respectively) and
#' can be referenced explicitly in the \code{\link{hzip}} formula (using
#' \code{-1} on that side to avoid R adding a duplicate automatic intercept),
#' or simply omitted from the formula (letting R's automatic intercept play
#' the same role).
#'
#' @section Dependencies:
#' This function calls \code{rgengamma()} to draw GLG random effects. Ensure that
#' such a function is available on the search path (e.g., from a package that
#' provides a generalized log-gamma RNG) or provide your own implementation
#' with the signature \code{rgengamma(n, mu, sigma, lambda)}. It also uses
#' \pkg{dplyr} and \pkg{tibble}.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 50
#' m <- rep(4, n)
#' N <- sum(m)
#'
#' # design matrices (with intercepts)
#' x1 <- cbind(1, rnorm(N))
#' x2 <- cbind(1, rnorm(N), rbinom(N, 1, 0.5))
#'
#' p1 <- ncol(x1); p2 <- ncol(x2)
#' lambda1 <- 0.7
#' beta1   <- c(-0.2, 0.6)
#' lambda2 <- 0.9
#' beta2   <- c( 0.3, 0.5, -0.4)
#'
#' para <- c(lambda1, beta1, lambda2, beta2)
#'
#' sim <- rhzip(n, m, para, x1, x2)
#' head(sim)
#' attr(sim, "propZeros")
#'
#' # Fit it back with hzip(). With the new convention NO column is dropped,
#' # so x1 (2 cols: intercept, slope) comes out as "x1,x2" and x2 (3 cols:
#' # intercept, slope, slope) comes out as "w1,w2,w3". Referencing only the
#' # real slopes ("x2" and "w2+w3") and letting R's formula add the
#' # intercept automatically reproduces the true model exactly (equivalently,
#' # one could write "x1+x2-1 | w1+w2+w3-1" using the explicit intercept
#' # columns and disabling the automatic one).
#' \donttest{
#' fit <- hzip(y ~ x2 | w2 + w3, data = sim)
#' summary(fit)
#' }
#'
#' # Example WITHOUT an intercept in the zero part (e.g. simulation Scenario 1
#' # style): x1 has no column of 1's, so every column of x1 is a real slope
#' # and must be fit with "-1" to avoid an automatic (spurious) intercept.
#' \donttest{
#' n2 <- 200; m2 <- rep(4, n2); N2 <- sum(m2)
#' w1b <- rep(rnorm(n2), times = m2)
#' w2b <- rep(runif(n2), times = m2)
#' x1b <- cbind(w1b, w2b)             # zero part: NO intercept
#' x2b <- cbind(1, w1b, w2b)          # count part: WITH intercept
#' parab <- c(0.8, 1, -1, 0.5, 1, 1, 1)   # (lambda1,beta11,beta12,lambda2,beta20,beta21,beta22)
#' simb <- rhzip(n2, m2, parab, x1b, x2b)
#' head(simb)   # columns: Ind, y, x1, x2 (zero part, both real slopes), w1, w2, w3 (count part)
#' fitb <- hzip(y ~ x1 + x2 - 1 | w2 + w3, data = simb)
#' summary(fitb)
#' }
#'
#' @seealso \code{\link{hzip}} for model fitting.
#'
# rhzip
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' @importFrom stats rbinom rpois
#' @export

rhzip <- function(n,m,para,x1,x2){
  p1 <- ncol(x1)
  p2 <- ncol(x2)

  lambda1 <- para[1]
  beta1 <- para[2:(p1+1)]

  lambda2 <- para[(p1+2)]
  beta2 <- para[(p1+3):(p1+p2+2)]

  ind <- rep(seq_len(n), times = m)

  sample1 <- data.frame(x1,ind)
  sample2 <- data.frame(x2,ind)

  y <- vector()
  conta1 <- vector()
  conta2 <- vector()
  i <- 1
  while(i <=n){
    s.aux1 <- filter(sample1,ind==i)
    b1 <-  rgengamma(n = 1, mu = 0, sigma = lambda1, lambda =lambda1)
    eta1 <- as.matrix(s.aux1[,1:ncol(x1)])%*%beta1+b1
    p <- 1-exp(-exp(eta1))

    s.aux2 <- filter(sample2,ind==i)
    b2 <- rgengamma(n = 1, mu = 0, sigma = lambda2, lambda =lambda2)
    eta2 <- as.matrix(s.aux2[,1:ncol(x2)])%*%beta2+b2
    mu <- exp(eta2)

    # y.auxi1 = 1 indica processo de contagem ativo (nao estrutural), prob = 1-p
    # (p = pi_ij = P(zero estrutural) via link cloglog; Y=0 ocorre se y.auxi1=0 OU y.auxi2=0)
    y.auxi1 <- rbinom(m[i],1,1-p)
    y.auxi2 <- rpois(m[i],mu)

    conta1 <- c(conta1,y.auxi1)
    conta2 <- c(conta2,y.auxi2)

    y <- c(y,y.auxi1*y.auxi2) #rZIP(mi,mu=mu,sigma=p)

    i <- i+1
  }

  colnames(x1) <- paste0("x", seq_len(ncol(x1)))
  colnames(x2) <- paste0("w", seq_len(ncol(x2)))

  out <- tibble(
    Ind = ind,
    y = y,
    !!!as.data.frame(x1),
    !!!as.data.frame(x2)
  )

  # Proportions

  pConta1 <- 100*(1-mean(conta1))
  pConta2 <- 100*mean(ifelse(conta2==0,1,0))
  pConta3 <- 100*mean(ifelse(y==0,1,0))

  propZeros <- data.frame(rbind(pConta1,(pConta3-pConta1),pConta3))
  colnames(propZeros) <- "Zero-proportions(%)"
  rownames(propZeros) <- c("Zeros","Count","Total")
  propZeros

  attr(out, "propZeros")  <- propZeros

  return(out)
}

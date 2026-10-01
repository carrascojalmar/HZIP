#' Salamanders data
#'
#' This dataset is adapted from the \pkg{glmmTMB} package and contains
#' salamander counts with information on mining status and species.
#' It is intended for illustrating zero-inflated Poisson models with
#' random effects using the \code{hzip()} function.
#'
#' @format A data frame with 644 rows and 5 variables:
#' \describe{
#'   \item{Ind}{Cluster identifier (integer, 1 to 161); each cluster
#'     groups the 4 repeated samples of one site-species combination.}
#'   \item{y}{Count response variable (integer).}
#'   \item{mined}{Mining status: \code{"yes"} or \code{"no"}.}
#'   \item{spp}{Species factor with multiple levels (e.g., GP, PR, DM, etc.).}
#'   \item{sample}{Repeated sampling occasion (factor with levels
#'     \code{"1"} to \code{"4"}).}
#' }
#'
#' @details
#' The dataset was originally included in the \pkg{glmmTMB} package
#' (Brooks et al., 2017), and has been slightly modified for testing
#' the \pkg{HZIP} package: observations are ordered by site, the
#' variables \code{count}, \code{mined}, \code{spp} and \code{sample}
#' are kept, \code{count} is renamed to \code{y}, and the cluster
#' identifier \code{Ind} is added.
#'
#' @source
#' Adapted from the \pkg{glmmTMB} package.
#'
#' @examples
#' \donttest{
#' data(salamanders, package = "HZIP")
#'
#' ## Fit zero-inflated Poisson with random effects
#' fit.salamander <- hzip(y ~ mined+spp+mined:spp | mined+spp+sample,
#'                        data = salamanders)
#' summary(fit.salamander)
#' }
"salamanders"

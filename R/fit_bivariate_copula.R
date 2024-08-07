#' Fit Bivariate Copula Model with Optional Marginal Parameter Estimation
#'
#' @description
#' The `fit_bivariate_copula` function fits a bivariate copula model to data with optional marginal parameter estimation.
#' It supports Gaussian, Clayton, and Joe copulas, and allows for specifying normal, lognormal, exponential, or beta marginal distributions for each variable.
#' This function utilizes Stan for Bayesian inference.
#'
#' @usage
#' fit_bivariate_copula(
#'   U,
#'   copula,
#'   marginals,
#'   iter = 2000,
#'   chains = 4,
#'   warmup = 1000,
#'   thin = 1,
#'   seed = NULL,
#'   control = list(adapt_delta = 0.8, max_treedepth = 10),
#'   cores = 1
#' )
#'
#' @param U Data matrix of observed marginals, with two variables.
#' @param copula The type of copula to fit. Options are "gaussian", "clayton", and "joe".
#' @param marginals A list specifying the marginal distributions for each variable. Options are "normal", "lognormal", "exponential", and "beta".
#' @param iter Number of iterations for each chain. Default is 2000.
#' @param chains Number of chains. Default is 4.
#' @param warmup Number of warmup iterations per chain. Default is 1000.
#' @param thin Thinning rate. Default is 1.
#' @param seed Random seed. Default is NULL.
#' @param control A list of parameters to control the sampler's behavior. Default is list(adapt_delta = 0.8, max_treedepth = 10).
#' @param cores Number of cores to use for parallel processing. Default is 1.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\code{fit}}{Stan model object.}
#' }
#'
#' @examples
#' \dontrun{
#' library(copulaStan)
#' library(copula)
#'
#' seed <- 2024
#' set.seed(seed)
#' true_rho <- 0.5
#' n <- 2000
#'
#' margins <- c("norm", "lnorm")
#' params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
#' cop <- normalCopula(param = true_rho, dim = 2)
#' mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
#' data <- rMvdc(n, mvdc_copula)
#'
#' fit <- fit_bivariate_copula(data,
#'                             copula = "gaussian", marginals = c("normal", "lognormal"),
#'                             seed = seed)
#' }
#'
#' @export
fit_bivariate_copula <- function(U, copula, marginals,
                                 iter = 2000, chains = 4, warmup = 1000, thin = 1,
                                 seed = NULL,
                                 control = list(adapt_delta = 0.8, max_treedepth = 10),
                                 cores = 1) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("'rstan' is required but is not installed.")
  }

  if (ncol(U) != 2) {
    stop("matrix U must have exactly two columns.")
  }

  N <- nrow(U)
  y1 <- U[, 1]
  y2 <- U[, 2]

  dist_map <- c("normal" = 1, "lognormal" = 2, "exponential" = 3, "beta" = 4)
  dist1 <- dist_map[marginals[[1]]]
  dist2 <- dist_map[marginals[[2]]]

  copula_map <- c("gaussian" = 1, "clayton" = 2, "joe" = 3)
  copula_type <- copula_map[copula]

  if (is.null(dist1) || is.null(dist2)) {
    stop("Invalid marginal distribution specified. Use 'normal', 'lognormal', 'exponential', or 'beta'.")
  }

  if (is.null(copula_type)) {
    stop("Invalid copula specified. Use 'gaussian', 'clayton', or 'joe'.")
  }

  stan_data <- list(N = N, y1 = y1, y2 = y2, dist1 = dist1, dist2 = dist2, copula_type = copula_type)

  stan_file <- system.file("stan", "fit_bivariate_copula.stan", package = "copulaStan")

  if (stan_file == "") {
    stop("Stan model file not found in the package.")
  }

  init_function <- function() {
    list(
      mu1 = if (dist1 == 1 || dist1 == 2) array(0, 1) else numeric(0),
      sigma1 = if (dist1 == 1 || dist1 == 2) array(1, 1) else numeric(0),
      lambda1 = if (dist1 == 3) array(1, 1) else numeric(0),
      alpha1 = if (dist1 == 4) array(1, 1) else numeric(0),
      beta1 = if (dist1 == 4) array(1, 1) else numeric(0),
      mu2 = if (dist2 == 1 || dist2 == 2) array(0, 1) else numeric(0),
      sigma2 = if (dist2 == 1 || dist2 == 2) array(1, 1) else numeric(0),
      lambda2 = if (dist2 == 3) array(1, 1) else numeric(0),
      alpha2 = if (dist2 == 4) array(1, 1) else numeric(0),
      beta2 = if (dist2 == 4) array(1, 1) else numeric(0),
      rho = if (copula_type == 1) array(0, 1) else numeric(0),
      theta = if (copula_type != 1) array(1, 1) else numeric(0)
    )
  }

  options(mc.cores = cores)

  fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    chains = chains,
    warmup = warmup,
    thin = thin,
    seed = seed,
    control = control,
    init = init_function
  )

  return(fit)
}

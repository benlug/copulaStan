# Fit Bivariate Copula Model

Fits a bivariate copula model to data with marginal parameter
estimation. Supports Gaussian, Clayton, and Joe copulas with normal,
lognormal, exponential, or beta marginal distributions. Uses CmdStan for
Bayesian inference via the `cmdstanr` package.

## Usage

``` r
fit_bivariate_copula(
  U,
  copula,
  marginals,
  iter = 1000,
  chains = 4,
  warmup = 1000,
  thin = 1,
  seed = NULL,
  adapt_delta = 0.8,
  max_treedepth = 10,
  parallel_chains = 1,
  refresh = 500
)
```

## Arguments

- U:

  A numeric matrix with exactly two columns containing the observed
  data. Each column corresponds to one variable.

- copula:

  Character string specifying the copula type. One of `"gaussian"`,
  `"clayton"`, or `"joe"`.

- marginals:

  A character vector of length 2 specifying the marginal distributions.
  Each element must be one of `"normal"`, `"lognormal"`,
  `"exponential"`, or `"beta"`.

- iter:

  Number of sampling iterations per chain (after warmup). Default is
  1000.

- chains:

  Number of MCMC chains. Default is 4.

- warmup:

  Number of warmup iterations per chain. Default is 1000.

- thin:

  Thinning rate. Default is 1.

- seed:

  Random seed for reproducibility. Default is `NULL`.

- adapt_delta:

  Target acceptance rate for NUTS. Default is 0.8.

- max_treedepth:

  Maximum tree depth for NUTS. Default is 10.

- parallel_chains:

  Number of chains to run in parallel. Default is 1.

- refresh:

  How often to print progress (in iterations). Set to 0 for silent.
  Default is 500.

## Value

A `copula_fit` object (S3 class) containing:

- `fit`:

  The CmdStanMCMC fit object from `cmdstanr`.

- `copula`:

  The copula type used.

- `marginals`:

  The marginal distributions used.

- `data_dim`:

  Dimensions of the input data (rows, columns).

## Examples

``` r
if (FALSE) { # \dontrun{
library(copula)
library(copulaStan)

set.seed(2024)
n <- 1000
cop <- normalCopula(param = 0.5, dim = 2)
margins <- c("norm", "lnorm")
params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "gaussian",
  marginals = c("normal", "lognormal"),
  seed = 2024
)
print(fit)
summary(fit)
} # }
```

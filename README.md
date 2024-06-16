# copulaStan

The `copulaStan` package provides functions to fit bivariate Gaussian and Clayton copula models using Stan. These models can handle both normal and non-normal marginals. 

## Installation

You can install the development version of `copulaStan` from [GitHub](https://github.com/benlug/copulaStan) with:

```r
# install.packages("devtools")
devtools::install_github("benlug/copulaStan")
```

## Example

This is a basic example which shows you how to fit a bivariate Gaussian copula model:

```r
library(copula)
library(copulaStan)

set.seed(123)
true_rho <- 0.5
cop <- normalCopula(param = true_rho, dim = 2)
U <- rCopula(1000, cop)

fit <- fit_bivariate_copula(U, copula = "gaussian", marginals = c("normal", "normal"), seed = 123)
fit
```

And an example for fitting a bivariate Clayton copula model:

```r
library(copula)
library(copulaStan)

seed <- 2024
set.seed(seed)
true_rho <- 0.5
n <- 2000  

margins <- c("norm", "lnorm")
params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
cop <- normalCopula(param = true_rho, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data, copula = "gaussian", marginals = margins,
                            seed = seed)
```

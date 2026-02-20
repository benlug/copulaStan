# copulaStan: Fitting Bivariate Copula Models

## Introduction

This vignette demonstrates how to fit bivariate copula models using the
`copulaStan` package. The package supports Gaussian, Clayton, and Joe
copulas with normal, lognormal, exponential, and beta marginal
distributions.

``` r
library(copulaStan)
library(copula)
```

### Example: Gaussian Copula

#### Normal + Lognormal Marginals

``` r
seed <- 2024
set.seed(seed)
true_rho <- 0.5
n <- 2000

margins <- c("norm", "lnorm")
params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
cop <- normalCopula(param = true_rho, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "gaussian", marginals = c("normal", "lognormal"),
  seed = seed
)
fit
summary(fit)
coef(fit)
```

#### Normal + Exponential Marginals

``` r
seed <- 2024
set.seed(seed)
true_rho <- 0.5
n <- 2000

margins <- c("norm", "exp")
params <- list(list(mean = 0.8, sd = 2), list(rate = 1))
cop <- normalCopula(param = true_rho, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "gaussian", marginals = c("normal", "exponential"),
  seed = seed
)
fit
```

#### Beta + Beta Marginals

``` r
seed <- 2024
set.seed(seed)
true_rho <- 0.5
n <- 2000

margins <- c("beta", "beta")
params <- list(list(shape1 = 2, shape2 = 5), list(shape1 = 3, shape2 = 4))
cop <- normalCopula(param = true_rho, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "gaussian", marginals = c("beta", "beta"),
  seed = seed
)
fit
```

### Example: Clayton Copula

#### Normal + Lognormal Marginals

``` r
seed <- 2024
set.seed(seed)
true_theta <- 2.0
n <- 2000

margins <- c("norm", "lnorm")
params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
cop <- claytonCopula(param = true_theta, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "clayton", marginals = c("normal", "lognormal"),
  seed = seed
)
fit
```

### Example: Joe Copula

#### Normal + Lognormal Marginals

``` r
seed <- 2024
set.seed(seed)
true_theta <- 2.0
n <- 2000

margins <- c("norm", "lnorm")
params <- list(list(mean = 0.8, sd = 2), list(meanlog = 0, sdlog = 0.8))
cop <- joeCopula(param = true_theta, dim = 2)
mvdc_copula <- mvdc(cop, margins = margins, paramMargins = params)
data <- rMvdc(n, mvdc_copula)

fit <- fit_bivariate_copula(data,
  copula = "joe", marginals = c("normal", "lognormal"),
  seed = seed
)
fit
```

## Mathematical Derivations

### Gaussian Copula Log-Likelihood

**Transformation to Uniform Marginals**:
$$U_{1} = F_{1}\left( X_{1} \right),\quad U_{2} = F_{2}\left( X_{2} \right)$$
where $F_{1}$ and $F_{2}$ are the CDFs of the marginals.

**Transformation to Standard Normal**:
$$Z_{1} = \Phi^{- 1}\left( U_{1} \right),\quad Z_{2} = \Phi^{- 1}\left( U_{2} \right)$$

**Copula Density** (ratio of bivariate to product of univariate
normals):
$$c\left( u_{1},u_{2};\rho \right) = \frac{1}{\sqrt{1 - \rho^{2}}}\exp\left( - \frac{1}{2\left( 1 - \rho^{2} \right)}\left( z_{1}^{2} + z_{2}^{2} - 2\rho z_{1}z_{2} \right) + \frac{z_{1}^{2}}{2} + \frac{z_{2}^{2}}{2} \right)$$

**Log-Likelihood**:
$$\log c\left( u_{1},u_{2};\rho \right) = - \frac{1}{2}\log\left( 1 - \rho^{2} \right) + \frac{\rho z_{1}z_{2} - \frac{1}{2}\rho^{2}\left( z_{1}^{2} + z_{2}^{2} \right)}{1 - \rho^{2}}$$

### Clayton Copula Log-Likelihood

The Clayton copula CDF:
$$C(u,v;\theta) = \left( u^{- \theta} + v^{- \theta} - 1 \right)^{- 1/\theta},\quad\theta > 0$$

The density (second mixed partial derivative):
$$c(u,v;\theta) = (\theta + 1)\, u^{- \theta - 1}\, v^{- \theta - 1}\left( u^{- \theta} + v^{- \theta} - 1 \right)^{- 1/\theta - 2}$$

Log-likelihood:
$$\log c(u,v;\theta) = \log(1 + \theta) + ( - \theta - 1)\left( \log u + \log v \right) - \frac{2\theta + 1}{\theta}\log\left( u^{- \theta} + v^{- \theta} - 1 \right)$$

### Joe Copula Log-Likelihood

The Joe copula CDF:
$$C(u,v;\theta) = 1 - \left\lbrack {\bar{u}}^{\theta} + {\bar{v}}^{\theta} - {\bar{u}}^{\theta}{\bar{v}}^{\theta} \right\rbrack^{1/\theta},\quad\theta \geq 1$$
where $\bar{u} = 1 - u$, $\bar{v} = 1 - v$.

Let $a = {\bar{u}}^{\theta}$, $b = {\bar{v}}^{\theta}$,
$S = a + b - ab$. The density:
$$c(u,v;\theta) = \theta\,{\bar{u}}^{\theta - 1}\,{\bar{v}}^{\theta - 1}\, S^{1/\theta - 2}\left\lbrack S + \frac{\theta - 1}{\theta}(1 - a)(1 - b) \right\rbrack$$

Log-likelihood:
$$\log c(u,v;\theta) = \log\theta + (\theta - 1)\left\lbrack \log\bar{u} + \log\bar{v} \right\rbrack + \left( \frac{1}{\theta} - 2 \right)\log S + \log\left\lbrack S + \frac{\theta - 1}{\theta}(1 - a)(1 - b) \right\rbrack$$

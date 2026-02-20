# Package index

## Fitting Models

Main function for fitting bivariate copula models.

- [`fit_bivariate_copula()`](https://benlug.github.io/copulaStan/reference/fit_bivariate_copula.md)
  : Fit Bivariate Copula Model

## Working with Results

S3 methods for inspecting and extracting results from fitted models.

- [`print(`*`<copula_fit>`*`)`](https://benlug.github.io/copulaStan/reference/print.copula_fit.md)
  : Print a copula_fit object
- [`summary(`*`<copula_fit>`*`)`](https://benlug.github.io/copulaStan/reference/summary.copula_fit.md)
  : Summarize a copula_fit object
- [`coef(`*`<copula_fit>`*`)`](https://benlug.github.io/copulaStan/reference/coef.copula_fit.md)
  : Extract point estimates from a copula_fit

## Internals

Internal constructor and utilities (not intended for direct use).

- [`copula_fit()`](https://benlug.github.io/copulaStan/reference/copula_fit.md)
  : Create a copula_fit object
- [`get_stan_model()`](https://benlug.github.io/copulaStan/reference/get_stan_model.md)
  : Get the compiled Stan model (compiles on first use, caches for
  reuse)

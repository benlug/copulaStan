# Summarize a copula_fit object

Summarize a copula_fit object

## Usage

``` r
# S3 method for class 'copula_fit'
summary(object, ...)
```

## Arguments

- object:

  A `copula_fit` object.

- ...:

  Additional arguments (unused).

## Value

A tibble of parameter summaries from
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

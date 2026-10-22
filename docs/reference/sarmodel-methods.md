# Tidy and glance methods for SITSAR models

Tidy and glance methods for SITSAR models

## Usage

``` r
# S3 method for class 'SARModel'
tidy(x, ...)

# S3 method for class 'SARModel'
glance(x, ...)

# S3 method for class 'SARModel'
logLik(object, ...)
```

## Arguments

- x:

  A fitted
  [SARModel](https://jodavid.github.io/SITSAR/reference/SARModel.md)
  object.

- ...:

  Not used.

- object:

  A fitted
  [SARModel](https://jodavid.github.io/SITSAR/reference/SARModel.md)
  object.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
tibble with one row per scalar parameter;
[`glance()`](https://generics.r-lib.org/reference/glance.html) a one-row
tibble; [`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of
class `logLik`.

## Examples

``` r
set.seed(1)
m <- G0I$new()$fit(rg0i(300, -2, 1, 4), L = 4)
generics::tidy(m)
#> # A tibble: 3 × 4
#>   term  estimate std.error fixed
#>   <chr>    <dbl>     <dbl> <lgl>
#> 1 alpha   -1.69      0.173 FALSE
#> 2 gamma    0.808     0.109 FALSE
#> 3 L        4        NA     TRUE 
generics::glance(m)
#> # A tibble: 1 × 7
#>   model     n logLik   AIC   BIC entropy converged
#>   <chr> <int>  <dbl> <dbl> <dbl>   <dbl> <lgl>    
#> 1 G0I     300  -274.  553.  560.   0.897 TRUE     
logLik(m)
#> 'log Lik.' -274.4602 (df=2)
```

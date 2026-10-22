# Estimation of G0I parameters

Maximum likelihood estimation of the \\\mathcal{G}\_I^0\\ parameters.
This is a convenience wrapper around
[G0I](https://jodavid.github.io/SITSAR/reference/G0I.md)`$fit()`, kept
for compatibility with previous versions of the package.

## Usage

``` r
g0iestimation(Theta, x, L_fixed = FALSE)
```

## Arguments

- Theta:

  Initial values `c(alpha, gamma, L)`. When `L_fixed = TRUE`, `Theta[3]`
  is the known number of looks.

- x:

  Vector of observations.

- L_fixed:

  Logical; if `TRUE`, \\L\\ is held fixed at `Theta[3]`.

## Value

Named vector `c(alpha, gamma, L)` of estimates.

## Examples

``` r
set.seed(1)
y <- rg0i(1000, alpha = -1.5, gamma = 2, L = 8)
g0iestimation(c(-1, 1, 1), y)
#>     alpha     gamma         L 
#> -1.500017  2.021340 10.188023 
g0iestimation(c(-1, 1, 8), y, L_fixed = TRUE)
#>     alpha     gamma         L 
#> -1.547246  2.139702  8.000000 
```

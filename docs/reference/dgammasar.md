# The Gamma distribution for SAR intensities

Density, distribution function, quantile function and random generation
for the multilook intensity Gamma law parametrised by its mean \\\mu\\
and number of looks \\L\\: \$\$f(z) = \frac{L^L}{\mu^L \Gamma(L)}
z^{L-1} \exp(-L z / \mu), \quad z \> 0.\$\$ This is the model for
homogeneous areas (fully developed speckle).

## Usage

``` r
dgammasar(x, mu, L, log = FALSE)

pgammasar(q, mu, L, lower.tail = TRUE, log.p = FALSE)

qgammasar(p, mu, L, lower.tail = TRUE, log.p = FALSE)

rgammasar(n, mu, L)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Mean (\\\mu \> 0\\).

- L:

  Number of looks (\\L \> 0\\).

- log, log.p:

  Logical; if `TRUE`, probabilities/densities are given as logarithms.

- lower.tail:

  Logical; if `TRUE` (default) probabilities are \\P(X \le x)\\,
  otherwise \\P(X \> x)\\.

- p:

  Vector of probabilities.

- n:

  Number of observations.

## Value

Density, distribution function, quantile function or random deviates.

## See also

[GammaSAR](https://jodavid.github.io/SITSAR/reference/GammaSAR.md) for
the R6 interface.

## Examples

``` r
dgammasar(c(0.5, 1, 2), mu = 1, L = 4)
#> [1] 0.7217882 0.7814673 0.1145046
rgammasar(3, mu = 1, L = 4)
#> [1] 0.3494577 1.4005644 0.6856437
```

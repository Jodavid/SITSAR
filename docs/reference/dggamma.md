# The generalized Gamma distribution for SAR intensities

Density, distribution function, quantile function and random generation
for the generalized Gamma (Stacy) law with power \\\alpha \> 0\\, shape
\\\phi \> 0\\ and inverse-scale \\\mu \> 0\\: \$\$f(z) = \frac{\alpha
\mu^{\alpha\phi}}{\Gamma(\phi)} z^{\alpha\phi - 1} \exp\\-(\mu
z)^\alpha\\, \quad z \> 0.\$\$ For \\\alpha = 1\\ it reduces to the
Gamma law with shape \\\phi\\ and rate \\\mu\\; in SAR applications
\\\phi\\ plays the role of the number of looks.

## Usage

``` r
dggamma(x, alpha, phi, mu, log = FALSE)

pggamma(q, alpha, phi, mu, lower.tail = TRUE, log.p = FALSE)

qggamma(p, alpha, phi, mu, lower.tail = TRUE, log.p = FALSE)

rggamma(n, alpha, phi, mu)
```

## Arguments

- x, q:

  Vector of quantiles.

- alpha:

  Power parameter (\\\alpha \> 0\\).

- phi:

  Shape parameter (\\\phi \> 0\\).

- mu:

  Inverse-scale parameter (\\\mu \> 0\\).

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

## References

Li, H.-C., Hong, W., Wu, Y.-R. and Fan, P.-Z. (2011). On the empirical-
statistical modeling of SAR images with generalized gamma distribution.
*IEEE Journal of Selected Topics in Signal Processing*, 5(3), 386–397.

## See also

[GGammaSAR](https://jodavid.github.io/SITSAR/reference/GGammaSAR.md) for
the R6 interface.

## Examples

``` r
dggamma(c(0.5, 1, 2), alpha = 1.5, phi = 2, mu = 1)
#> [1] 0.2633207 0.5518192 0.3546345
rggamma(3, alpha = 1.5, phi = 2, mu = 1)
#> [1] 2.148747 1.881254 1.194354
```

# The G0I distribution for SAR intensities

Density, distribution function, quantile function and random generation
for the \\\mathcal{G}\_I^0(\alpha, \gamma, L)\\ distribution, the
multiplicative model for multilook SAR intensities with reciprocal gamma
backscatter and gamma speckle (Frery et al., 1997).

## Usage

``` r
dg0i(x, alpha, gamma, L, log = FALSE)

pg0i(q, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE)

qg0i(p, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE)

rg0i(n, alpha, gamma, L)
```

## Arguments

- x, q:

  Vector of quantiles.

- alpha:

  Roughness parameter (\\\alpha \< 0\\).

- gamma:

  Scale parameter (\\\gamma \> 0\\).

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

`dg0i()` gives the density, `pg0i()` the distribution function, `qg0i()`
the quantile function and `rg0i()` generates random deviates.

## Details

The density is \$\$f(z) = \frac{L^L \Gamma(L - \alpha)}{\gamma^{\alpha}
\Gamma(-\alpha) \Gamma(L)} \frac{z^{L - 1}}{(\gamma + L z)^{L -
\alpha}}, \quad z \> 0,\$\$ with roughness \\\alpha \< 0\\, scale
\\\gamma \> 0\\ and number of looks \\L \> 0\\. If \\Z \sim
\mathcal{G}\_I^0(\alpha, \gamma, L)\\ then \\-\alpha Z / \gamma\\
follows a Snedecor F distribution with \\2L\\ and \\-2\alpha\\ degrees
of freedom, which is used to compute the distribution and quantile
functions exactly.

## References

Frery, A. C., Muller, H.-J., Yanasse, C. C. F. and Sant'Anna, S. J. S.
(1997). A model for extremely heterogeneous clutter. *IEEE Transactions
on Geoscience and Remote Sensing*, 35(3), 648–659.

Nascimento, A. D. C., Cintra, R. J. and Frery, A. C. (2010). Hypothesis
testing in speckled data with stochastic distances. *IEEE Transactions
on Geoscience and Remote Sensing*, 48(1), 373–385.

## See also

[G0I](https://jodavid.github.io/SITSAR/reference/G0I.md) for the R6
interface (estimation and entropy).

## Examples

``` r
set.seed(1)
z <- rg0i(5, alpha = -1.5, gamma = 2, L = 8)
dg0i(z, alpha = -1.5, gamma = 2, L = 8)
#> [1] 0.41574253 0.19938366 0.44241797 0.28971117 0.04380056
pg0i(c(0.5, 1, 2), alpha = -1.5, gamma = 2, L = 8)
#> [1] 0.08289216 0.29852294 0.58466711
qg0i(c(0.1, 0.5, 0.9), alpha = -1.5, gamma = 2, L = 8)
#> [1] 0.5416068 1.6194935 6.9285395
```

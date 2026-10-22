# Compound complex Wishart distributions

Density and random generation for the compound truncated Poisson complex
Wishart (CTPCW) and compound geometric complex Wishart (CGCW) laws
(Ferreira, Nascimento and Frery, 2022). They describe the PolSAR return
as a random sum \\S = \sum\_{i=1}^N Z_i\\ of independent (unscaled)
complex Wishart matrices \\Z_i \sim \mathcal{W}^C_m(\Sigma, L)\\ with
\\E(Z_i) = L\Sigma\\, where the number of scatterers \\N\\ follows a
zero-truncated Poisson law with mean parameter \\\lambda\\ (CTPCW) or a
geometric law on \\\\1, 2, \ldots\\\\ with success probability \\p\\
(CGCW). Their densities are \$\$f(S) = \sum\_{k=1}^{\infty} P(N = k)
\frac{\|S\|^{kL - m}}{\|\Sigma\|^{kL}\Gamma_m(kL)}
\exp\\-\mathrm{tr}(\Sigma^{-1}S)\\,\$\$ which are able to describe
multimodal intensities. The series is evaluated in log-scale and
truncated when the terms become negligible.

## Usage

``` r
dctpcw(x, lambda, Sigma, L, log = FALSE, kmax = 500L)

rctpcw(n, lambda, Sigma, L, output = c("array", "vec"))

dcgcw(x, p, Sigma, L, log = FALSE, kmax = 500L)

rcgcw(n, p, Sigma, L, output = c("array", "vec"))
```

## Arguments

- x:

  PolSAR observations: a complex array `m x m x n`, a single Hermitian
  matrix or an `n x m(m+1)/2` matrix (see
  [`as_polsar_vec()`](https://jodavid.github.io/SITSAR/reference/polsar-format.md)).

- lambda:

  Mean parameter of the zero-truncated Poisson law (\\\lambda \> 0\\).

- Sigma:

  Hermitian positive definite `m x m` matrix (the mean).

- L:

  Number of looks of each Wishart component (\\L \> m - 1\\).

- log:

  Logical; if `TRUE` the log-density is returned.

- kmax:

  Maximum number of terms of the series.

- n:

  Number of observations.

- output:

  Format of the generated sample: `"array"` (`m x m x n`) or `"vec"`
  (`n x m(m+1)/2`).

- p:

  Success probability of the geometric law (\\0 \< p \< 1\\).

## Value

`dctpcw()`/`dcgcw()` return (log-)densities; `rctpcw()`/`rcgcw()` return
the simulated sample, with the simulated numbers of components in the
attribute `"N"`.

## References

Ferreira, J. A., Nascimento, A. D. C. and Frery, A. C. (2022). PolSAR
models with multimodal intensities. *Remote Sensing*, 14(20), 5083.
[doi:10.3390/rs14205083](https://doi.org/10.3390/rs14205083)

## See also

[CompoundWishart](https://jodavid.github.io/SITSAR/reference/CompoundWishart.md)
for the R6 interface and EM estimation.

## Examples

``` r
Sigma <- polsar_example_sigma("forest")
set.seed(1)
S <- rctpcw(3, lambda = 2, Sigma = Sigma, L = 4)
attr(S, "N")
#> [1] 1 2 2
dctpcw(S, lambda = 2, Sigma = Sigma, L = 4, log = TRUE)
#> [1] -4.110388 -2.837832 -3.709631
S2 <- rcgcw(3, p = 0.6, Sigma = Sigma, L = 4)
dcgcw(S2, p = 0.6, Sigma = Sigma, L = 4, log = TRUE)
#> [1] -7.2704475 -0.7213841 -6.2591386
```

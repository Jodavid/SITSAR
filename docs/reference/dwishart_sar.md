# The scaled complex Wishart distribution

Density and random generation for the scaled complex Wishart law
\\\mathcal{W}\_m(\Sigma, L)\\, the reference model for multilook PolSAR
covariance (coherency) matrices in homogeneous areas: \$\$f(Z) =
\frac{L^{mL} \|Z\|^{L-m}}{\|\Sigma\|^L \Gamma_m(L)} \exp\\-L\\
\mathrm{tr}(\Sigma^{-1} Z)\\,\$\$ where \\\Gamma_m(L) =
\pi^{m(m-1)/2}\prod\_{i=0}^{m-1}\Gamma(L-i)\\ and \\E(Z) = \Sigma\\.

## Usage

``` r
dwishart_sar(x, Sigma, L, log = FALSE)

rwishart_sar(n, Sigma, L, output = c("array", "vec"))
```

## Arguments

- x:

  PolSAR observations: a complex array `m x m x n`, a single Hermitian
  matrix or an `n x m(m+1)/2` matrix (see
  [`as_polsar_vec()`](https://jodavid.github.io/SITSAR/reference/polsar-format.md)).

- Sigma:

  Hermitian positive definite `m x m` matrix (the mean).

- L:

  Number of looks (\\L \> m - 1\\).

- log:

  Logical; if `TRUE` the log-density is returned.

- n:

  Number of observations.

- output:

  Format of the generated sample: `"array"` (`m x m x n`) or `"vec"`
  (`n x m(m+1)/2`).

## Value

`dwishart_sar()` returns a vector of (log-)densities; `rwishart_sar()`
returns the simulated sample.

## Details

Random matrices are generated with the complex Bartlett decomposition,
so non-integer numbers of looks \\L \> m - 1\\ are allowed.

## References

Goodman, N. R. (1963). Statistical analysis based on a certain
multivariate complex Gaussian distribution. *The Annals of Mathematical
Statistics*, 34(1), 152–177.

Lee, J.-S., Grunes, M. R. and Kwok, R. (1994). Classification of
multi-look polarimetric SAR imagery based on complex Wishart
distribution. *International Journal of Remote Sensing*, 15(11),
2299–2311.

## See also

[WishartSAR](https://jodavid.github.io/SITSAR/reference/WishartSAR.md)
for the R6 interface.

## Examples

``` r
Sigma <- polsar_example_sigma("forest")
set.seed(1)
Z <- rwishart_sar(3, Sigma, L = 4)
dim(Z)
#> [1] 3 3 3
dwishart_sar(Z, Sigma, L = 4, log = TRUE)
#> [1] 11.42539 14.85491 13.54976
```

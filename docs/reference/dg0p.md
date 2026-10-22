# The G0 polarimetric distribution

Density and random generation for the polarimetric
\\\mathcal{G}^0_m(\alpha, \Sigma, L)\\ law (also known as
\\\mathcal{G}^0_P\\), obtained from the multiplicative model \\Z = X Y\\
with \\Y \sim \mathcal{W}\_m(\Sigma, L)\\ and reciprocal gamma texture
\\X \sim \Gamma^{-1}(-\alpha, -\alpha - 1)\\ with unit mean (Freitas et
al., 2005): \$\$f(Z) = \frac{L^{mL}\Gamma(mL - \alpha)\|Z\|^{L-m}}
{\|\Sigma\|^L \Gamma_m(L) \Gamma(-\alpha)(-\alpha-1)^{\alpha}}
\left\[L\\\mathrm{tr}(\Sigma^{-1}Z) - \alpha - 1\right\]^{\alpha -
mL},\$\$ with \\\alpha \< -1\\. When \\\alpha \to -\infty\\ it converges
to the scaled complex Wishart law.

## Usage

``` r
dg0p(x, alpha, Sigma, L, log = FALSE)

rg0p(n, alpha, Sigma, L, output = c("array", "vec"))
```

## Arguments

- x:

  PolSAR observations: a complex array `m x m x n`, a single Hermitian
  matrix or an `n x m(m+1)/2` matrix (see
  [`as_polsar_vec()`](https://jodavid.github.io/SITSAR/reference/polsar-format.md)).

- alpha:

  Roughness parameter (\\\alpha \< -1\\).

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

`dg0p()` returns (log-)densities and `rg0p()` the simulated sample.

## References

Freitas, C. C., Frery, A. C. and Correia, A. H. (2005). The polarimetric
G distribution for SAR data analysis. *Environmetrics*, 16(1), 13–31.

Nascimento, A. D. C., Ferreira, J. A. and Frery, A. C. (2023).
Unsupervised segmentation of PolSAR data with complex Wishart and G0m
distributions and Shannon entropy. *Statistics and Computing*, 33, 129.
[doi:10.1007/s11222-023-10298-0](https://doi.org/10.1007/s11222-023-10298-0)

## See also

[G0Polarimetric](https://jodavid.github.io/SITSAR/reference/G0Polarimetric.md)
for the R6 interface.

## Examples

``` r
Sigma <- polsar_example_sigma("urban")
set.seed(1)
Z <- rg0p(3, alpha = -3, Sigma = Sigma, L = 4)
dg0p(Z, alpha = -3, Sigma = Sigma, L = 4, log = TRUE)
#> [1]  7.655394 12.416744  6.848470
```

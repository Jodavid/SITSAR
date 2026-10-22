# Entropy-based test for two SAR samples

Tests whether two samples come from the same distribution by comparing
their estimated Shannon entropies. Under the null hypothesis of equal
entropies, and by the asymptotic normality of the maximum likelihood
estimator of the entropy, \$\$S = \frac{(\hat H_1 - \hat
H_2)^2}{\widehat{\mathrm{Var}}(\hat H_1) + \widehat{\mathrm{Var}}(\hat
H_2)}\$\$ follows asymptotically a \\\chi^2_1\\ distribution. This
statistic is the basis of the entropy-based change detectors of the
thesis of J. A. Ferreira (2021).

## Usage

``` r
entropy_test(
  x,
  y,
  model = c("gamma", "g0i", "ggamma", "wishart", "g0p"),
  L = NULL
)
```

## Arguments

- x, y:

  Samples: numeric vectors of intensities (intensity models) or PolSAR
  observations (matrix models).

- model:

  One of `"gamma"`, `"g0i"`, `"ggamma"`, `"wishart"`, `"g0p"`.

- L:

  Known number of looks; estimated when `NULL`.

## Value

An object of class `htest`.

## References

Ferreira, J. A. (2021). *Entropias e somas compostas para matrizes
estocásticas aplicadas a sistemas PolSAR*. PhD thesis, Universidade
Federal de Pernambuco.

## Examples

``` r
set.seed(1)
x <- rg0i(200, alpha = -3, gamma = 2, L = 4)
y <- rg0i(200, alpha = -3, gamma = 2, L = 4)
z <- rg0i(200, alpha = -1.5, gamma = 1, L = 4)
entropy_test(x, y, model = "g0i", L = 4)
#> 
#>  Shannon entropy test (G0I model)
#> 
#> data:  x and y
#> S = 2.989, df = 1, p-value = 0.08383
#> alternative hypothesis: true difference in entropies is not equal to 0
#> sample estimates:
#> entropy of x entropy of y 
#>    0.7837866    0.9918088 
#> 
entropy_test(x, z, model = "g0i", L = 4)$p.value
#> [1] 0.001877851
```

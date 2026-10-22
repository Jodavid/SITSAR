# Scaled complex Wishart model for PolSAR data (R6)

R6 class for the scaled complex Wishart law \\\mathcal{W}\_m(\Sigma,
L)\\ (see
[`dwishart_sar()`](https://jodavid.github.io/SITSAR/reference/dwishart_sar.md)).
The maximum likelihood estimator of \\\Sigma\\ is the sample mean and
\\L\\ solves \\m\log L - \psi^{(0)}\_m(L) = \log\|\bar Z\| -
\overline{\log\|Z\|}\\ (Anfinsen et al., 2009). The entropy variance is
available in closed form (Nascimento, Ferreira and Frery, 2023).

## References

Anfinsen, S. N., Doulgeris, A. P. and Eltoft, T. (2009). Estimation of
the equivalent number of looks in polarimetric synthetic aperture radar
imagery. *IEEE Transactions on Geoscience and Remote Sensing*, 47(11),
3795–3809.

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`WishartSAR`

## Methods

### Public methods

- [`WishartSAR$new()`](#method-WishartSAR-initialize)

- [`WishartSAR$density()`](#method-WishartSAR-density)

- [`WishartSAR$random()`](#method-WishartSAR-random)

- [`WishartSAR$fit()`](#method-WishartSAR-fit)

- [`WishartSAR$entropy()`](#method-WishartSAR-entropy)

- [`WishartSAR$entropy_variance()`](#method-WishartSAR-entropy_variance)

- [`WishartSAR$clone()`](#method-WishartSAR-clone)

Inherited methods

- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `WishartSAR$new()`

Create a Wishart model.

#### Usage

    WishartSAR$new(Sigma = diag(3) + (0+0i), L = 4)

#### Arguments

- `Sigma`:

  Hermitian positive definite matrix (default identity 3x3).

- `L`:

  Number of looks.

------------------------------------------------------------------------

### `WishartSAR$density()`

(Log-)density.

#### Usage

    WishartSAR$density(x, log = FALSE)

#### Arguments

- `x`:

  PolSAR observations (array `m x m x n` or vec matrix).

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `WishartSAR$random()`

Simulate PolSAR observations.

#### Usage

    WishartSAR$random(n, output = "array")

#### Arguments

- `n`:

  Number of observations.

- `output`:

  `"array"` or `"vec"`.

------------------------------------------------------------------------

### `WishartSAR$fit()`

Fit by maximum likelihood.

#### Usage

    WishartSAR$fit(x, L = NULL)

#### Arguments

- `x`:

  PolSAR observations (array `m x m x n` or vec matrix).

- `L`:

  Known number of looks. If `NULL` (default) it is estimated.

------------------------------------------------------------------------

### `WishartSAR$entropy()`

Shannon entropy.

#### Usage

    WishartSAR$entropy()

------------------------------------------------------------------------

### `WishartSAR$entropy_variance()`

Asymptotic variance of the entropy estimator (closed form).

#### Usage

    WishartSAR$entropy_variance(n = self$n)

#### Arguments

- `n`:

  Sample size.

------------------------------------------------------------------------

### `WishartSAR$clone()`

The objects of this class are cloneable with this method.

#### Usage

    WishartSAR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
Z <- rwishart_sar(300, polsar_example_sigma("forest"), L = 4)
m <- WishartSAR$new()$fit(Z)
m
#> <Scaled complex Wishart> fitted to n = 300 observations
#>   Sigma:
#>                  [,1]              [,2]              [,3]
#> [1,] 0.34230+0.00000i 0.013160+0.00643i 0.059910-0.00602i
#> [2,] 0.01316-0.00643i 0.101100+0.00000i 0.001538+0.00261i
#> [3,] 0.05991+0.00602i 0.001538-0.00261i 0.212300+0.00000i
#>   L = 3.97893
#>   logLik = 3618.44
#>   Shannon entropy = -12.0615
m$entropy_variance()
#> [1] 0.02448245
```

# G0 polarimetric model for PolSAR data (R6)

R6 class for the \\\mathcal{G}^0_m(\alpha, \Sigma, L)\\ law (see
[`dg0p()`](https://jodavid.github.io/SITSAR/reference/dg0p.md)).
Following Nascimento, Ferreira and Frery (2023), \\\Sigma\\ is estimated
by the sample mean (since \\E(Z) = \Sigma\\) and \\(\alpha, L)\\ by
maximising the likelihood (analytic score). The variance of the entropy
estimator combines the inverse observed information of \\(\alpha, L)\\
with a moment estimate of the variance of \\m\log\|\hat\Sigma\|\\.

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`G0Polarimetric`

## Methods

### Public methods

- [`G0Polarimetric$new()`](#method-G0Polarimetric-initialize)

- [`G0Polarimetric$density()`](#method-G0Polarimetric-density)

- [`G0Polarimetric$random()`](#method-G0Polarimetric-random)

- [`G0Polarimetric$fit()`](#method-G0Polarimetric-fit)

- [`G0Polarimetric$entropy()`](#method-G0Polarimetric-entropy)

- [`G0Polarimetric$clone()`](#method-G0Polarimetric-clone)

Inherited methods

- [`SARModel$entropy_variance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-entropy_variance)
- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `G0Polarimetric$new()`

Create a G0 polarimetric model.

#### Usage

    G0Polarimetric$new(alpha = -3, Sigma = diag(3) + (0+0i), L = 4)

#### Arguments

- `alpha`:

  Roughness (\< -1).

- `Sigma`:

  Hermitian positive definite matrix.

- `L`:

  Number of looks.

------------------------------------------------------------------------

### `G0Polarimetric$density()`

(Log-)density.

#### Usage

    G0Polarimetric$density(x, log = FALSE)

#### Arguments

- `x`:

  PolSAR observations.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `G0Polarimetric$random()`

Simulate PolSAR observations.

#### Usage

    G0Polarimetric$random(n, output = "array")

#### Arguments

- `n`:

  Number of observations.

- `output`:

  `"array"` or `"vec"`.

------------------------------------------------------------------------

### `G0Polarimetric$fit()`

Fit by maximum likelihood.

#### Usage

    G0Polarimetric$fit(x, L = NULL)

#### Arguments

- `x`:

  PolSAR observations.

- `L`:

  Known number of looks. If `NULL` (default) it is estimated.

------------------------------------------------------------------------

### `G0Polarimetric$entropy()`

Shannon entropy.

#### Usage

    G0Polarimetric$entropy()

------------------------------------------------------------------------

### `G0Polarimetric$clone()`

The objects of this class are cloneable with this method.

#### Usage

    G0Polarimetric$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
Z <- rg0p(400, alpha = -4, Sigma = polsar_example_sigma("urban"), L = 4)
m <- G0Polarimetric$new()$fit(Z, L = 4)
m
#> <G0 polarimetric> fitted to n = 400 observations
#>   alpha = -4.56299
#>   Sigma:
#>                   [,1]              [,2]             [,3]
#> [1,]  0.88150+0.00000i  0.02508+0.00207i -0.15380+0.1836i
#> [2,]  0.02508-0.00207i  0.05535+0.00000i -0.00867+0.0141i
#> [3,] -0.15380-0.18360i -0.00867-0.01410i  0.47100+0.0000i
#>   L = 4 (fixed)
#>   logLik = 3778.55
#>   Shannon entropy = -9.51015
m$entropy_variance()
#> [1] 0.0976793
```

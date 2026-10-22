# Gamma model for SAR intensities (R6)

R6 class for the multilook Gamma law \\\Gamma(\mu, L)\\; see
[`dgammasar()`](https://jodavid.github.io/SITSAR/reference/dgammasar.md).
The maximum likelihood estimator of \\\mu\\ is the sample mean and \\L\\
solves \\\log L - \psi^{(0)}(L) = \log\bar z - \overline{\log z}\\.

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`GammaSAR`

## Methods

### Public methods

- [`GammaSAR$new()`](#method-GammaSAR-initialize)

- [`GammaSAR$density()`](#method-GammaSAR-density)

- [`GammaSAR$random()`](#method-GammaSAR-random)

- [`GammaSAR$fit()`](#method-GammaSAR-fit)

- [`GammaSAR$entropy()`](#method-GammaSAR-entropy)

- [`GammaSAR$entropy_variance()`](#method-GammaSAR-entropy_variance)

- [`GammaSAR$clone()`](#method-GammaSAR-clone)

Inherited methods

- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `GammaSAR$new()`

Create a Gamma model.

#### Usage

    GammaSAR$new(mu = 1, L = 1)

#### Arguments

- `mu`:

  Mean.

- `L`:

  Number of looks.

------------------------------------------------------------------------

### `GammaSAR$density()`

(Log-)density.

#### Usage

    GammaSAR$density(x, log = FALSE)

#### Arguments

- `x`:

  Numeric vector of intensities.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `GammaSAR$random()`

Simulate intensities.

#### Usage

    GammaSAR$random(n)

#### Arguments

- `n`:

  Number of observations.

------------------------------------------------------------------------

### `GammaSAR$fit()`

Fit by maximum likelihood.

#### Usage

    GammaSAR$fit(x, L = NULL)

#### Arguments

- `x`:

  Numeric vector of positive intensities.

- `L`:

  Known number of looks. If `NULL` (default) it is estimated.

------------------------------------------------------------------------

### `GammaSAR$entropy()`

Shannon entropy.

#### Usage

    GammaSAR$entropy()

------------------------------------------------------------------------

### `GammaSAR$entropy_variance()`

Asymptotic variance of the entropy estimator (closed form).

#### Usage

    GammaSAR$entropy_variance(n = self$n)

#### Arguments

- `n`:

  Sample size.

------------------------------------------------------------------------

### `GammaSAR$clone()`

The objects of this class are cloneable with this method.

#### Usage

    GammaSAR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
z <- rgammasar(500, mu = 2, L = 4)
m <- GammaSAR$new()$fit(z)
m
#> <Gamma> fitted to n = 500 observations
#>   mu = 1.97692
#>   L = 3.86696
#>   logLik = -666.161
#>   Shannon entropy = 1.33232
m$entropy()
#> [1] 1.332323
m$entropy_variance()
#> [1] 0.001113651
generics::tidy(m)
#> # A tibble: 2 × 4
#>   term  estimate std.error fixed
#>   <chr>    <dbl>     <dbl> <lgl>
#> 1 mu        1.98    0.0450 FALSE
#> 2 L         3.87    0.235  FALSE
```

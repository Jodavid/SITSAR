# Generalized Gamma model for SAR intensities (R6)

R6 class for the generalized Gamma law (see
[`dggamma()`](https://jodavid.github.io/SITSAR/reference/dggamma.md)).
Estimation is carried out by profiling: for a given power \\\alpha\\,
the maximum likelihood estimates of \\\phi\\ and \\\mu\\ are available
in closed form (up to a one-dimensional root), so the profile likelihood
is maximised over \\\alpha\\ only.

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`GGammaSAR`

## Methods

### Public methods

- [`GGammaSAR$new()`](#method-GGammaSAR-initialize)

- [`GGammaSAR$density()`](#method-GGammaSAR-density)

- [`GGammaSAR$random()`](#method-GGammaSAR-random)

- [`GGammaSAR$fit()`](#method-GGammaSAR-fit)

- [`GGammaSAR$entropy()`](#method-GGammaSAR-entropy)

- [`GGammaSAR$clone()`](#method-GGammaSAR-clone)

Inherited methods

- [`SARModel$entropy_variance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-entropy_variance)
- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `GGammaSAR$new()`

Create a generalized Gamma model.

#### Usage

    GGammaSAR$new(alpha = 1, phi = 1, mu = 1)

#### Arguments

- `alpha`:

  Power (\> 0).

- `phi`:

  Shape (\> 0).

- `mu`:

  Inverse scale (\> 0).

------------------------------------------------------------------------

### `GGammaSAR$density()`

(Log-)density.

#### Usage

    GGammaSAR$density(x, log = FALSE)

#### Arguments

- `x`:

  Numeric vector of intensities.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `GGammaSAR$random()`

Simulate intensities.

#### Usage

    GGammaSAR$random(n)

#### Arguments

- `n`:

  Number of observations.

------------------------------------------------------------------------

### `GGammaSAR$fit()`

Fit by (profile) maximum likelihood.

#### Usage

    GGammaSAR$fit(x, phi = NULL)

#### Arguments

- `x`:

  Numeric vector of positive intensities.

- `phi`:

  Known shape (number of looks). If `NULL` it is estimated.

------------------------------------------------------------------------

### `GGammaSAR$entropy()`

Shannon entropy.

#### Usage

    GGammaSAR$entropy()

------------------------------------------------------------------------

### `GGammaSAR$clone()`

The objects of this class are cloneable with this method.

#### Usage

    GGammaSAR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
z <- rggamma(1000, alpha = 1.5, phi = 3, mu = 2)
m <- GGammaSAR$new()$fit(z)
m
#> <Generalized Gamma> fitted to n = 1000 observations
#>   alpha = 1.34428
#>   phi = 3.4875
#>   mu = 2.47093
#>   logLik = -458.979
#>   Shannon entropy = 0.458979
```

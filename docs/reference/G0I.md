# G0I model for SAR intensities (R6)

R6 class for the \\\mathcal{G}\_I^0(\alpha, \gamma, L)\\ law (see
[`dg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md)) with
maximum likelihood estimation (analytic score, quasi-Newton optimisation
on a transformed scale with several starting points), closed-form
Shannon entropy (Ferreira and Nascimento, 2020) and the asymptotic
variance of its estimator.

## References

Ferreira, J. A. and Nascimento, A. D. C. (2020). Shannon entropy for the
G0I model: a new segmentation approach. *IEEE Journal of Selected Topics
in Applied Earth Observations and Remote Sensing*, 13, 2547–2553.
[doi:10.1109/JSTARS.2020.2997666](https://doi.org/10.1109/JSTARS.2020.2997666)

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`G0I`

## Methods

### Public methods

- [`G0I$new()`](#method-G0I-initialize)

- [`G0I$density()`](#method-G0I-density)

- [`G0I$random()`](#method-G0I-random)

- [`G0I$fit()`](#method-G0I-fit)

- [`G0I$entropy()`](#method-G0I-entropy)

- [`G0I$clone()`](#method-G0I-clone)

Inherited methods

- [`SARModel$entropy_variance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-entropy_variance)
- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `G0I$new()`

Create a G0I model.

#### Usage

    G0I$new(alpha = -2, gamma = 1, L = 1)

#### Arguments

- `alpha`:

  Roughness (\< 0).

- `gamma`:

  Scale (\> 0).

- `L`:

  Number of looks (\> 0).

------------------------------------------------------------------------

### `G0I$density()`

(Log-)density.

#### Usage

    G0I$density(x, log = FALSE)

#### Arguments

- `x`:

  Numeric vector of intensities.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `G0I$random()`

Simulate intensities.

#### Usage

    G0I$random(n)

#### Arguments

- `n`:

  Number of observations.

------------------------------------------------------------------------

### `G0I$fit()`

Fit by maximum likelihood.

#### Usage

    G0I$fit(x, L = NULL, start = NULL)

#### Arguments

- `x`:

  Numeric vector of positive intensities.

- `L`:

  Known number of looks. If `NULL` (default) it is estimated.

- `start`:

  Optional named vector of starting values (`alpha`, `gamma` and, if
  estimated, `L`).

------------------------------------------------------------------------

### `G0I$entropy()`

Shannon entropy.

#### Usage

    G0I$entropy()

------------------------------------------------------------------------

### `G0I$clone()`

The objects of this class are cloneable with this method.

#### Usage

    G0I$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
z <- rg0i(1000, alpha = -3, gamma = 2, L = 4)
m <- G0I$new()$fit(z, L = 4)
m
#> <G0I> fitted to n = 1000 observations
#>   alpha = -3.13702
#>   gamma = 2.0585
#>   L = 4 (fixed)
#>   logLik = -838.747
#>   Shannon entropy = 0.824063
m$entropy_variance()
#> [1] 0.001294807
G0I$new()$fit(z)$params # L estimated
#> $alpha
#> [1] -3.957879
#> 
#> $gamma
#> [1] 2.803738
#> 
#> $L
#> [1] 3.006068
#> 
```

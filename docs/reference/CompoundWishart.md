# Compound complex Wishart models (R6)

R6 class for the compound truncated Poisson (CTPCW) and compound
geometric (CGCW) complex Wishart laws of Ferreira, Nascimento and Frery
(2022); see
[`dctpcw()`](https://jodavid.github.io/SITSAR/reference/compound_wishart.md).
Parameters are estimated by the expectation-maximisation (EM) algorithm,
treating the number of scatterers \\N\\ as missing data:

- E-step: \\\hat N_i = E(N \mid S_i)\\;

- M-step: \\\hat\Sigma = \bar S / (L \bar N)\\ and \\\hat\lambda/(1 -
  e^{-\hat\lambda}) = \bar N\\ (CTPCW) or \\\hat p = 1/\bar N\\ (CGCW).

The Shannon entropy has no closed form and is approximated by Monte
Carlo.

## References

Ferreira, J. A., Nascimento, A. D. C. and Frery, A. C. (2022). PolSAR
models with multimodal intensities. *Remote Sensing*, 14(20), 5083.
[doi:10.3390/rs14205083](https://doi.org/10.3390/rs14205083)

## Super class

[`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) -\>
`CompoundWishart`

## Public fields

- `counting`:

  Counting law: `"poisson"` (zero-truncated) or `"geometric"`.

- `iterations`:

  Number of EM iterations of the last fit.

## Methods

### Public methods

- [`CompoundWishart$new()`](#method-CompoundWishart-initialize)

- [`CompoundWishart$density()`](#method-CompoundWishart-density)

- [`CompoundWishart$random()`](#method-CompoundWishart-random)

- [`CompoundWishart$fit()`](#method-CompoundWishart-fit)

- [`CompoundWishart$entropy()`](#method-CompoundWishart-entropy)

- [`CompoundWishart$entropy_variance()`](#method-CompoundWishart-entropy_variance)

- [`CompoundWishart$clone()`](#method-CompoundWishart-clone)

Inherited methods

- [`SARModel$glance()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-glance)
- [`SARModel$log_likelihood()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-log_likelihood)
- [`SARModel$n_parameters()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-n_parameters)
- [`SARModel$print()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-print)
- [`SARModel$tidy()`](https://jodavid.github.io/SITSAR/reference/SARModel.html#method-tidy)

------------------------------------------------------------------------

### `CompoundWishart$new()`

Create a compound Wishart model.

#### Usage

    CompoundWishart$new(
      counting = c("poisson", "geometric"),
      Sigma = diag(3) + (0+0i),
      L = 4,
      lambda = 1,
      p = 0.5
    )

#### Arguments

- `counting`:

  `"poisson"` (CTPCW) or `"geometric"` (CGCW).

- `Sigma`:

  Hermitian positive definite matrix.

- `L`:

  Number of looks of each component.

- `lambda, p`:

  Parameter of the counting law.

------------------------------------------------------------------------

### `CompoundWishart$density()`

(Log-)density.

#### Usage

    CompoundWishart$density(x, log = FALSE)

#### Arguments

- `x`:

  PolSAR observations.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `CompoundWishart$random()`

Simulate PolSAR observations.

#### Usage

    CompoundWishart$random(n, output = "array")

#### Arguments

- `n`:

  Number of observations.

- `output`:

  `"array"` or `"vec"`.

------------------------------------------------------------------------

### `CompoundWishart$fit()`

Fit by the EM algorithm.

#### Usage

    CompoundWishart$fit(x, L, tol = 1e-06, max_iter = 500L)

#### Arguments

- `x`:

  PolSAR observations.

- `L`:

  Known number of looks of each component.

- `tol`:

  Convergence tolerance (relative change of the parameters).

- `max_iter`:

  Maximum number of iterations.

------------------------------------------------------------------------

### `CompoundWishart$entropy()`

Monte Carlo approximation of the Shannon entropy.

#### Usage

    CompoundWishart$entropy(nsim = 2000L)

#### Arguments

- `nsim`:

  Number of simulated observations.

------------------------------------------------------------------------

### `CompoundWishart$entropy_variance()`

Not available for compound models.

#### Usage

    CompoundWishart$entropy_variance(n = self$n)

#### Arguments

- `n`:

  Ignored.

------------------------------------------------------------------------

### `CompoundWishart$clone()`

The objects of this class are cloneable with this method.

#### Usage

    CompoundWishart$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
S <- rctpcw(200, lambda = 1.5, Sigma = polsar_example_sigma("forest"), L = 4)
m <- CompoundWishart$new(counting = "poisson")$fit(S, L = 4)
m$params$lambda
#> [1] 1.462612
m$iterations
#> [1] 33
```

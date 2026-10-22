# Base R6 class for SAR and PolSAR statistical models

`SARModel` is the abstract parent of all model classes of `SITSAR`
([GammaSAR](https://jodavid.github.io/SITSAR/reference/GammaSAR.md),
[G0I](https://jodavid.github.io/SITSAR/reference/G0I.md),
[GGammaSAR](https://jodavid.github.io/SITSAR/reference/GGammaSAR.md),
[WishartSAR](https://jodavid.github.io/SITSAR/reference/WishartSAR.md),
[G0Polarimetric](https://jodavid.github.io/SITSAR/reference/G0Polarimetric.md)
and
[CompoundWishart](https://jodavid.github.io/SITSAR/reference/CompoundWishart.md)).
Every model can evaluate densities, simulate data, be fitted by maximum
likelihood (`$fit()`) and report its Shannon entropy and the asymptotic
variance of the entropy estimator.

Fitted models work with the `generics` verbs
[tidy()](https://generics.r-lib.org/reference/tidy.html),
[glance()](https://generics.r-lib.org/reference/glance.html), as well as
with [`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html).

## Public fields

- `name`:

  Model name.

- `params`:

  Named list of parameters.

- `fixed`:

  Names of parameters held fixed during estimation.

- `n`:

  Sample size used in the last fit (`NULL` if not fitted).

- `loglik`:

  Maximised log-likelihood (`NULL` if not fitted).

- `converged`:

  Logical convergence flag of the last fit.

## Methods

### Public methods

- [`SARModel$new()`](#method-SARModel-initialize)

- [`SARModel$density()`](#method-SARModel-density)

- [`SARModel$random()`](#method-SARModel-random)

- [`SARModel$fit()`](#method-SARModel-fit)

- [`SARModel$log_likelihood()`](#method-SARModel-log_likelihood)

- [`SARModel$entropy()`](#method-SARModel-entropy)

- [`SARModel$entropy_variance()`](#method-SARModel-entropy_variance)

- [`SARModel$n_parameters()`](#method-SARModel-n_parameters)

- [`SARModel$tidy()`](#method-SARModel-tidy)

- [`SARModel$glance()`](#method-SARModel-glance)

- [`SARModel$print()`](#method-SARModel-print)

- [`SARModel$clone()`](#method-SARModel-clone)

------------------------------------------------------------------------

### `SARModel$new()`

Abstract constructor.

#### Usage

    SARModel$new(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `SARModel$density()`

(Log-)density of observations.

#### Usage

    SARModel$density(x, log = FALSE)

#### Arguments

- `x`:

  Observations.

- `log`:

  Return log-densities?

------------------------------------------------------------------------

### `SARModel$random()`

Generate random observations.

#### Usage

    SARModel$random(n)

#### Arguments

- `n`:

  Number of observations.

------------------------------------------------------------------------

### `SARModel$fit()`

Fit the model by maximum likelihood.

#### Usage

    SARModel$fit(x, ...)

#### Arguments

- `x`:

  Observations.

- `...`:

  Model specific arguments.

------------------------------------------------------------------------

### `SARModel$log_likelihood()`

Log-likelihood of a sample.

#### Usage

    SARModel$log_likelihood(x)

#### Arguments

- `x`:

  Observations.

------------------------------------------------------------------------

### `SARModel$entropy()`

Shannon entropy at the current parameters.

#### Usage

    SARModel$entropy()

------------------------------------------------------------------------

### `SARModel$entropy_variance()`

Asymptotic variance of the entropy estimator (delta method with the
inverse observed Fisher information).

#### Usage

    SARModel$entropy_variance(n = self$n)

#### Arguments

- `n`:

  Sample size. Defaults to the size of the fitted sample.

------------------------------------------------------------------------

### `SARModel$n_parameters()`

Number of free parameters.

#### Usage

    SARModel$n_parameters()

------------------------------------------------------------------------

### `SARModel$tidy()`

Parameters as a tibble.

#### Usage

    SARModel$tidy()

------------------------------------------------------------------------

### `SARModel$glance()`

One-row summary of the fit.

#### Usage

    SARModel$glance()

------------------------------------------------------------------------

### `SARModel$print()`

Print the model.

#### Usage

    SARModel$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `SARModel$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SARModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

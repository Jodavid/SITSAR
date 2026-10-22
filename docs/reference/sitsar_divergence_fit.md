# Engine fit function for `divergence_classifier()`

Internal fit function used by the 'parsnip' engines; it is exported only
so that 'parsnip' can call it.

## Usage

``` r
sitsar_divergence_fit(
  x,
  y,
  family = "gaussian",
  divergence = "renyi",
  order = 0.5,
  kind = 1
)
```

## Arguments

- x:

  Data frame of predictors.

- y:

  Factor of classes.

- family:

  `"gaussian"` or `"wishart"`.

- divergence, order, kind:

  See
  [`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md).

## Value

A fitted
[DivergenceClassifier](https://jodavid.github.io/SITSAR/reference/DivergenceClassifier.md).

# Divergence-based classifier as a 'parsnip' model

`divergence_classifier()` defines the divergence-based classifier of
Ferreira, Coelho and Nascimento (2021) as a 'parsnip' model
specification, so it can be used with 'workflows', 'recipes', 'rsample',
'tune' and 'yardstick'. The engine selects the family of features:
`"gaussian"` (multivariate normal intensity features, the default) or
`"wishart"` (complex Wishart PolSAR matrix features). The predictors
must be in the format of
[`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md);
extra columns such as `row`/`col` are ignored by the engine (it is good
practice to give them an `"id"` role in a recipe).

## Usage

``` r
divergence_classifier(
  mode = "classification",
  divergence = NULL,
  order = NULL,
  kind = NULL,
  engine = "gaussian"
)

# S3 method for class 'divergence_classifier'
update(
  object,
  parameters = NULL,
  divergence = NULL,
  order = NULL,
  kind = NULL,
  fresh = FALSE,
  ...
)
```

## Arguments

- mode:

  A single character string for the type of model. The only possible
  value is `"classification"`.

- divergence:

  Divergence: `"renyi"` (default), `"kl"`, `"jeffreys"` or
  `"bhattacharyya"`.

- order:

  Order of the Rényi divergence, in (0, 1). Tunable with
  [`renyi_order()`](https://jodavid.github.io/SITSAR/reference/renyi_order.md).

- kind:

  Form of the Rényi divergence (1 to 4); see
  [divergences](https://jodavid.github.io/SITSAR/reference/divergences.md).

- engine:

  `"gaussian"` or `"wishart"`.

- object:

  A `divergence_classifier` model specification.

- parameters:

  A 1-row tibble or named list with the main parameters.

- fresh:

  Replace the arguments instead of updating them?

- ...:

  Not used.

## Value

A model specification object with classes `divergence_classifier` and
`model_spec`.

## See also

[DivergenceClassifier](https://jodavid.github.io/SITSAR/reference/DivergenceClassifier.md),
[`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md)

## Examples

``` r
spec <- divergence_classifier(divergence = "renyi", order = 0.5) |>
  parsnip::set_engine("wishart")
spec
#> Divergence Classifier Model Specification (classification)
#> 
#> Main Arguments:
#>   divergence = renyi
#>   order = 0.5
#> 
#> Computational engine: wishart 
#> 

feats <- sar_features(SARImage$new(SanFrancisco150)$crop(1:60, 1:60), window = 5)
feats$class <- factor(as.vector(SanFrancisco150_gt[1:60, 1:60]))
fit <- parsnip::fit(spec, class ~ ., data = feats)
head(predict(fit, feats))
#> # A tibble: 6 × 1
#>   .pred_class
#>   <fct>      
#> 1 1          
#> 2 1          
#> 3 1          
#> 4 1          
#> 5 1          
#> 6 1          
head(predict(fit, feats, type = "prob"))
#> # A tibble: 6 × 1
#>   .pred_1
#>     <dbl>
#> 1       1
#> 2       1
#> 3       1
#> 4       1
#> 5       1
#> 6       1
```

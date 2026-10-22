# Divergence-based classifier for SAR/PolSAR features (R6)

Implements the family of divergence-based classifiers of Ferreira,
Coelho and Nascimento (2021). Each pixel is represented by the model
fitted to the observations of a moving window (see
[`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md))
and each class by the model fitted to its training pixels. A pixel is
assigned to the class that minimises the chosen divergence between the
two fitted models: \$\$\hat\imath = \arg\min\_{i \in \\1, \ldots, g\\}
d(\hat\theta_h, \hat\theta_i).\$\$

Two families of features are available:

- `"gaussian"`:

  multivariate normal laws for the intensity vector, using the columns
  `mean_*` and `cov_*` of
  [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md).
  The class covariance combines the average local covariance and the
  covariance of the local means (law of total variance).

- `"wishart"`:

  scaled complex Wishart laws for the full PolSAR covariance matrix,
  using the columns `c11`, `c12_re`, `c12_im`, ... of
  [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md).
  Since all laws share the number of looks, it only rescales the
  divergences and does not affect the classification.

The classifier is also available as the 'parsnip' model
[`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md),
for use in 'tidymodels' workflows.

## References

Ferreira, J. A., Coelho, H. and Nascimento, A. D. C. (2021). A family of
divergence-based classifiers for polarimetric synthetic aperture radar
(PolSAR) imagery vector and matrix features. *International Journal of
Remote Sensing*, 42(4), 1201–1229.
[doi:10.1080/01431161.2020.1826060](https://doi.org/10.1080/01431161.2020.1826060)

## Public fields

- `family`:

  `"gaussian"` or `"wishart"`.

- `divergence`:

  Divergence used (see
  [`divergence_gaussian()`](https://jodavid.github.io/SITSAR/reference/divergences.md)).

- `order`:

  Order of the Rényi divergence.

- `kind`:

  Form (1-4) of the Rényi divergence.

- `classes`:

  Class labels.

- `prototypes`:

  List with the fitted model of each class.

- `channels`:

  Channels (Gaussian) or matrix order (Wishart) used.

## Methods

### Public methods

- [`DivergenceClassifier$new()`](#method-DivergenceClassifier-initialize)

- [`DivergenceClassifier$fit()`](#method-DivergenceClassifier-fit)

- [`DivergenceClassifier$distances()`](#method-DivergenceClassifier-distances)

- [`DivergenceClassifier$predict()`](#method-DivergenceClassifier-predict)

- [`DivergenceClassifier$print()`](#method-DivergenceClassifier-print)

- [`DivergenceClassifier$clone()`](#method-DivergenceClassifier-clone)

------------------------------------------------------------------------

### `DivergenceClassifier$new()`

Create a classifier.

#### Usage

    DivergenceClassifier$new(
      family = c("gaussian", "wishart"),
      divergence = c("renyi", "kl", "jeffreys", "bhattacharyya"),
      order = 0.5,
      kind = 1
    )

#### Arguments

- `family`:

  `"gaussian"` or `"wishart"`.

- `divergence`:

  `"kl"`, `"jeffreys"`, `"bhattacharyya"` or `"renyi"`.

- `order`:

  Order of the Rényi divergence, in (0, 1).

- `kind`:

  Form of the Rényi divergence (1 to 4).

------------------------------------------------------------------------

### `DivergenceClassifier$fit()`

Fit the class prototypes.

#### Usage

    DivergenceClassifier$fit(x, y)

#### Arguments

- `x`:

  Data frame of features (see
  [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md));
  other columns are ignored.

- `y`:

  Factor (or vector) of class labels.

------------------------------------------------------------------------

### `DivergenceClassifier$distances()`

Divergence between every row of `x` and every class.

#### Usage

    DivergenceClassifier$distances(x)

#### Arguments

- `x`:

  Data frame of features.

------------------------------------------------------------------------

### `DivergenceClassifier$predict()`

Predict classes, pseudo-probabilities or distances.

#### Usage

    DivergenceClassifier$predict(x, type = c("class", "prob", "distance"))

#### Arguments

- `x`:

  Data frame of features.

- `type`:

  `"class"`, `"prob"` (softmax of minus the divergences) or
  `"distance"`.

------------------------------------------------------------------------

### `DivergenceClassifier$print()`

Print the classifier.

#### Usage

    DivergenceClassifier$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `DivergenceClassifier$clone()`

The objects of this class are cloneable with this method.

#### Usage

    DivergenceClassifier$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
img <- SARImage$new(SanFrancisco150)$crop(1:80, 1:80)
feats <- sar_features(img, window = 5)
feats$class <- factor(as.vector(SanFrancisco150_gt[1:80, 1:80]),
                      labels = c("ocean", "forest", "urban"))
set.seed(1)
train <- feats[sample(nrow(feats), 300), ]

clf <- DivergenceClassifier$new(family = "wishart", divergence = "renyi",
                                order = 0.5, kind = 3)
clf$fit(train, train$class)
pred <- clf$predict(feats)
mean(pred == feats$class)
#> [1] 0.9920313
head(clf$predict(feats, type = "distance"))
#> # A tibble: 6 × 3
#>   ocean forest urban
#>   <dbl>  <dbl> <dbl>
#> 1  1.72   6.56  9.48
#> 2  1.68   6.55  9.50
#> 3  1.65   6.53  9.50
#> 4  1.69   6.59  9.58
#> 5  1.71   6.66  9.68
#> 6  1.62   6.51  9.51
```

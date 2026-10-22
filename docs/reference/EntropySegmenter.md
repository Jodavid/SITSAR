# Entropy-based segmentation by (stochastic) EM (R6)

Unsupervised segmentation of SAR/PolSAR images from maps of estimated
Shannon entropies (Ferreira and Nascimento, 2020; Nascimento, Ferreira
and Frery, 2023). By the asymptotic normality of the entropy estimators,
the entropy vector of a pixel of the \\g\\-th region is approximately
Gaussian, so the image is modelled as a finite mixture of Gaussian laws
(with diagonal covariance matrices) fitted by the
expectation-maximisation algorithm (`method = "em"`) or by its
stochastic version (`method = "sem"`), in which the E-step is followed
by a random allocation of every pixel to a class drawn from its
posterior probabilities. Pixels are then labelled by the maximum a
posteriori rule.

## References

Ferreira, J. A. and Nascimento, A. D. C. (2020). Shannon entropy for the
G0I model: a new segmentation approach. *IEEE Journal of Selected Topics
in Applied Earth Observations and Remote Sensing*, 13, 2547–2553.
[doi:10.1109/JSTARS.2020.2997666](https://doi.org/10.1109/JSTARS.2020.2997666)

Nascimento, A. D. C., Ferreira, J. A. and Frery, A. C. (2023).
Unsupervised segmentation of PolSAR data with complex Wishart and G0m
distributions and Shannon entropy. *Statistics and Computing*, 33, 129.
[doi:10.1007/s11222-023-10298-0](https://doi.org/10.1007/s11222-023-10298-0)

## Public fields

- `k`:

  Number of classes.

- `method`:

  `"em"` or `"sem"`.

- `weights`:

  Mixture weights.

- `means`:

  Matrix (`k x features`) of class means.

- `variances`:

  Matrix (`k x features`) of class variances.

- `posterior`:

  Matrix of posterior probabilities.

- `labels`:

  Integer vector of labels (classes ordered by the mean of the first
  feature).

- `loglik`:

  Trace of the log-likelihood.

- `iterations`:

  Number of iterations.

- `converged`:

  Logical convergence flag.

- `dim`:

  Spatial dimensions when fitted to a map (`NULL` otherwise).

## Methods

### Public methods

- [`EntropySegmenter$new()`](#method-EntropySegmenter-initialize)

- [`EntropySegmenter$fit()`](#method-EntropySegmenter-fit)

- [`EntropySegmenter$predict()`](#method-EntropySegmenter-predict)

- [`EntropySegmenter$label_map()`](#method-EntropySegmenter-label_map)

- [`EntropySegmenter$print()`](#method-EntropySegmenter-print)

- [`EntropySegmenter$clone()`](#method-EntropySegmenter-clone)

------------------------------------------------------------------------

### `EntropySegmenter$new()`

Create a segmenter.

#### Usage

    EntropySegmenter$new(
      k = 2,
      method = c("em", "sem"),
      max_iter = 200L,
      tol = 1e-08
    )

#### Arguments

- `k`:

  Number of classes.

- `method`:

  `"em"` or `"sem"`.

- `max_iter`:

  Maximum number of iterations.

- `tol`:

  Tolerance for the relative change of the log-likelihood.

------------------------------------------------------------------------

### `EntropySegmenter$fit()`

Fit the mixture.

#### Usage

    EntropySegmenter$fit(x, init = NULL, is_map = TRUE)

#### Arguments

- `x`:

  Entropy features: a matrix (one map, or pixels in rows and features in
  columns when `is_map = FALSE`), an array `nrow x ncol x features`
  (several maps), or a numeric vector.

- `init`:

  Initial labels (integer vector) or `NULL` for k-means.

- `is_map`:

  Treat a matrix `x` as a single map? Default `TRUE`.

------------------------------------------------------------------------

### `EntropySegmenter$predict()`

Labels of new observations (maximum a posteriori).

#### Usage

    EntropySegmenter$predict(x, is_map = TRUE)

#### Arguments

- `x`:

  Entropy features in the same format used in `$fit()`.

- `is_map`:

  Treat a matrix `x` as a single map?

------------------------------------------------------------------------

### `EntropySegmenter$label_map()`

Labels as a matrix with the spatial dimensions of the fitted map.

#### Usage

    EntropySegmenter$label_map()

------------------------------------------------------------------------

### `EntropySegmenter$print()`

Print a summary.

#### Usage

    EntropySegmenter$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `EntropySegmenter$clone()`

The objects of this class are cloneable with this method.

#### Usage

    EntropySegmenter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
H <- entropy_map(SanFrancisco150, model = "wishart", window = 7, L = 4)
set.seed(1)
seg <- EntropySegmenter$new(k = 3)$fit(H)
seg
#> <EntropySegmenter> k = 3 | method: EM 
#>   iterations: 20 | converged: TRUE 
#>   weights: 0.265 0.268 0.467 
#>   class sizes: 5933 5773 10794 
plot_map(seg$label_map(), main = "Entropy-based segmentation")

kappa_statistics(SanFrancisco150_gt, seg$label_map())
#> # A tibble: 1 × 5
#>   kappa   variance std_error overall_accuracy     n
#>   <dbl>      <dbl>     <dbl>            <dbl> <dbl>
#> 1 0.914 0.00000566   0.00238            0.945 22500
```

# Kappa coefficient and overall accuracy

Computes Cohen's kappa coefficient of agreement, its asymptotic variance
(Congalton and Green, 2009) and the overall accuracy from a confusion
matrix or from two vectors (or matrices) of labels. Unlike
[`yardstick::kap()`](https://yardstick.tidymodels.org/reference/kap.html),
the variance of the estimator is also returned, which allows comparing
classifiers with Z tests.

## Usage

``` r
kappa_statistics(truth, estimate = NULL, match = FALSE)
```

## Arguments

- truth:

  Either a square confusion matrix (estimates in rows, truth in columns)
  or a vector/matrix of reference labels.

- estimate:

  Vector/matrix of estimated labels (when `truth` is not a confusion
  matrix).

- match:

  Relabel the estimated classes to maximise agreement?

## Value

A tibble with `kappa`, `variance`, `std_error`, `overall_accuracy` and
the number of observations `n`.

## Details

Labels of unsupervised methods are arbitrary: use `match = TRUE` to
relabel the estimated classes so that they best agree with the truth (an
exhaustive search over permutations for up to 8 classes).

## References

Congalton, R. G. and Green, K. (2009). *Assessing the Accuracy of
Remotely Sensed Data: Principles and Practices*. CRC Press.

## Examples

``` r
truth <- c(1, 1, 2, 2, 3, 3, 3)
est <- c(1, 2, 2, 2, 3, 3, 1)
kappa_statistics(truth, est)
#> # A tibble: 1 × 5
#>   kappa variance std_error overall_accuracy     n
#>   <dbl>    <dbl>     <dbl>            <dbl> <dbl>
#> 1 0.576   0.0583     0.242            0.714     7
kappa_statistics(table(est, truth))
#> # A tibble: 1 × 5
#>   kappa variance std_error overall_accuracy     n
#>   <dbl>    <dbl>     <dbl>            <dbl> <dbl>
#> 1 0.576   0.0583     0.242            0.714     7
```

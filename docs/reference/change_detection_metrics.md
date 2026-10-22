# Accuracy measures for change detection

Computes the usual accuracy measures of binary change maps: detection
rate (DR, percentage of true changes that were detected), false alarm
rate (FA, percentage of unchanged pixels flagged as changed), and the
percentages of false positives (FP) and false negatives (FN) with
respect to all pixels.

## Usage

``` r
change_detection_metrics(truth, estimate)
```

## Arguments

- truth:

  Logical (or 0/1) vector/matrix of true changes.

- estimate:

  Logical (or 0/1) vector/matrix of detected changes.

## Value

A tibble with `DR`, `FA`, `FP`, `FN` (in percent) and the overall
accuracy `OA` (in percent).

## Examples

``` r
truth <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
est <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
change_detection_metrics(truth, est)
#> # A tibble: 1 × 5
#>      DR    FA    FP    FN    OA
#>   <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1    50  33.3    20    20    60
```

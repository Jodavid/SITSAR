# Ground truth of the San Francisco image

Reference map of the
[SanFrancisco150](https://jodavid.github.io/SITSAR/reference/SanFrancisco150.md)
image with three classes.

## Usage

``` r
SanFrancisco150_gt
```

## Format

An integer matrix `150 x 150` with values `1` (ocean), `2` (vegetation)
and `3` (urban area).

## Examples

``` r
table(SanFrancisco150_gt)
#> SanFrancisco150_gt
#>     1     2     3 
#>  5927  6433 10140 
plot_map(SanFrancisco150_gt)
```

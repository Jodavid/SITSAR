# Example PolSAR covariance matrices

Typical \\3 \times 3\\ covariance matrices (scaled to unit magnitude)
for urban, forest and pasture areas, taken from the simulation studies
of the authors' dissertation and thesis. They are handy to simulate
PolSAR data.

## Usage

``` r
polsar_example_sigma(type = c("urban", "forest", "pasture"))
```

## Arguments

- type:

  One of `"urban"`, `"forest"` or `"pasture"`.

## Value

A complex Hermitian positive definite \\3 \times 3\\ matrix.

## Examples

``` r
polsar_example_sigma("forest")
#>                    [,1]               [,2]               [,3]
#> [1,] 0.360932+0.000000i 0.011050+0.003759i 0.063896+0.001581i
#> [2,] 0.011050-0.003759i 0.098960+0.000000i 0.006593+0.006868i
#> [3,] 0.063896-0.001581i 0.006593-0.006868i 0.208843+0.000000i
```

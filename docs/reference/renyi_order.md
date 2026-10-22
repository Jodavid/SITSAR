# Tuning parameters for the divergence classifier

'dials' parameter objects for
[`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md):
`renyi_order()` for the order of the Rényi divergence and
`divergence_type()` for the divergence measure. Both require the 'dials'
package.

## Usage

``` r
renyi_order(range = c(0.05, 0.95), trans = NULL)

divergence_type(values = c("renyi", "kl", "jeffreys", "bhattacharyya"))

renyi_kind(range = c(1L, 4L))
```

## Arguments

- range:

  Range of the Rényi order.

- trans:

  A transformation object (or `NULL`).

- values:

  Candidate divergences.

## Value

A 'dials' parameter object.

## Examples

``` r
if (requireNamespace("dials", quietly = TRUE)) {
  renyi_order()
  divergence_type()
}
#> Divergence (qualitative)
#> 4 possible values include:
#> 'renyi', 'kl', 'jeffreys', and 'bhattacharyya'
```

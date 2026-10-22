# Entropy-based segmentation of a SAR/PolSAR image

Convenience wrapper that computes the entropy map(s) with
[`entropy_map()`](https://jodavid.github.io/SITSAR/reference/entropy_map.md)
and segments them with
[EntropySegmenter](https://jodavid.github.io/SITSAR/reference/EntropySegmenter.md).

## Usage

``` r
segment_sar(
  image,
  k = 2,
  model = c("wishart", "g0p", "gamma", "g0i", "ggamma"),
  window = 7,
  L = NULL,
  channels = NULL,
  method = c("em", "sem"),
  ...
)
```

## Arguments

- image:

  A [SARImage](https://jodavid.github.io/SITSAR/reference/SARImage.md)
  or an array accepted by `SARImage$new()`.

- k:

  Number of classes.

- model:

  One of `"gamma"`, `"g0i"`, `"ggamma"` (intensity models, applied to
  each channel) or `"wishart"`, `"g0p"` (matrix models, for PolSAR
  images).

- window:

  Odd side of the moving window.

- L:

  Known number of looks. If `NULL` (default) it is estimated in every
  window.

- channels:

  Channels used by intensity models (default: all).

- method:

  `"em"` or `"sem"`.

- ...:

  Passed to `EntropySegmenter$new()`.

## Value

A fitted
[EntropySegmenter](https://jodavid.github.io/SITSAR/reference/EntropySegmenter.md);
the labels are available with `$label_map()`.

## Examples

``` r
set.seed(1)
seg <- segment_sar(SanFrancisco150, k = 3, model = "gamma", window = 7, L = 4)
table(seg$label_map(), SanFrancisco150_gt)
#>    SanFrancisco150_gt
#>        1    2    3
#>   1 5845   72    0
#>   2   81 5904  350
#>   3    1  457 9790
```

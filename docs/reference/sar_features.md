# Local features of SAR images for pixel classification

Computes moving-window statistics around every pixel, in the format used
by
[DivergenceClassifier](https://jodavid.github.io/SITSAR/reference/DivergenceClassifier.md)
and the
[`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md)
'parsnip' model. Window statistics are computed with integral images, so
the cost does not depend on the window size. Windows are truncated at
the image borders.

## Usage

``` r
sar_features(image, window = 3, type = c("all", "gaussian", "wishart"))
```

## Arguments

- image:

  A [SARImage](https://jodavid.github.io/SITSAR/reference/SARImage.md)
  or an array accepted by `SARImage$new()`.

- window:

  Odd window side.

- type:

  Which features: `"all"`, `"gaussian"` or `"wishart"`.

## Value

A tibble with the columns `row`, `col` and the features.

## Details

Two groups of columns are produced:

- Gaussian features:

  `mean_<ch>`: local mean of each intensity channel; `cov_<ch1>_<ch2>`:
  local (unbiased) covariance between the intensities. Together they
  define the multivariate normal fitted to the intensity vector in the
  window.

- Wishart features (PolSAR only):

  `c<i><i>` and `c<i><j>_re`/`c<i><j>_im`: entries of the local mean
  covariance matrix, i.e. the maximum likelihood estimate of \\\Sigma\\
  for the scaled complex Wishart law in the window.

## Examples

``` r
img <- SARImage$new(SanFrancisco150)$crop(1:20, 1:20)
f <- sar_features(img, window = 3)
names(f)
#>  [1] "row"       "col"       "mean_HH"   "mean_HV"   "mean_VV"   "cov_HH_HH"
#>  [7] "cov_HH_HV" "cov_HH_VV" "cov_HV_HV" "cov_HV_VV" "cov_VV_VV" "c11"      
#> [13] "c22"       "c33"       "c12_re"    "c12_im"    "c13_re"    "c13_im"   
#> [19] "c23_re"    "c23_im"   
```

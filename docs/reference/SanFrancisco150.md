# San Francisco PolSAR image (150 x 150 pixels)

A \\150 \times 150\\ crop of the fully polarimetric L-band AIRSAR image
of the San Francisco Bay area (nominal 4 looks), used in the
dissertation and thesis of J. A. Ferreira. The scene contains ocean,
vegetation (forest/park) and urban areas; see
[SanFrancisco150_gt](https://jodavid.github.io/SITSAR/reference/SanFrancisco150_gt.md).

## Usage

``` r
SanFrancisco150
```

## Format

A complex array of dimension `150 x 150 x 6`. For each pixel the six
layers are the distinct entries of the \\3 \times 3\\ Hermitian
covariance matrix: \\\|S\_{HH}\|^2\\, \\\|S\_{HV}\|^2\\,
\\\|S\_{VV}\|^2\\ (real intensities), \\S\_{HH}S\_{HV}^\*\\,
\\S\_{HH}S\_{VV}^\*\\ and \\S\_{HV}S\_{VV}^\*\\ (complex products).

## Source

NASA/JPL AIRSAR mission.

## Examples

``` r
dim(SanFrancisco150)
#> [1] 150 150   6
img <- SARImage$new(SanFrancisco150, looks = 4)
img$plot()
```

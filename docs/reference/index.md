# Package index

## Package overview

- [`SITSAR`](https://jodavid.github.io/SITSAR/reference/SITSAR-package.md)
  [`SITSAR-package`](https://jodavid.github.io/SITSAR/reference/SITSAR-package.md)
  : SITSAR: Information Theory, Geometry and Statistics for SAR Data

## Intensity distributions

Univariate models for SAR intensities.

- [`dg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md)
  [`pg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md)
  [`qg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md)
  [`rg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md) : The
  G0I distribution for SAR intensities
- [`dgammasar()`](https://jodavid.github.io/SITSAR/reference/dgammasar.md)
  [`pgammasar()`](https://jodavid.github.io/SITSAR/reference/dgammasar.md)
  [`qgammasar()`](https://jodavid.github.io/SITSAR/reference/dgammasar.md)
  [`rgammasar()`](https://jodavid.github.io/SITSAR/reference/dgammasar.md)
  : The Gamma distribution for SAR intensities
- [`dggamma()`](https://jodavid.github.io/SITSAR/reference/dggamma.md)
  [`pggamma()`](https://jodavid.github.io/SITSAR/reference/dggamma.md)
  [`qggamma()`](https://jodavid.github.io/SITSAR/reference/dggamma.md)
  [`rggamma()`](https://jodavid.github.io/SITSAR/reference/dggamma.md) :
  The generalized Gamma distribution for SAR intensities

## Matrix (PolSAR) distributions

Models for Hermitian positive definite PolSAR matrices.

- [`dwishart_sar()`](https://jodavid.github.io/SITSAR/reference/dwishart_sar.md)
  [`rwishart_sar()`](https://jodavid.github.io/SITSAR/reference/dwishart_sar.md)
  : The scaled complex Wishart distribution
- [`dg0p()`](https://jodavid.github.io/SITSAR/reference/dg0p.md)
  [`rg0p()`](https://jodavid.github.io/SITSAR/reference/dg0p.md) : The
  G0 polarimetric distribution
- [`dctpcw()`](https://jodavid.github.io/SITSAR/reference/compound_wishart.md)
  [`rctpcw()`](https://jodavid.github.io/SITSAR/reference/compound_wishart.md)
  [`dcgcw()`](https://jodavid.github.io/SITSAR/reference/compound_wishart.md)
  [`rcgcw()`](https://jodavid.github.io/SITSAR/reference/compound_wishart.md)
  : Compound complex Wishart distributions
- [`as_polsar_array()`](https://jodavid.github.io/SITSAR/reference/polsar-format.md)
  [`as_polsar_vec()`](https://jodavid.github.io/SITSAR/reference/polsar-format.md)
  : Convert between PolSAR storage formats
- [`polsar_example_sigma()`](https://jodavid.github.io/SITSAR/reference/polsar_example_sigma.md)
  : Example PolSAR covariance matrices

## R6 statistical models

Estimation, Shannon entropy and asymptotic variances.

- [`SARModel`](https://jodavid.github.io/SITSAR/reference/SARModel.md) :
  Base R6 class for SAR and PolSAR statistical models
- [`GammaSAR`](https://jodavid.github.io/SITSAR/reference/GammaSAR.md) :
  Gamma model for SAR intensities (R6)
- [`G0I`](https://jodavid.github.io/SITSAR/reference/G0I.md) : G0I model
  for SAR intensities (R6)
- [`GGammaSAR`](https://jodavid.github.io/SITSAR/reference/GGammaSAR.md)
  : Generalized Gamma model for SAR intensities (R6)
- [`WishartSAR`](https://jodavid.github.io/SITSAR/reference/WishartSAR.md)
  : Scaled complex Wishart model for PolSAR data (R6)
- [`G0Polarimetric`](https://jodavid.github.io/SITSAR/reference/G0Polarimetric.md)
  : G0 polarimetric model for PolSAR data (R6)
- [`CompoundWishart`](https://jodavid.github.io/SITSAR/reference/CompoundWishart.md)
  : Compound complex Wishart models (R6)
- [`tidy(`*`<SARModel>`*`)`](https://jodavid.github.io/SITSAR/reference/sarmodel-methods.md)
  [`glance(`*`<SARModel>`*`)`](https://jodavid.github.io/SITSAR/reference/sarmodel-methods.md)
  [`logLik(`*`<SARModel>`*`)`](https://jodavid.github.io/SITSAR/reference/sarmodel-methods.md)
  : Tidy and glance methods for SITSAR models
- [`g0iestimation()`](https://jodavid.github.io/SITSAR/reference/g0iestimation.md)
  : Estimation of G0I parameters

## Shannon entropy

- [`entropy_gamma()`](https://jodavid.github.io/SITSAR/reference/entropy.md)
  [`entropy_g0i()`](https://jodavid.github.io/SITSAR/reference/entropy.md)
  [`entropy_ggamma()`](https://jodavid.github.io/SITSAR/reference/entropy.md)
  [`entropy_wishart()`](https://jodavid.github.io/SITSAR/reference/entropy.md)
  [`entropy_g0p()`](https://jodavid.github.io/SITSAR/reference/entropy.md)
  : Shannon entropies of SAR and PolSAR models
- [`entropy_map()`](https://jodavid.github.io/SITSAR/reference/entropy_map.md)
  : Shannon entropy maps
- [`entropy_test()`](https://jodavid.github.io/SITSAR/reference/entropy_test.md)
  : Entropy-based test for two SAR samples

## Images

- [`SARImage`](https://jodavid.github.io/SITSAR/reference/SARImage.md) :
  SAR and PolSAR images (R6)
- [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md)
  : Local features of SAR images for pixel classification
- [`plot_map()`](https://jodavid.github.io/SITSAR/reference/plot_map.md)
  : Plot a map

## Segmentation and change detection

- [`EntropySegmenter`](https://jodavid.github.io/SITSAR/reference/EntropySegmenter.md)
  : Entropy-based segmentation by (stochastic) EM (R6)
- [`segment_sar()`](https://jodavid.github.io/SITSAR/reference/segment_sar.md)
  : Entropy-based segmentation of a SAR/PolSAR image
- [`change_map()`](https://jodavid.github.io/SITSAR/reference/change_map.md)
  : Entropy-based change detection between two images

## Divergences and classification

- [`divergence_gaussian()`](https://jodavid.github.io/SITSAR/reference/divergences.md)
  [`divergence_wishart()`](https://jodavid.github.io/SITSAR/reference/divergences.md)
  : Divergences between Gaussian and complex Wishart models
- [`distancesg0i()`](https://jodavid.github.io/SITSAR/reference/distancesg0i.md)
  : Stochastic distances between G0I distributions
- [`distancesg0itable()`](https://jodavid.github.io/SITSAR/reference/distancesg0itable.md)
  : Table of G0I stochastic distances between two image regions
- [`DivergenceClassifier`](https://jodavid.github.io/SITSAR/reference/DivergenceClassifier.md)
  : Divergence-based classifier for SAR/PolSAR features (R6)
- [`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md)
  [`update(`*`<divergence_classifier>`*`)`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md)
  : Divergence-based classifier as a 'parsnip' model
- [`renyi_order()`](https://jodavid.github.io/SITSAR/reference/renyi_order.md)
  [`divergence_type()`](https://jodavid.github.io/SITSAR/reference/renyi_order.md)
  [`renyi_kind()`](https://jodavid.github.io/SITSAR/reference/renyi_order.md)
  : Tuning parameters for the divergence classifier

## Accuracy assessment

- [`kappa_statistics()`](https://jodavid.github.io/SITSAR/reference/kappa_statistics.md)
  : Kappa coefficient and overall accuracy
- [`change_detection_metrics()`](https://jodavid.github.io/SITSAR/reference/change_detection_metrics.md)
  : Accuracy measures for change detection

## Datasets

- [`SanFrancisco150`](https://jodavid.github.io/SITSAR/reference/SanFrancisco150.md)
  : San Francisco PolSAR image (150 x 150 pixels)
- [`SanFrancisco150_gt`](https://jodavid.github.io/SITSAR/reference/SanFrancisco150_gt.md)
  : Ground truth of the San Francisco image

## Re-exported functions

- [`reexports`](https://jodavid.github.io/SITSAR/reference/reexports.md)
  [`as_tibble`](https://jodavid.github.io/SITSAR/reference/reexports.md)
  [`tidy`](https://jodavid.github.io/SITSAR/reference/reexports.md)
  [`glance`](https://jodavid.github.io/SITSAR/reference/reexports.md) :
  Objects exported from other packages

# Changelog

## SITSAR 1.0.0

- Complete rewrite of the package, gathering the methods of the master’s
  dissertation and PhD thesis of Jodavid A. Ferreira.
- Intensity distributions: Gamma, generalized Gamma and G0I (exact
  distribution and quantile functions through the F distribution).
- Matrix distributions: scaled complex Wishart, G0 polarimetric,
  compound truncated Poisson and compound geometric complex Wishart.
- R6 classes `GammaSAR`, `G0I`, `GGammaSAR`, `WishartSAR`,
  `G0Polarimetric` and `CompoundWishart` with maximum likelihood / EM
  estimation, Shannon entropy and asymptotic variance of the entropy
  estimator, plus
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) methods.
- `SARImage` R6 class,
  [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md),
  fast vectorised
  [`entropy_map()`](https://jodavid.github.io/SITSAR/reference/entropy_map.md)
  and
  [`plot_map()`](https://jodavid.github.io/SITSAR/reference/plot_map.md).
- Entropy-based segmentation (`EntropySegmenter`,
  [`segment_sar()`](https://jodavid.github.io/SITSAR/reference/segment_sar.md)),
  hypothesis tests
  ([`entropy_test()`](https://jodavid.github.io/SITSAR/reference/entropy_test.md))
  and change detection
  ([`change_map()`](https://jodavid.github.io/SITSAR/reference/change_map.md)).
- Divergences between Gaussian and complex Wishart laws, the
  `DivergenceClassifier` R6 class and the
  [`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md)
  ‘parsnip’ model for ‘tidymodels’ workflows.
- Accuracy measures:
  [`kappa_statistics()`](https://jodavid.github.io/SITSAR/reference/kappa_statistics.md)
  and
  [`change_detection_metrics()`](https://jodavid.github.io/SITSAR/reference/change_detection_metrics.md).
- New dataset `SanFrancisco150_gt` (ground truth) and new octagonal
  logo.
- Functions of previous versions
  ([`dg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md),
  [`pg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md),
  [`qg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md),
  [`rg0i()`](https://jodavid.github.io/SITSAR/reference/dg0i.md),
  [`g0iestimation()`](https://jodavid.github.io/SITSAR/reference/g0iestimation.md),
  [`distancesg0i()`](https://jodavid.github.io/SITSAR/reference/distancesg0i.md)
  and
  [`distancesg0itable()`](https://jodavid.github.io/SITSAR/reference/distancesg0itable.md))
  are kept.

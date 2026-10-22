# SITSAR: Information Theory, Geometry and Statistics for SAR Data

Statistical models, information-theoretic measures and image analysis
tools for Synthetic Aperture Radar (SAR) and polarimetric SAR (PolSAR)
data, gathering the methods developed in the master's dissertation
(divergence-based classifiers) and the PhD thesis (Shannon entropies,
entropy-based segmentation and change detection, and compound Wishart
models) of Jodavid A. Ferreira, supervised by Abraão D. C. Nascimento.

## Main tools

- Distributions:

  [dg0i](https://jodavid.github.io/SITSAR/reference/dg0i.md),
  [dgammasar](https://jodavid.github.io/SITSAR/reference/dgammasar.md),
  [dggamma](https://jodavid.github.io/SITSAR/reference/dggamma.md),
  [dwishart_sar](https://jodavid.github.io/SITSAR/reference/dwishart_sar.md),
  [dg0p](https://jodavid.github.io/SITSAR/reference/dg0p.md),
  [compound_wishart](https://jodavid.github.io/SITSAR/reference/compound_wishart.md).

- R6 models:

  [GammaSAR](https://jodavid.github.io/SITSAR/reference/GammaSAR.md),
  [G0I](https://jodavid.github.io/SITSAR/reference/G0I.md),
  [GGammaSAR](https://jodavid.github.io/SITSAR/reference/GGammaSAR.md),
  [WishartSAR](https://jodavid.github.io/SITSAR/reference/WishartSAR.md),
  [G0Polarimetric](https://jodavid.github.io/SITSAR/reference/G0Polarimetric.md),
  [CompoundWishart](https://jodavid.github.io/SITSAR/reference/CompoundWishart.md)
  (estimation, entropy and asymptotic variance).

- Images:

  [SARImage](https://jodavid.github.io/SITSAR/reference/SARImage.md),
  [`sar_features()`](https://jodavid.github.io/SITSAR/reference/sar_features.md),
  [`entropy_map()`](https://jodavid.github.io/SITSAR/reference/entropy_map.md),
  [`plot_map()`](https://jodavid.github.io/SITSAR/reference/plot_map.md).

- Segmentation and change detection:

  [EntropySegmenter](https://jodavid.github.io/SITSAR/reference/EntropySegmenter.md),
  [`segment_sar()`](https://jodavid.github.io/SITSAR/reference/segment_sar.md),
  [`entropy_test()`](https://jodavid.github.io/SITSAR/reference/entropy_test.md),
  [`change_map()`](https://jodavid.github.io/SITSAR/reference/change_map.md).

- Classification:

  [divergences](https://jodavid.github.io/SITSAR/reference/divergences.md),
  [DivergenceClassifier](https://jodavid.github.io/SITSAR/reference/DivergenceClassifier.md)
  and the 'parsnip' model
  [`divergence_classifier()`](https://jodavid.github.io/SITSAR/reference/divergence_classifier.md).

- Accuracy:

  [`kappa_statistics()`](https://jodavid.github.io/SITSAR/reference/kappa_statistics.md),
  [`change_detection_metrics()`](https://jodavid.github.io/SITSAR/reference/change_detection_metrics.md).

## See also

Useful links:

- <https://github.com/jodavid/SITSAR>

- <https://jodavid.github.io/SITSAR/>

- Report bugs at <https://github.com/jodavid/SITSAR/issues>

## Author

**Maintainer**: Jodavid Ferreira <jodavid@protonmail.com>
([ORCID](https://orcid.org/0000-0002-2131-6464))

Authors:

- Jodavid Ferreira <jodavid@protonmail.com>
  ([ORCID](https://orcid.org/0000-0002-2131-6464))

- Abraao Nascimento <abraao@de.ufpe.br>
  ([ORCID](https://orcid.org/0000-0003-2673-219X))

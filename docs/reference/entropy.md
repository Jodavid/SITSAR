# Shannon entropies of SAR and PolSAR models

Closed-form Shannon entropies \\H = E\\-\log f(Z)\\\\ for the models
used in SAR (intensity) and PolSAR (matrix) image analysis. All
functions are vectorised in the scalar parameters.

## Usage

``` r
entropy_gamma(mu, L)

entropy_g0i(alpha, gamma, L)

entropy_ggamma(alpha, phi, mu)

entropy_wishart(Sigma, L)

entropy_g0p(alpha, Sigma, L)
```

## Arguments

- mu:

  Mean of the Gamma law, or inverse scale of the generalized Gamma law.

- L:

  Number of looks.

- alpha:

  Roughness (G0I, G0 polarimetric) or power (generalized Gamma)
  parameter.

- gamma:

  Scale parameter of the G0I law.

- phi:

  Shape parameter of the generalized Gamma law.

- Sigma:

  Hermitian positive definite covariance matrix.

## Value

Numeric vector of entropies.

## Details

- Gamma:

  \\H = L - \log L + \log\mu + \log\Gamma(L) + (1 - L)\psi^{(0)}(L)\\.

- G0I:

  (Ferreira and Nascimento, 2020) \\H =
  -\log\frac{L^L\Gamma(L-\alpha)}{\Gamma(-\alpha)\gamma^\alpha
  \Gamma(L)} + (1 - L)\[\psi^{(0)}(L) - \log L + \log\gamma -
  \psi^{(0)}(-\alpha)\] + (L - \alpha)\[\psi^{(0)}(L - \alpha) +
  \log\gamma - \psi^{(0)}(-\alpha)\]\\.

- Generalized Gamma:

  \\H = -\log\alpha + \log\Gamma(\phi) + \phi - (\phi -
  1/\alpha)\psi^{(0)}(\phi) - \log\mu\\.

- Scaled complex Wishart:

  (Nascimento, Ferreira and Frery, 2023) \\H = m\log\|\Sigma\| - m^2\log
  L + \log\Gamma_m(L) + mL + (m - L)\psi^{(0)}\_m(L)\\, with
  \\\psi^{(0)}\_m(L) = \sum\_{i=0}^{m-1}\psi^{(0)}(L - i)\\.

- G0 polarimetric:

  (Nascimento, Ferreira and Frery, 2023) obtained from \\E\log\|Z\| =
  m\[\log(-\alpha-1) - \psi^{(0)}(-\alpha)\] + \log\|\Sigma\| - m\log
  L + \psi^{(0)}\_m(L)\\ and \\E\log\[L\\\mathrm{tr}(\Sigma^{-1}Z) -
  \alpha - 1\] = \log(-\alpha-1) + \psi^{(0)}(mL - \alpha) -
  \psi^{(0)}(-\alpha)\\.

## References

Ferreira, J. A. and Nascimento, A. D. C. (2020). Shannon entropy for the
G0I model: a new segmentation approach. *IEEE Journal of Selected Topics
in Applied Earth Observations and Remote Sensing*, 13, 2547–2553.
[doi:10.1109/JSTARS.2020.2997666](https://doi.org/10.1109/JSTARS.2020.2997666)

Nascimento, A. D. C., Ferreira, J. A. and Frery, A. C. (2023).
Unsupervised segmentation of PolSAR data with complex Wishart and G0m
distributions and Shannon entropy. *Statistics and Computing*, 33, 129.
[doi:10.1007/s11222-023-10298-0](https://doi.org/10.1007/s11222-023-10298-0)

## Examples

``` r
entropy_gamma(mu = 1, L = 4)
#> [1] 0.6371121
entropy_g0i(alpha = -3, gamma = 2, L = 4)
#> [1] 0.8625083
entropy_ggamma(alpha = 1, phi = 4, mu = 4) # equals the Gamma entropy
#> [1] 0.6371121
S <- polsar_example_sigma("forest")
entropy_wishart(S, L = 4)
#> [1] -12.04464
entropy_g0p(alpha = -5, Sigma = S, L = 4)
#> [1] -12.48772
```

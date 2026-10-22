# Divergences between Gaussian and complex Wishart models

Closed-form divergences between two multivariate normal laws
\\N_p(\mu_1, \Sigma_1)\\ and \\N_p(\mu_2, \Sigma_2)\\, and between two
scaled complex Wishart laws \\\mathcal{W}\_m(\Sigma_1, L)\\ and
\\\mathcal{W}\_m(\Sigma_2, L)\\ with a common number of looks. They are
the building blocks of the divergence-based classifiers of Ferreira,
Coelho and Nascimento (2021).

## Usage

``` r
divergence_gaussian(
  mu1,
  Sigma1,
  mu2,
  Sigma2,
  type = c("kl", "jeffreys", "bhattacharyya", "renyi"),
  order = 0.5,
  kind = 1
)

divergence_wishart(
  Sigma1,
  Sigma2,
  L,
  type = c("kl", "jeffreys", "bhattacharyya", "renyi"),
  order = 0.5,
  kind = 1
)
```

## Arguments

- mu1, mu2:

  Mean vectors.

- Sigma1, Sigma2:

  Covariance matrices (real symmetric for Gaussian laws, Hermitian for
  Wishart laws).

- type:

  Divergence: `"kl"`, `"jeffreys"`, `"bhattacharyya"` or `"renyi"`.

- order:

  Order \\\alpha \in (0, 1)\\ of the Rényi divergence.

- kind:

  Form of the Rényi divergence (1 to 4).

- L:

  Number of looks (common to both Wishart laws).

## Value

A single non-negative number.

## Details

Available measures (`type`):

- `"kl"`:

  Kullback-Leibler divergence \\D\_{KL}(\theta_1 \\ \theta_2)\\.

- `"jeffreys"`:

  Jeffreys divergence \\D\_{KL}(\theta_1 \\ \theta_2) + D\_{KL}(\theta_2
  \\ \theta_1)\\.

- `"bhattacharyya"`:

  Bhattacharyya distance \\-\log \int \sqrt{f_1 f_2}\\.

- `"renyi"`:

  Rényi divergence of order \\\alpha\\ (`order`), in one of four forms
  (`kind`): 1, \\D^\alpha_R(\theta_1\\\theta_2)\\; 2,
  \\D^\alpha_R(\theta_2\\\theta_1)\\; 3, the symmetrised mean
  \\\[D^\alpha_R(\theta_1\\\theta_2) +
  D^\alpha_R(\theta_2\\\theta_1)\]/2\\; and 4, \\c
  \log\\\[K(\theta_1\\\theta_2) + K(\theta_2\\\theta_1)\]/2\\\\, where
  \\K(\theta_1\\\theta_2) = \int f_1^\alpha f_2^{1-\alpha}\\.

For Gaussian laws the Rényi divergence uses the normalisation of Pardo
(2006), \\c = 1/\[\alpha(\alpha - 1)\]\\, which gives \$\$D^\alpha_R =
\frac{1}{2}\Delta^\top M^{-1}\Delta -
\frac{1}{2\alpha(\alpha-1)}\log\frac{\|M\|}{\|\Sigma_1\|^{1-\alpha}
\|\Sigma_2\|^{\alpha}},\$\$ with \\\Delta = \mu_1 - \mu_2\\ and \\M =
\alpha\Sigma_2 + (1-\alpha)\Sigma_1\\. For Wishart laws the
normalisation of Frery, Nascimento and Cintra (2014), \\c = 1/(\alpha -
1)\\, is used: \$\$D^\alpha_R = \frac{L}{\alpha -
1}\log\frac{\|(\alpha\Sigma_1^{-1} +
(1-\alpha)\Sigma_2^{-1})^{-1}\|}{\|\Sigma_1\|^\alpha\|\Sigma_2\|^{1-\alpha}}.\$\$
Both converge to the Kullback-Leibler divergence when \\\alpha \to 1\\.

## References

Ferreira, J. A., Coelho, H. and Nascimento, A. D. C. (2021). A family of
divergence-based classifiers for polarimetric synthetic aperture radar
(PolSAR) imagery vector and matrix features. *International Journal of
Remote Sensing*, 42(4), 1201–1229.
[doi:10.1080/01431161.2020.1826060](https://doi.org/10.1080/01431161.2020.1826060)

Frery, A. C., Nascimento, A. D. C. and Cintra, R. J. (2014). Analytic
expressions for stochastic distances between relaxed complex Wishart
distributions. *IEEE Transactions on Geoscience and Remote Sensing*,
52(2), 1213–1226.

Pardo, L. (2006). *Statistical Inference Based on Divergence Measures*.
Chapman & Hall/CRC.

## Examples

``` r
S1 <- diag(3); S2 <- diag(c(2, 1, 0.5))
divergence_gaussian(c(0, 0, 0), S1, c(1, 0, 0), S2, type = "kl")
#> [1] 0.5
divergence_gaussian(c(0, 0, 0), S1, c(1, 0, 0), S2, type = "renyi",
                    order = 0.5, kind = 3)
#> [1] 0.5688994

W1 <- polsar_example_sigma("forest"); W2 <- polsar_example_sigma("urban")
divergence_wishart(W1, W2, L = 4, type = "jeffreys")
#> [1] 11.28573
divergence_wishart(W1, W2, L = 4, type = "renyi", order = 0.9)
#> [1] 3.368488
```

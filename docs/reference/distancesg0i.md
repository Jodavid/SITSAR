# Stochastic distances between G0I distributions

Computes stochastic distances between two \\\mathcal{G}\_I^0\\ laws by
numerical integration (Nascimento, Cintra and Frery, 2010).

## Usage

``` r
distancesg0i(
  par_1,
  par_2,
  p = 1/2,
  type = c("KL", "Renyi", "Bhattacharyya", "Hellinger", "Arithmetic-Geometric",
    "Triangular", "Harmonic-Mean")
)
```

## Arguments

- par_1, par_2:

  Parameter vectors `c(alpha, gamma, L)`.

- p:

  Order of the Rényi distance, in (0, 1).

- type:

  One of `"KL"`, `"Renyi"`, `"Bhattacharyya"`, `"Hellinger"`,
  `"Arithmetic-Geometric"`, `"Triangular"` or `"Harmonic-Mean"`.

## Value

An object of class `distancesg0i` with the parameters, the type and the
`distance_value`.

## References

Nascimento, A. D. C., Cintra, R. J. and Frery, A. C. (2010). Hypothesis
testing in speckled data with stochastic distances. *IEEE Transactions
on Geoscience and Remote Sensing*, 48(1), 373–385.

## Examples

``` r
distancesg0i(c(-1.5, 2, 4), c(-1.5, 2, 4))
#> 
#> G0I stochastic distance
#> 
#> Type:      KL 
#> Distance:  0 
res <- distancesg0i(c(-1.5, 2, 8), c(-3, 2, 4), type = "Hellinger")
res
#> 
#> G0I stochastic distance
#> 
#> Type:      Hellinger 
#> Distance:  0.129977 
```

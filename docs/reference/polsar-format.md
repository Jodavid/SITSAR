# Convert between PolSAR storage formats

PolSAR observations are Hermitian positive definite matrices. `SITSAR`
accepts them either as a complex array of dimension `m x m x n` or as an
`n x p` complex matrix in which each row stores the \\p = m(m+1)/2\\
distinct entries in the order \\C\_{11}, \ldots, C\_{mm}\\ (the
intensities) followed by the upper triangle \\C\_{12}, C\_{13}, \ldots,
C\_{23}, \ldots\\. For \\m = 3\\ this is the usual layer order
`HH, HV, VV, HH.HV*, HH.VV*, HV.VV*` of the
[SanFrancisco150](https://jodavid.github.io/SITSAR/reference/SanFrancisco150.md)
image.

## Usage

``` r
as_polsar_array(x)

as_polsar_vec(x)
```

## Arguments

- x:

  A complex array (`m x m x n`), a single `m x m` matrix, or an `n x p`
  complex matrix (or vector of length `p`).

## Value

`as_polsar_array()` returns a complex array `m x m x n`;
`as_polsar_vec()` returns an `n x p` complex matrix.

## Examples

``` r
S <- diag(3) + 0i
S[1, 2] <- 0.2 + 0.1i
S[2, 1] <- Conj(S[1, 2])
v <- as_polsar_vec(S)
v
#>      [,1] [,2] [,3]     [,4] [,5] [,6]
#> [1,] 1+0i 1+0i 1+0i 0.2+0.1i 0+0i 0+0i
as_polsar_array(v)[, , 1]
#>          [,1]     [,2] [,3]
#> [1,] 1.0+0.0i 0.2+0.1i 0+0i
#> [2,] 0.2-0.1i 1.0+0.0i 0+0i
#> [3,] 0.0+0.0i 0.0+0.0i 1+0i
```

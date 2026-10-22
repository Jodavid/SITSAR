# Table of G0I stochastic distances between two image regions

Fits the \\\mathcal{G}\_I^0\\ law to each intensity channel of two
rectangular regions of an image and computes all the stochastic
distances of
[`distancesg0i()`](https://jodavid.github.io/SITSAR/reference/distancesg0i.md)
between them.

## Usage

``` r
distancesg0itable(image, marks1, marks2, p = 1/2, L = NULL)
```

## Arguments

- image:

  A [SARImage](https://jodavid.github.io/SITSAR/reference/SARImage.md)
  or an array (`nrow x ncol x layers`).

- marks1, marks2:

  Region limits as `c(left, right, up, down)`, i.e. first and last
  column, first and last row.

- p:

  Order of the Rényi distance.

- L:

  Known number of looks; estimated when `NULL`.

## Value

An object of class `distancesg0itable` with the marks and the table of
distances (`distancestable`, channels in rows).

## Examples

``` r
# \donttest{
res <- distancesg0itable(SanFrancisco150, c(10, 30, 10, 30),
                         c(10, 30, 120, 140), L = 4)
res
#> 
#> G0I distances table
#> 
#>           KL    Renyi Bhattacharyya Hellinger Arithmetic-Geometric Triangular
#> HH 11.385316 4.271043     2.1355217 0.8818171             5.047965   1.895534
#> HV 37.876250 8.330567     4.1652833 0.9844747            18.248424   1.993597
#> VV  3.211437 1.179735     0.5898677 0.4455994             1.239628   1.178514
#>    Harmonic-Mean
#> HH     2.9520368
#> HV     5.7440612
#> VV     0.8897873
# }
```

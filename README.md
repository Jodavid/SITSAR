
# SITSAR <img src="man/figures/logo.png" style="float: right" height="139"/>

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/SITSAR)](https://cran.r-project.org/package=SITSAR)
[![CRAN
Download](https://cranlogs.r-pkg.org/badges/grand-total/SITSAR)](https://cran.r-project.org/package=SITSAR)
<!-- badges: end -->

Last update: 31-05-2024

## Information Theory, Geometry and Statistics for Synthetic Aperture Radar (SAR) Data

### Installation

``` r
# Installation
install.packages("devtools")
devtools::install_github("jodavid/SITSAR")
```

### Use

``` r
# import package
library(SITSAR)
```

### Example

``` r

library(SITSAR)

set.seed(1) # determining a seed
alpha <- -1.5
gama <- 2
L <- 3
Xobs<- rg0i(100,alpha, gama, L)
```

``` r

dg0i(Xobs,alpha, gama, L)
```

``` r

library(SITSAR)

set.seed(1) # seed fixed
alpha <- -1.5
gama <- 2
L <- 1
Xobs2<- rg0i(100,alpha,gama,L)

# The parameter L is estimated.
g0iestimation(c(alpha,gama,L), Xobs2)

# The parameter L is fixed.
g0iestimation(c(-1.1,1,1),Xobs2, L_fixed = TRUE)



L <- 3
Xobs2<- rg0i(100,alpha,gama,L)

# The parameter L is estimated.
g0iestimation(c(alpha,gama,L), Xobs2)

# The parameter L is fixed.
g0iestimation(c(-1.1,1,3),Xobs2, L_fixed = TRUE)
```

``` r

set.seed(1) # seed fixed
alpha <- -1.5
gama <- 2
L <- 4

distancesg0i(c(alpha,gama,L),c(alpha,gama,L))
res <- distancesg0i(c(alpha,gama,L),c(alpha,gama,4))
print(res)
```

``` r

library(SITSAR)

set.seed(1) # determining a seed
data("SanFrancisco150")
image <- SanFrancisco150
# left,right, up,down)
marks1 <- c(10,30,50,100)
marks2 <- c(20,50,80,130)

teste <- distancesg0itable(image,marks1, marks2, p = 1/2)
teste
```

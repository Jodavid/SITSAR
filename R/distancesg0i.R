#' G0I Distances
#'
#' \code{distancesg0i} G0I Distances
#'
#'
#' @param par_1 parameters vector for sample 1
#' @param par_2 parameters vector for sample 2
#' @param p     distance Renyi parameter
#' @param type type of distance for calculate
#'
#' @return Distance value between samples.
#'
#' @references
#' \itemize{
#'   \item Nascimento, Abraao DC, Renato J. Cintra, and Alejandro C. Frery. "Hypothesis testing in speckled data with stochastic distances." IEEE Transactions on geoscience and remote sensing 48.1 (2009): 373-385.
#' }
#'
#' @examples
#'
#' set.seed(1) # determining a seed
#'
#' alpha <- -1.5
#' gama <- 2
#' L <- 8
#'
#'
#' distancesg0i(c(alpha,gama,L),c(alpha,gama,L))
#' res <- distancesg0i(c(alpha,gama,L),c(alpha,gama,4))
#' print(res)
#'
#'
#' @importFrom stats integrate
#'
#' @export
distancesg0i = function(par_1, par_2, p = 1/2, type = "KL"){
  UseMethod("distancesg0i")
}

#' @export
distancesg0i.default = function(par_1, par_2, p = 1/2, type = "KL"){

  # --------------------------------------------------------
  if(par_1[1] > 0 | par_2[1] > 0 | par_1[2] < 0 | par_2[2] < 0 |
     par_1[2] < 0 | par_2[2] < 0){
    stop("Your need to redefine the parameters")
    }
  # --------------------------------------------------------

  # --------------------------------------------------------
  alpha1 <- par_1[1]
  gama1 <- par_1[2]
  L1 <- par_1[3]
  alpha2 <- par_2[1]
  gama2 <- par_2[2]
  L2 <- par_2[3]
  # --------------------------------------------------------

  # --------------------------------------------------------
  distance <- switch(type,
                # ------------
                "KL" = {
                  # ------------
                  integrate(
                    function(x)
                      ( dg0i(x,alpha1,gama1,L1)-dg0i(x,alpha2,gama2,L2) ) *
                      log( dg0i(x,alpha1,gama1,L1)/dg0i(x,alpha2,gama2,L2) ),
                    0, Inf )$value
                  # ------------
                },
                "Renyi" = {
                  # ------------
                  ( log(integrate(
                    function(x)
                      ( dg0i(x,alpha1,gama1,L1)^p ) * ( dg0i(x,alpha2,gama2,L2)^(1-p) ),
                    0, Inf )$value) +
                      log(integrate(
                        function(x)
                          ( dg0i(x,alpha2,gama2,L2)^p ) * ( dg0i(x,alpha1,gama1,L1)^(1-p) ),
                        0, Inf )$value) )/(2*(p-1))
                  # ------------
                },
                "Bhattacharyya" = {
                  # ------------
                  I <- integrate(
                    function(x)
                      sqrt( dg0i(x,alpha1,gama1,L1) * dg0i(x,alpha2,gama2,L2) ),
                    0, Inf )$value
                  -log(I)
                  # ------------
                },
                "Hellinger" = {
                  # ------------
                  I <- integrate(
                    function(x)
                      sqrt( dg0i(x,alpha1,gama1,L1) * dg0i(x,alpha2,gama2,L2) ),
                    0, Inf )$value
                  1 - I
                  # ------------
                },
                "Arithmetic-Geometric" = {
                  # ------------
                  integrate(
                    function(x)
                      ( dg0i(x,alpha1,gama1,L1)+dg0i(x,alpha2,gama2,L2) ) *
                      log( ( dg0i(x,alpha1,gama1,L1)+dg0i(x,alpha2,gama2,L2) )/
                             ( 2 * sqrt(dg0i(x,alpha1,gama1,L1)*dg0i(x,alpha2,gama2,L2)) )
                      ),
                    0, Inf )$value
                  # ------------
                },
                "Triangular" = {
                  # ------------
                  integrate(
                    function(x)
                      ( dg0i(x,alpha1,gama1,L1) - dg0i(x,alpha2,gama2,L2) )^2/
                      ( dg0i(x,alpha1,gama1,L1) + dg0i(x,alpha2,gama2,L2) ),
                    0, Inf )$value
                  # ------------
                },
                "Harmonic-Mean" = {
                  # ------------
                  I <- integrate(
                    function(x)
                      ( dg0i(x,alpha1,gama1,L1) - dg0i(x,alpha2,gama2,L2) )^2/
                      ( dg0i(x,alpha1,gama1,L1) + dg0i(x,alpha2,gama2,L2) ),
                    0, Inf )$value
                  -log(1-I/2)
                  # ------------
                },
  )
  # --------------------------------------------------------

  #-------------------------------------------------------
  structure(list(
    par_1 = par_1,
    par_2 = par_2,
    type = type,
    distance_value = distance
  ),
  class = "distancesg0i"
  )

}


#' @export
print.distancesg0i <- function(x, ...) {
  # -----------------
  cat("\nDistance Calculation Summary\n\n")
  # -----------------
  cat("Type:\n")
  print(x$type)
  cat("Distance:\n")
  print(x$distance_value)
  # -----------------
}

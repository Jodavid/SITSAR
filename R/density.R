#' The G0I Distribution
#'
#' \code{dg0i} function related to the G0I Distribution.
#'
#'
#' @param x vector of quantiles.
#' @param alpha value for alpha parameter, this value is < 0.
#' @param gama value for gamma parameter, this value is > 0.
#' @param L value for L parameter, this value is > 0.
#'
#' @return \code{dg0i} gives the density
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
#' Xobs<- rg0i(100,alpha, gama, L)
#' dg0i(Xobs,alpha, gama, L)
#'
#'
#' @importFrom stats rgamma
#' @importFrom invgamma rinvgamma
#'
#' @export
dg0i = function(x,alpha, gama, L){

  # ..> Parameters
  #alpha = theta[1]; gama = theta[2]; L = theta[3];
  # ..> Normalizing constant
  const <- L^L/(beta(-alpha,L)*(gama^alpha))
  # ..> Density
  res <- const * x^(L-1)/(gama+L*x)^(L-alpha)
  return(res)

}

#' The G0I Distribution
#'
#' \code{rg0i} function related to the G0I Distribution.
#'
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param alpha value for alpha parameter, this value is < 0.
#' @param gama value for gamma parameter, this value is > 0.
#' @param L value for L parameter, this value is > 0.
#'
#' @return \code{rg0i} generates random deviates.
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
#' Xobs<- rg0i(100,alpha, gama, L)
#' dg0i(Xobs,alpha, gama, L)
#'
#'
#' @importFrom stats rgamma
#' @importFrom invgamma rinvgamma
#'
#' @export
rg0i = function(n,alpha,gama,L){
  res = stats::rgamma(n,shape=L,rate=L) * invgamma::rinvgamma(n, shape=-alpha, rate = gama)
  return(res)
}

# ----------------

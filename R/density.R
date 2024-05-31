#' The G0I Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the g0i distribution with parameters alpha, gamma and L.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param alpha value for alpha parameter, this value is < 0.
#' @param gamma value for gamma parameter, this value is > 0.
#' @param L value for L parameter, this value is > 0.
#' @param log.p Logical; if TRUE, probabilities p are given as \eqn{log(p)}.
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X ≤ x]} otherwise, \eqn{P[X > x]}.
#'
#' @return The density, distribution function, quantile function or random generation for the G0I distribution.
#'
#' @name g0i
#'
#' @aliases dg0i pg0i qg0i rg0i
#'
#' @references
#' \itemize{
#'   \item Nascimento, Abraao DC, Renato J. Cintra, and Alejandro C. Frery. "Hypothesis testing in speckled data with stochastic distances." IEEE Transactions on geoscience and remote sensing 48.1 (2009): 373-385.
#' }
#'
#'
#' @examples
#' set.seed(1) # determining a seed
#'
#' alpha <- -1.5
#' gama <- 2
#' L <- 8
#'
#' Xobs<- rg0i(100,alpha, gama, L)
#' dg0i(Xobs,alpha, gama, L)
#' qg0i(Xobs,alpha, gama, L)
#'
#' qobs<- runif(10)
#' pg0i(qobs,alpha, gama, L)
#'
#' @useDynLib SITSAR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rgamma
#' @importFrom invgamma rinvgamma

NULL

#' @rdname g0i
#' @export
dg0i <- function(x, alpha, gamma, L) {
    return(.Call('_SITSAR_dg0i', PACKAGE = 'SITSAR', x, alpha, gamma, L))
}

#' @rdname g0i
#' @export
rg0i = function(n,alpha,gamma,L){
  res = stats::rgamma(n,shape=L,rate=L) * invgamma::rinvgamma(n, shape=-alpha, rate = gamma)
  return(res)
}

# ----------------


#' @rdname g0i
#' @export
pg0i <- function(q, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE) {

  if (lower.tail == TRUE) {
    res <- .Call('_SITSAR_pg0i', PACKAGE = 'SITSAR', q, alpha, gamma, L, step = 0.001)
  }else{
    res <- 1 - .Call('_SITSAR_pg0i', PACKAGE = 'SITSAR', q, alpha, gamma, L, step = 0.001)
  }

  if (log.p == TRUE) {
    res_log <- log(res)
  }else{
    res_log <- res
  }
  return (res_log)
}

#' @rdname g0i
#' @export
qg0i <- function(p, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE) {

  if (lower.tail == TRUE) {
    res <- .Call('_SITSAR_qg0i', PACKAGE = 'SITSAR', p, alpha, gamma, L, tol = 1e-6, step = 0.001)
  }else{
    res <- .Call('_SITSAR_qg0i', PACKAGE = 'SITSAR', 1-p, alpha, gamma, L, tol = 1e-6, step = 0.001)
  }

  if (log.p == TRUE) {
    res_log <- log(res)
  }else{
    res_log <- res
  }
  return (res_log)
}

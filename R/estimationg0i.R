#' G0I Estimation
#'
#' \code{g0iestimation} Estimation of G0I distribution parameters
#'
#'
#' @param Theta vector of random initial values considering the assumptions of the parameters, \eqn{alpha < 0}, \eqn{gamma > 0 } and \eqn{L > 0}.
#' @param x vector of observations.
#' @param L_fixed logical; if TRUE, the parameter \eqn{L} is fixed. \eqn{L} is not fixed by default, \eqn{L = FALSE}.
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
#' y <- rg0i(1000, alpha, gama, L)
#'
#' # Estimation with all parameters
#' g0iestimation(c(-1, 1, 1), y)
#'
#' # Estimation with fixed L
#' g0iestimation(c(-1, 1, 1), y, L_fixed = TRUE)
#'
#' @importFrom stats integrate
#' @importFrom fitdistrplus fitdist
#'
#' @export
g0iestimation <- function(Theta,x, L_fixed = FALSE) {

  if ( L_fixed == TRUE){
    fit.params <- suppressWarnings(
      # fitdistrplus::fitdist(x, "g0i",
      #                       start = list(alpha = Theta[1], gamma = Theta[2]),
      #                       lower = c(-Inf, 1e-5),upper = c(-1e-1, Inf),
      #                       fix.arg = list(L =  Theta[3]),
      #                       method = "mle"
      # )
      fitdistrplus::mledist(x, "g0i",
                            start = list(alpha = Theta[1], gamma = Theta[2]),
                            optim.method = "default",
                            fix.arg = list(L =  Theta[3]),
                            gradient = grad_g0i_Lfixed,
                            lower = c(-Inf, 1e-5, .1),
                            upper = c(-1e-1, Inf, Inf)
      )
    )

    res <- c(fit.params$estimate,Theta[3])

    return(res)

  }else{
    fit.params <- suppressWarnings(
      # fitdistrplus::fitdist(x, "g0i",
      #                       start = list(alpha = Theta[1], gamma = Theta[2], L = Theta[3]),
      #                       lower = c(-Inf, 1e-5, .1),upper = c(-1e-1, Inf, 30),
      #                       method = "mle"
      # )
      fitdistrplus::mledist(x, "g0i",
                            start = list(alpha = Theta[1], gamma = Theta[2], L = Theta[3]),
                            gradient = grad_g0i,
                            lower = c(-Inf, 1e-5, .1),
                            upper = c(-1e-1, Inf, 35),
                            optim.method = "default",
      )
    )
    return(fit.params$estimate)
  }
}

#' A function to return the gradient of the log-likelihood for the "BFGS"
#' @noRd
grad_g0i <-  function(par, fix.arg, obs, ddistnam){

  alpha = par[1]
  gama = par[2]
  L = par[3]

  n <- length(obs);

  res1 <- -n * log(gama)  - n * digamma(L - alpha) +
          n * digamma(-alpha) + sum( log(gama + L*obs) )

  res2 <- -n * alpha / gama + (alpha - L) * sum( 1 / (gama + L * obs) )

  res3 <- n * (1 + log(L)) + n * digamma(L - alpha) - n * digamma(L) +
          sum(log(obs)) - sum(log(gama + L * obs)) +
          (alpha - L) * sum(obs/(gama + L*obs))

  res <- c(res1, res2, res3)

  return(res)

}


#' A function to return the gradient of the log-likelihood for the "BFGS"
#' L fixed
#' @noRd
grad_g0i_Lfixed <-  function(par, fix.arg, obs, ddistnam){

  alpha = par[1]
  gama = par[2]

  # Defina L
  L = fix.arg$L

  n <- length(obs);

  res1 <- -n * log(gama)  - n * digamma(L-alpha) + n * digamma(-alpha) +
    sum( log(gama + L*obs) )
  res2 <- -n * alpha / gama + (alpha-L) * sum( 1 / (gama + L * obs) )

  res <- c(res1, res2)

  return(res)

}


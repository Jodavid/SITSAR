#' Parameter estimators for G0I
#'
#' \code{estimatorg0i} Parameter estimators for G0I
#'
#'
#' @param initial vector of initial values.
#' @param Xobs vector of sample for parameters estimators
#' @param method obtained of \code{maxLik}, maximization method, currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. The default method is "NR" for unconstrained problems, and "NM" or "BFGS" for constrained problems, depending on if the grad argument was provided. "BHHH" is a good alternative given the likelihood is returned observation-wise (see maxBHHH).
#'
#' @return \code{dg0i} gives the density, and \code{rg0i} generates random deviates.
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
#' Xobs2<- rg0i(100,alpha,gama,L)
#' estimatorg0i(c(alpha,gama,L), Xobs2)
#' estimatorg0i(c(-1.1,1,1),Xobs2)
#'
#' @importFrom maxLik maxLik
#'
#' @export
estimatorg0i = function(initial, Xobs, method ="NM"){

  A <- rbind(
    c(-1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1))
  B <- c(0, 0, 0)

  rr = maxLik::maxLik(loglikeg0i,
              start=initial,
              Xobs = Xobs,
              constraints=list(ineqA=A, ineqB=B),
              method = method);

  estimate <- stats::coef(rr)
  return(estimate)

}



loglikeg0i = function(theta, Xobs){
  # ..> GI0 log-likelihood
  res <- mean(log(dg0i(Xobs,theta[1],theta[2],theta[3])),na.rm=T)
  return(res)
}

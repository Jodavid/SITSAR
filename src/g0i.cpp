#include <RcppArmadillo.h>
using namespace Rcpp;

// Função para calcular a densidade da distribuição G0I sem Armadillo
// [[Rcpp::export]]
NumericVector dg0i(NumericVector x, double alpha, double gamma, double L) {
  int n = x.size();
  NumericVector densidade(n);
  
  if (alpha >= 0 || gamma <= 0 || L <= 0) {
    stop("Os parâmetros devem satisfazer: alpha < 0, gamma > 0, L > 0.");
  }
  
  double coef = pow(L, L) * tgamma(L - alpha) / (pow(gamma, alpha) * tgamma(-alpha) * tgamma(L));
  
  for (int i = 0; i < n; ++i) {
    if (x[i] > 0) {
      densidade[i] = coef * pow(x[i], L - 1) / pow((gamma + L * x[i]), L - alpha);
    } else {
      densidade[i] = 0;
    }
  }
  
  return densidade;
}

// Função para calcular a CDF da distribuição G0I usando Armadillo para integração
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec pg0i(const arma::vec& q, double alpha, double gamma, double L, double step = 0.001) {
  int n = q.size();
  arma::vec cdf(n, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    arma::vec z_vals = arma::regspace(0.0, step, q[i]);
    NumericVector z_vals_cpp = as<NumericVector>(wrap(z_vals));
    NumericVector dens_vals_cpp = dg0i(z_vals_cpp, alpha, gamma, L);
    arma::vec dens_vals = as<arma::vec>(dens_vals_cpp);
    cdf[i] = arma::sum(dens_vals) * step;
  }
  
  return cdf;
}

// Função para calcular o quantil da distribuição G0I
// [[Rcpp::export]]
NumericVector qg0i(NumericVector p, double alpha, double gamma, double L, double tol = 1e-6, double step = 0.001) {
  int n = p.size();
  NumericVector quantiles(n);
  
  for (int i = 0; i < n; ++i) {
    double lower = 0.0;
    double upper = 10.0; // valor inicial arbitrário para o limite superior
    double mid;
    
    while (upper - lower > tol) {
      mid = (lower + upper) / 2.0;
      arma::vec mid_vec = arma::vec({mid});
      arma::vec cdf_val = pg0i(mid_vec, alpha, gamma, L, step);
      if (cdf_val[0] < p[i]) {
        lower = mid;
      } else {
        upper = mid;
      }
    }
    
    quantiles[i] = mid;
  }
  
  return quantiles;
}

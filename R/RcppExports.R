# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dg0i <- function(x, alpha, gamma, L) {
    .Call(`_SITSAR_dg0i`, x, alpha, gamma, L)
}

pg0i <- function(q, alpha, gamma, L, step = 0.001) {
    .Call(`_SITSAR_pg0i`, q, alpha, gamma, L, step)
}

qg0i <- function(p, alpha, gamma, L, tol = 1e-6, step = 0.001) {
    .Call(`_SITSAR_qg0i`, p, alpha, gamma, L, tol, step)
}


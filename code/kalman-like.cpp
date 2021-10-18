#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//
//
struct KFOUT{
  arma::mat att; // a_t|t
  arma::mat Ptt; // same deal
  double lik;
};


KFOUT kf1step(arma::mat a0, arma::mat P0, arma::mat dt,
              arma::mat ct, arma::mat Tt,
              arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
  
  // Make predictions using a_t
  arma::mat pred = ct + Zt * a0;
  arma::mat vt = yt - pred;
  arma::mat Ft = GGt + Zt * P0 * Zt.t();
  Ft /= 2;
  Ft += Ft.t(); // ensure symmetric
  arma::mat Ftinv = arma::inv(Ft);
  arma::mat Kt = Tt * P0 * Zt.t() * Ftinv;
  
  
  // Incorporate and update information into a_t+1
  arma::mat Lt = Tt - Kt * Zt;
  arma::mat a1 = Tt * a0 + Kt * vt + dt;
  arma::mat P1 = Tt * P0 * Lt.t() + HHt;
  
  
  // Calculate likelihood
  double ftdet = arma::log_det_sympd(Ft);
  double mahalanobis = arma::as_scalar(vt.t() * Ftinv * vt);
  mahalanobis += ftdet;
  double lik = - 0.5 * Ft.n_rows * log(2*M_PI) - 0.5 * mahalanobis;
  KFOUT output = {a1, P1, lik};
  return output;
}

// [[Rcpp::export]]
double kalmanll(arma::mat a0, arma::mat P0, arma::mat dt,
                  arma::mat ct, arma::mat Tt,
                  arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
  
  arma::uword n = yt.n_cols;
  
  double llik = 0;
  double liktmp;
  
  for (arma::uword iter=0; iter<n; iter++) {
    KFOUT step = kf1step(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt.col(iter));
    a0 = step.att;
    P0 = step.Ptt;
    liktmp = step.lik;
    llik += liktmp / n;
  }
  return llik;
}

/*** R
ssmodel_ll <- function(y, simsout, kill = 1e-10) {
  eu = simsout$eu
  if (eu[1] == -2 || eu[2] == -2) return(-Inf)
  G = simsout$G
  C = simsout$C
  M = simsout$M
  nx = simsout$nx
  ny = simsout$ny
  ns = simsout$ns
  
  G[abs(G) < kill] = 0
  C[abs(C) < kill] = 0
  M[abs(M) < kill] = 0
  
  
  ## write the system in terms of
  ## yo_t+1 = c+ Z*x_t + G*v_t+1
  ## x_t+1 = d+ T*x_t   + H*w_t+1
  Zt = cbind(diag(ny), matrix(0, ny, nx - ny))   # it is not clear to me that the observations are in the first 7 (last 7?)
  dt = as.matrix(C)
  Tt = G
  HH = M %*% t(M)
  GG = diag(0, ny) # note: no observation errors?
  ct = matrix(0, ny, 1)
  
  
  a0 = matrix(0, nx, 1)
  P0 = 10 * diag(nx) # as in the matlab code
  ## P0 = matrix(solve(diag(nx^2) - T %x% T) %*% c(HH), nx, nx)
  ## prior is the steady state (way too slow, c(HH) is wrong)
  out = kalmanll(a0, P0, dt, ct, Tt, Zt, HH, GG, y)
  out
}

*/

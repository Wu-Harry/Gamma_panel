#include <Rcpp.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;

// C++ function for estimating beta and theta

// [[Rcpp::export]]
NumericMatrix M_step_1(NumericMatrix E_W, NumericVector E_xi, NumericMatrix X, NumericVector grid, 
                     NumericVector C, NumericVector beta_0, double theta_0,
                     int m, int n) {
  // Update theta
  double theta = theta_0;
  
  // Update beta
  int len = beta_0.size();
  NumericVector grad(len);
  NumericMatrix H(len, len);
  NumericMatrix l_sum(m, 2);
  NumericMatrix l_sum_x(m, len);
  NumericVector Exp_Xbeta(0);
  for (int i = 0; i < n; i++){
    double temp = 0;
    for (int j = 0; j < len; j++){
      temp += X(i,j) * beta_0[j];
    }
    Exp_Xbeta.push_back(exp(temp));
  }
  
  for (int l = 0; l < m; l++) {
    double temp = 0;
    for (int i = 0; i < n; i++) {
      if (C[i] >= grid[l]) temp += E_xi[i] * Exp_Xbeta[i];
    }
    l_sum(l,0) += temp;
    
    NumericVector tempV(len);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < len; j++){
        if (C[i] >= grid[l]) tempV[j] += E_xi[i] * Exp_Xbeta[i] * X(i,j);
      }
    }
    for (int j = 0; j < len; j++) {
      l_sum_x(l,j) += tempV[j];
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int l = 0; l < m; l++) {
      DoubleVector gr1(len);
      if (C[i]>=grid[l]) {
        for (int j = 0; j < len; j++) {
          gr1[j] = E_W(i,l) * (X(i,j) - l_sum_x(l,j) / l_sum(l,0));
        }
      }
      for (int j = 0; j < len; j++) {
        grad[j] += gr1[j];
        for (int k = 0; k < len; k++) {
          H(j,k) -= gr1[j] * gr1[k];
        }
      }
    }
  }
  
  NumericMatrix Out(len, len + 2);
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      Out(i,j) = H(i,j);
    }
    Out(i, len) = grad[i];
  }
  Out(0, len + 1) = theta;
  
  return Out;
}




/*** R

*/

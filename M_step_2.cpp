#include <Rcpp.h>
using namespace Rcpp;

// C++ function for estimating Lambda

// [[Rcpp::export]]
NumericVector M_step_2(NumericMatrix E_W, NumericVector E_xi, NumericMatrix X, NumericVector grid, 
                       NumericVector C, NumericVector beta_0, int m, int n) {
  
  // Update lambda
  int len = beta_0.size();
  int len_g = grid.size();
  NumericVector lam_new(len_g);
  NumericMatrix l_sum(m, 2);
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
    
    temp = 0;
    for (int i = 0; i < n; i++) {
      if (C[i] >= grid[l]) temp += E_W(i,l);
    }
    l_sum(l,1) += temp;
    
    lam_new[l] = l_sum(l,1) / l_sum(l,0);
  }
  
  return lam_new;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


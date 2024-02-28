#include <Rcpp.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix E_step(NumericMatrix Time, LogicalMatrix delt, IntegerVector K,
                     NumericVector lam_0, NumericVector para_0,
                     NumericVector grid, NumericVector Exp,
                     NumericVector nodes, NumericVector whts, double n, bool SE) {
  NumericMatrix E_W(n, grid.size() + 3); // matrix S to store function values
  NumericVector E_xi(n), E_log_xi(n);
  double Loglike = 0;
  
  for (int i = 0; i < n; i++) {
    NumericMatrix S(3 + grid.size(), nodes.size());
    double Pi1 = 0;
    NumericVector delta_ij_vec(0);
    //printf("1 ");
    
    for (int j = 0; j < K[i]; j++) {
      double delta_ij = 0;
      NumericVector l_index(0);
      for (int l = 0; l < grid.size(); l++) {
        if (grid[l] > Time(i,j) & grid[l] <= Time(i,j + 1)) {
          l_index.push_back(l);
          delta_ij += Exp[i] * lam_0[l];
        }
      }
      //printf("%1f ", delta_ij);
      delta_ij_vec.push_back(delta_ij);
      Pi1 += (1 - delt(i,j)) * delta_ij;
    }
    Pi1 += para_0[0];
    //printf("Pi=%1f ", Pi1);
    //printf("del=%1f,%1f,%1f,%1f,%1f", delta_ij_vec[1],delta_ij_vec[2], delta_ij_vec[3],delta_ij_vec[4],delta_ij_vec[5]);
    //printf("2 ");
    
    for (int k = 0; k < nodes.size(); k++) {
      double s = 0;
      for (int j = 0; j < K[i]; j++) {
        s += log(1 - delt(i,j) * exp(-delta_ij_vec[j]*nodes[k]/Pi1));
        }
      //printf("2.1 ");
      //printf("s=%1f ", s);
      S(1,k) = exp((para_0[0]-1)*log(nodes[k]/Pi1)+ s) / Pi1;
      S(2,k) = exp((para_0[0])*log(nodes[k]/Pi1)+ s) / Pi1;
      S(3,k) = log(nodes[k]/Pi1)*S(1,k);
      //printf("2.2 ");
      
      for (int j = 0; j < K[i]; j++) {
        NumericVector l_index(0);
        for (int l = 0; l < grid.size(); l++) {
          if (grid[l] > Time(i,j) & grid[l] <= Time(i,j + 1)) {
            l_index.push_back(l);
            }
          }
        for (int m = 0; m < l_index.size(); m++) {
          S(l_index[m] + 3,k) = delt(i,j) * nodes[k] * pow(Pi1,-1) * lam_0[l_index[m]] * Exp[i] * S(1,k) / (1 - exp(-nodes[k] * delta_ij_vec[j] * pow(Pi1,-1)));
          }
        }
      //printf("2.3 ");
      }
    //printf("3 ");
    
    double denominator = 0;
    for (int k = 0; k < nodes.size(); k++) {
      denominator += S(1,k) * whts[k];
    }
    
    for (int k = 0; k < nodes.size(); k++) {
      E_xi[i] += S(2,k) / denominator * whts[k];
      E_log_xi[i] += S(3,k) / denominator * whts[k];
    }
    for (int l = 0; l < grid.size(); l++) {
      if (grid[l] > Time(i,0) && grid[l] <= Time(i,K[i])) {
        for (int k = 0; k < nodes.size(); k++) {
          E_W(i,l) += S(l+3,k) / denominator * whts[k];
          }
        }
      }
    //printf("4 ");
    //printf("%d\n",i);
    if (SE == 1) Loglike = Loglike + ( log(para_0[0])*para_0[0] - log(tgamma(para_0[0])) + log(denominator) );
    }
  
  for (int i = 0; i < n; i++) {
    E_W(i, grid.size()) = E_xi[i];
    E_W(i, grid.size() + 1) = E_log_xi[i];
    }
  E_W(0, grid.size() + 2) = Loglike;
  //printf("5");
  
  return E_W;
  }



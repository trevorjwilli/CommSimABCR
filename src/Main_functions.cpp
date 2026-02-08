#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int samp_com_for_birth_cpp(NumericMatrix x) {
  int nrow = x.nrow(); int ncol = x.ncol();
  NumericVector sums(nrow);
  NumericVector probabilities(nrow);
  
  double csize = 0.0;
  
  for (int i = 0; i < nrow; i++) {
    double total = 0.0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
      csize += x(i, j);
    }
    sums[i] = total;
  }
  
  for (int i = 0; i < nrow; i++) {
    probabilities[i] = sums[i] / csize;
  }
  
  IntegerVector vec(nrow);
  for (int i = 0; i < nrow; i++) {
    vec[i] = i;
  }
  
  IntegerVector out_vec = Rcpp::sample(vec, 1, false, probabilities);
  return out_vec[0];
}
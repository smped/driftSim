#include <Rcpp.h>
using namespace Rcpp;

//' @title  Simulate the progeny
//'
//' @description Simulates the progeny for a given breeding pair, to a specified litter size
//'
//' @param pair A vector of alleles, of length 2.
//' Alleles must be represented as 0, 1 or 2.
//' @param litter The litter size from the breeding pair.
//'
//' @return The alleles of progeny, also using the 0, 1 or 2 representation
//'
//' @examples
//' breed(c(1, 1))
//'
//' # Sampling a population
//' pop <-matrix(sample(c(0, 1, 2), 100, TRUE), ncol = 2)
//' progeny <- apply(pop, 1, breed)
// [[Rcpp::export]]
NumericVector breed(IntegerVector pair, int litter){

  // Check the pair only has two members
  int n = pair.size();
  // Check all values are with 0, 1 or 2
  if (n != 2) stop("Pairs can only be specified with 2 values");
  if (max(pair) > 2 || min(pair) < 0) stop("Alleles can only be between 0 and 2");
  NumericVector m = rbinom(litter, 1, pair[0]/2.0);
  NumericVector f = rbinom(litter, 1, pair[1]/2.0);
  return m + f;
}

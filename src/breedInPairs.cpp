#include <Rcpp.h>
using namespace Rcpp;

//' @title  Simulate the progeny
//'
//' @description Simulates the progeny for a set of breeding pairs, to a specified litter size
//'
//' @param pairs A two-column matrix of alleles.
//' Alleles must be represented as 0, 1 or 2.
//' @param litter The litter size from the breeding pair.
//'
//' @return The alleles of progeny, also using the 0, 1 or 2 representation. Returned as a vector
//'
//' @examples
//' breed(c(1, 1), litter = 6)
//'
//' # Sampling a population
//' pop <-matrix(sample(c(0, 1, 2), 10, TRUE), ncol = 2)
//' apply(pop, 1, breed, litter = 6)
//' breedInPairs(pop, litter = 6)
//' @export
// [[Rcpp::export]]
NumericVector breedInPairs(IntegerMatrix pairs, int litter){

  int npairs = pairs.nrow(), ncol = pairs.ncol();

  // Check the pairs only has two members
  if (ncol != 2) stop("Pairs can only be specified with 2 values");

  // Check all values are with 0, 1 or 2
  if (max(pairs) > 2 || min(pairs) < 0) stop("Alleles can only be between 0 and 2");

  NumericVector out(npairs*litter);
  int counter = 0;

  for (int i = 0; i < npairs; i++){ // For each pair (row) in the matrix
    for (int j = 0; j < litter; j++){ // Generate a child for each member of the litter
      NumericVector m = rbinom(1, 1, pairs(i, 0)/2.0);
      NumericVector f = rbinom(1, 1, pairs(i, 1)/2.0);
      out(counter) = m(0) + f(0);
      counter++;
      //std::cout<<i<<" "<<j<<" "<<counter<<" "<<m(0)<<" "<<f(0)<<" "<<std::endl;
    }
  }

  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

//' @title Set migration Probabilities
//'
//' @description Defines the migration rate using a radial model with Pop1 as the central population.
//' All others form a ring aorund this population.
//'
//' @param n The number of populations
//' @param mig The migration rate
//'
//' @return A matrix with migration probabilities
//'
//' @examples
//' setMigProbs(6, 0.1)
//' @export
// [[Rcpp::export]]
NumericMatrix setMigProbs(int n, double mig){
  // Create the matrix
  NumericMatrix out(n, n);

  // Setup Row & Colnames
  CharacterVector nm(n);

  // Fill in migration probs for self pops
  for (int i = 0; i <n; i++){
    // Generate row & column names
    std::stringstream s;
    s<<(i+1);
    nm(i) = "Pop" + s.str();
    out(i, i) = 1 - mig;
  }

  for (int i = 0; i < (n - 1); i++){

    // Add the central population
    out(i + 1, 0) = mig / (n - 1);
    out(0, i + 1) = mig / 3;
    // The clockwise neighbours
    int fwd = (i + 1) % (n - 1) + 1;
    out(fwd, i + 1) = mig / 3;
    // The anti-clockwise neighbours
    int back = (i + n - 2) % (n - 1) + 1;
    out(back, i + 1) = mig / 3;

  }

  rownames(out) = nm;
  colnames(out) = nm;

  return out;
}

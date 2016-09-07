#' @title  Simulate the progeny
#'
#' @description Simulates the progeny for a given breeding pair, to a specified litter size
#'
#' @param pair A vector of alleles, of length 2.
#' Alleles must be represented as 0, 1 or 2.
#' @param litter The litter size from the breeding pair.
#'
#' @return The alleles of progeny, also using the 0, 1 or 2 representation
#'
#' @examples
#' breed(c(1, 1))
#'
#' # Sampling a population
#' pop <-matrix(sample(c(0, 1, 2), 100, TRUE), ncol = 2)
#' progeny <- apply(pop, 1, breed)
#'
#' @export
breed <- function(pair, litter = 6){

  pair <- as.integer(pair)
  stopifnot(any(pair %in% c(0, 1, 2)))
  stopifnot(length(pair) ==2)
  litter <- as.integer(litter)
  m <- rbinom(litter, 1, pair[1]/2) # Gametes from individual 1, nominally m
  f <- rbinom(litter, 1, pair[2]/2) # Gametes from individual 2, nominally f
  m + f # Progeny

}

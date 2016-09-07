#' @title Simulate Drift
#'
#' @description Simulate drift across a number of populations and generations, for a single allele
#'
#' @param f0 The allele frqeuency pre-bottleneck
#' @param N0 The effective population size pre-bottleneck
#' @param Nt The effective population size after 't' generations
#' @param t The number of generations
#' @param n The number of populations
#' @param mig The migration rate
#' @param surv The survival rate at the original population bottleneck
#' @param litter The litter size at each generation
#'
#' @return A list with components \code{ft} and \code{nEff}.
#' These denote the final allele frequency, and the effective populations sizes at each geenration respectively.
#'
#' @examples
#' simDrift(f0 = 0.8, N0 = 100, Nt = 200, t = 10, n = 6, mig = 0.01, surv = 0.1, litter= 6)
#'
#' @import magrittr
#'
#' @export
simDrift <- function(f0, N0, Nt, t, n, mig, surv, litter, ...){

  # Convert all to integers where required
  N0 <- as.integer(N0)
  Nt <- as.integer(Nt)
  t <- as.integer(t)
  n <- as.integer(n)
  litter <- as.integer(litter)

  # Check all variables
  stopifnot(c(f0 > 0, f0 < 1)) # The original frequency must be between 0 & 1
  stopifnot(!is.na(c(N0, Nt, t, n, litter))) # Ensure all required as integers were specified as numeric
  stopifnot(c(t, n) > 0) # There must be at least one generation & one population
  stopifnot(c(mig >= 0, mig < 1)) # Zero migration is acceptable
  stopifnot(c(surv > 0, surv <= 1)) # 100% survival is aceptable if there is no bottleneck
  stopifnot(litter >= 3) # Litters must be greater than 2 for population growth

  # Simulate the starting populations
  pop0 <- replicate(n, rbinom(2*N0, 2, 1 - f0), simplify = FALSE)
  # Introduce the bottleneck
  keep <- replicate(n, as.logical(rbinom(2*N0, 1, surv)), simplify = FALSE)
  pop0 <- mapply(function(x, y){x[y]}, x = pop0, y = keep, SIMPLIFY = FALSE)
  pop0
  # Check that all populations are viable (i.e have >= 2 members)
  viable <- vapply(pop0, function(x){length(x)>=2}, logical(1))
  if (any(!viable)) stop("One or more populations are too small to be viable post-bottleneck")

  # Estimate population growth for each population
  genSizes <- lapply(pop0, function(x){
    seq(log(length(x)/2), log(Nt), length.out = t + 1) %>%
      exp %>%
      round(0) %>%
      extract(-1)
  })

  # Form into breeding pairs.
  # If an odd number of individuals is in a population, assume polygamy and use recursion
  pairs <- vector("list", t + 1)
  progeny <- vector("list", t)
  pairs[[1]] <- lapply(pop0, function(x){suppressWarnings(matrix(x, ncol = 2))})

  for (i in 1:t){
    # Breed at every iteration
    progeny[[i]] <- lapply(pairs[[i]],
                           function(x){
                             apply(x, 1, breed, litter = litter) %>% as.vector()
                           })

    # Allow migration pre-survival
    progeny[[i]] <- migrate(progeny[[i]], mig, ...)

    # Check population sizes are appropriate
    nKeep <- lapply(genSizes, magrittr::extract, i)
    popVsKeep <- mapply(function(pop, n){length(pop) > 2*n},
                        pop = progeny[[i]],
                        n = nKeep,
                        SIMPLIFY = FALSE)
    if (any(!unlist(popVsKeep))) stop("Breeding rates unable to give required final population size")

    # Select the population for breeding pairs in the next generation
    pairs[[i + 1]] <- mapply(function(pop, n){
      pairMat <- suppressWarnings(matrix(pop, ncol = 2))
      ind <- sample.int(nrow(pairMat), n)
      pairMat[ind,]
    },
    pop = progeny[[i]],
    n = nKeep,
    SIMPLIFY = FALSE)

  }

  # Return population 1 as the population of interest
  list(ft = mean(2 - progeny[[t]][[1]])/2,
       nEff = genSizes[[1]])

}

#' @title Simulate unbounded migration
#'
#' @description Simulate migration between populations.
#' Equal migration between populations is assumed
#' The size of the first population is used to determine the number of migrants
#'
#' @param pops A list of populations with each containing a vector of individuals as 0, 2 or 2
#' @param rate The upper bound on the migration rate
#' @param method Does migration happen randomly between populations (method = "random")
#' or in a circular direction (method = "circular").
#'
#' @return A list of populations with the identical structure as supplied in pops.
#' The rate of shuffling between populations will have been performed according
#' to the migration rate and method specified in the function call.
#'
#' @examples
#' # Simulate 3 populations
#' pops <- lapply(c(10, 12, 14), function(x){sample(0:2, x, replace = TRUE)})
#' names(pops) <- paste0("Pop", LETTERS[seq_along(pops)])
#' # Add names so the individuals can be tracked
#' pops <- sapply(names(pops),
#'                function(x){
#'                  pop <- pops[[x]]
#'                  names(pop) <- paste(x, letters[seq_along(pop)], sep=":")
#'                  pop},
#'                simplify = FALSE)
#' # Now perform simulated migration
#' migrate(pops, 0.1, "random")
#' migrate(pops, 0.1, "circular")
#'
#' @export
migrate <- function(pops, rate, method = "random"){

  stopifnot(is.list(pops))
  stopifnot(length(pops) > 1)
  stopifnot(vapply(pops, is.vector, logical(1)))
  stopifnot(method %in% c("random", "circular"))

  # Check each population only has 0, 1, 2 values
  all012 <- vapply(pops, function(x){all(x %in% 0:2)}, logical(1))
  stopifnot(all012)

  nMig <- round(rate*length(pops[[1]]), 0)
  if (nMig > 0 ){ # Only migrate if the first population is big enough

    mvSamps <- lapply(pops, function(x){
      sample.int(length(x), nMig)
    })
    # Remove the individuals from each population
    rmPops <- mapply(function(x, rm){x[-rm]},
                     x = pops,
                     rm = mvSamps,
                     SIMPLIFY = FALSE)
    # Collect the individuals to move
    mvPops <- mapply(function(x, rm){x[rm]},
                     x = pops,
                     rm = mvSamps,
                     SIMPLIFY = FALSE)
    if (method == "random") {
      # Shuffle the populations
      mvPops <- mvPops[sample.int(length(mvPops))]
    }
    if (method == "circular") {
      # Move in a step-wise manner
      mvPops <- mvPops[c(seq(2, length(mvPops)), 1)]
    }
    # Reform the populations
    pops <- mapply(function(a, b){c(a, b)},
                   a = rmPops,
                   b = mvPops,
                   SIMPLIFY = FALSE)
  }

  pops

}

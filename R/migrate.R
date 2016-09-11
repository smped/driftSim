#' @title Simulate migration
#'
#' @description Simulate migration between populations.
#' Equal migration between populations is assumed
#' The size of the first population is used to determine the number of migrants
#'
#' @param pops A list of populations with each containing a vector of individuals as 0, 2 or 2
#' @param rate The upper bound on the migration rate
#'
#' @return A list of populations with the same length as supplied in pops.
#' The rate of shuffling between populations will have been performed according
#' to the migration rate specified in the function call.
#' Population sizes may vary slightly post migration
#'
#' @details To model migration, the first population is taken as the centre of a circle.
#' All populations are in contact with this population, but only have two neighbours.
#' Migration can thus occur between the central population and all others, but only between neighbours in the outer ring.
#' This is done by assigning probabilites of population membership to each individual based on it's original population.
#'
#' @examples
#' # Simulate 5 populations
#' pops <- lapply(c(10, 12, 14, 10, 10), function(x){sample(0:2, x, replace = TRUE)})
#' names(pops) <- paste0("Pop", LETTERS[seq_along(pops)])
#' # Now perform simulated migration
#' migrate(pops, 0.1)
#'
#' @import dplyr
#' @import magrittr
#'
#' @export
migrate <- function(pops, rate){

  stopifnot(is.list(pops))
  stopifnot(length(pops) > 1)
  stopifnot(vapply(pops, is.vector, logical(1)))
  stopifnot(rate >= 0, rate <= 1)

  # Check each population only has 0, 1, 2 values
  all012 <- vapply(pops, function(x){all(x %in% 0:2)}, logical(1))
  stopifnot(all012)

  # Collect the initial population info for each individual into a data.frame
  nPops <- length(pops)
  popDf <- lapply(seq_along(pops), function(x){
    dplyr::data_frame(AlleleCount = pops[[x]],
               SourcePop = paste0("Pop", x))
  }) %>%
    dplyr::bind_rows()

  # Assign migration probablities for individuals
  nIndividuals <- nrow(popDf)
  probs <- matrix(0, nrow = nIndividuals, ncol = nPops)
  # Probabilities of membership for Pop1
  probs[, 1] <- rate / 3
  probs[popDf$SourcePop == "Pop1", 1] <- 1 - rate
  probs[popDf$SourcePop == "Pop1", 2:nPops] <- rate/(nPops - 1)
  neighbourMat <- cbind(curPop = paste0("Pop", 2:nPops),
                        rightPop  = paste0("Pop", 2:nPops)[c(2:(nPops -1),1)],
                        leftPop = paste0("Pop", 2:nPops)[c((nPops - 1),1:(nPops -2))])
  # Assign migration probabilities for other populations
  for (i in 2:nPops){
    curPop <- paste0("Pop", i)
    neighbours <- neighbourMat[i - 1, c("rightPop", "leftPop")]
    probs[popDf$SourcePop == curPop, i] <- 1 - rate
    probs[popDf$SourcePop %in% neighbours, i] <- rate / 3
  }

  # Assign post-migration population membership for each individual
  postPops <- apply(probs, MARGIN = 1, FUN = function(x){
    sample.int(n = nPops, size = 1, prob = x)
  })

  # Produce the output
  popDf %>%
    dplyr::mutate(PostPop = postPops) %>%
    split(f = .$PostPop) %>%
    lapply(magrittr::extract2, "AlleleCount") %>%
    magrittr::set_names(names(pops))

}

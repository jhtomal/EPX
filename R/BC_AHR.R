#' Calculate AHR (or expected AHR)
#'
#' Calculates average hitrate (AHR), which is equivalent to average precision.
#' When there are ties in ranking (not all values in phat are unique), the
#' result is the expection of AHR. The algorithm that produces this analytic
#' result assumes that the items in any tied group are in an arbitrary order
#' within the group.
#'
#' Implementation adapted from Wang (2005, Chapter 3). Please also see Chapter 3
#' of Wang (2005) for AHR and expected AHR formulas.
#'
#' @references
#' Wang, M. (2015).
#' \emph{Statistical Methods for High Throughput Screening Drug Discovery Data}.
#' University of Waterloo.
#' Retrieved from \url{http://etd.uwaterloo.ca/etd/y32wang2005.pdf}
#'
#' @param y True (binary) response vector where 1 is the rare/relevant class.
#' @param phat Numeric vector of estimated probabilities of relevance.
#' @param ... Further arguments passed to or from other methods.
#' @return Numeric value of average hitrate; expected average hitrate when there
#'   are ties.
#' @examples
#' ## AHR when there are no ties in phat:
#' resp <- c(1, 0, 0, 0, 1)
#' prob <- (1:5)*0.1
#' AHR(y = resp, phat = prob)
#' # expect answer: 1/2 * (1 + 0 + 0 + 0 + 2/5)
#'
#' ## (Expected) AHR when there are ties in phat:
#' resp <- c(1, 1, 0,   0,   0,   0,   0,    1,   0, 0)
#' prob <- c(1, 1, 1, 0.4, 0.4, 0.3, 0.2, 0.15, 0.1, 0)
#' AHR(y = resp, phat = prob)
#' # expect answer: 1/3 * (2/3 + 1/2 * (1/3 + 2/3) + 1/3 * 4/3 +
#' #                       1/8 * (2/3 + 2/3 + 2/3 + 1))
#' @useDynLib EPX preci
#' @export

AHR <- function(y, phat, ...) {

  precision <- 0
  n <- length(y)

  pre.sum <- .C("preci",
                n = as.integer(n),
                score =  as.double(phat),
                trueclass = as.integer(y),
                precision = as.double(precision))

  return(pre.sum$precision)

}

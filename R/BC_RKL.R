#' Calculate rank last
#'
#' The peformance measure rank last (RKL) is calculated as follows: after
#' ranking the observations in decreasing order via \code{phat}, RKL is the rank
#' of the last truly relevant observation. Hence, RKL can take on integer values
#' from 1 to \eqn{n}, where \eqn{n} is the total number of observations. If
#' there are ties, the last object in the tied group determines RKL. That is, if
#' all \eqn{n} objects are tied at the first rank but only one object is truly
#' relevant, RKL will have a value of \eqn{n}.
#'
#' @param y True (binary) response vector where 1 is the rare/relevant class.
#' @param phat Numeric vector of estimated probabilities of relevance.
#' @param ... Further arguments passed to or from other methods.
#' @return Numeric value of RKL.
#' @examples
#' ## without ties in phat
#' resp <- c(rep(1, 50), rep(0, 50))
#' prob <- (1:100)*0.01
#' RKL(y = resp, phat = prob) # expect 100
#'
#' resp <- c(rep(0, 50), rep(1, 50))
#' RKL(y = resp, phat = prob) # expect 50
#'
#' ## with ties in phat
#' resp <- sample(c(1, 0), 100, replace = TRUE)
#' prob <- rep(1, 100)
#' RKL(y = resp, phat = prob) # expect 100
#' @export

RKL <- function(y, phat, ...) {

  y <- data.matrix(y)
  phat <- data.matrix(phat)

  phat.sorted <- rev(sort(phat)) # largest to smallest so indices = rank
  y.sorted <- rev(y[order(phat)])

  rkl.noties <- max(which(y.sorted == 1))
  rkl.noties.ind <- phat.sorted[rkl.noties]
  max.of.ties <- max(which(phat.sorted == rkl.noties.ind))

  return(max.of.ties)

}

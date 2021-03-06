#' Calculate Initial Enhancement
#'
#' Calculates initial enhancement (IE), which is the precision at one specific
#' shortlist length (cutoff) normalised by the proportion of relevants in the
#' total sample size (Tomal et al. 2015). Since IE is a rescaling of
#' precision, we expect IE and AHR to lead to similar conclusions as an
#' assessment metric for the EPX algorithm.
#'
#' Let \eqn{c} be the cutoff and \eqn{h(c)} be the hitrate at \eqn{c}. Let also
#' \eqn{A} be the total number of relevants and \eqn{N} be the total number of
#' observations. IE is defined as \deqn{IE = h(c) / (A / N)}
#' IE calculation does not change whether there are ties in \code{phat} or not.
#'
#' @references
#' Tomal, J. H., Welch, W. J., & Zamar, R. H. (2015).
#' Ensembling classification models based on phalanxes of variables with
#' applications in drug discovery.
#' \emph{The Annals of Applied Statistics},
#' \emph{9}(1), 69-93.
#' \url{http://doi.org/10.1214/14-AOAS778}
#'
#' @param y True (binary) response vector where 1 is the rare/relevant class.
#' @param phat Numeric vector of estimated probabilities of relevance.
#' @param cutoff Shortlist cutoff length, and so must not exceed length of
#'   \code{y}; default is half the sample size.
#' @param ... Further arguments passed to or from other methods.
#' @return Numeric value of IE.
#' @examples
#' ## IE when there are no ties in phat:
#' \donttest{
#' resp <- c(1, 1, 0,   0,   0,   0,   0,    1,   0, 0)
#' prob <- (10:1) * 0.1
#' IE(y = resp, phat = prob, cutoff = 3)
#' # expect answer: (2/3) / (3/10)
#'
#' ## IE when there are ties
#' resp <- c(1, 1, 0,   0,   0,   0,   0,    1,   0, 0)
#' prob <- c(1, 1, 1, 0.4, 0.4, 0.3, 0.2, 0.15, 0.1, 0)
#' IE(y = resp, phat = prob, cutoff = 3)
#' }
#' # expect answer: same as above
#' @export
IE <- function(y, phat, cutoff = length(y) / 2, ...) {

  y <- data.matrix(y)
  phat <- data.matrix(phat)  # otherwise get.nhit.all complains

  ranked.phat <- rev(sort(phat))  # largest to smallest so index = rank
  ranked.y <- rev(y[order(phat)])

  hn <- sum(ranked.y[1:cutoff]) / cutoff # hit rate for top n relevant = Hn/n
  return(hn / (sum(y) / length(y)))

}

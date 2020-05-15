#' Plot hit curve for an "\code{epx}" object
#'
#' Plots the hit curve for the fitted values of an "\code{\link{epx}}" object.
#'
#' Order the cases by decreasing \code{phat} (predicted probabilities of
#' relevance) values, and plot the expected number and actual number of hits as
#' cases are selected. Cases with tied \code{phat} values are grouped together.
#' See \link{hit.curve} in order to plot a hit curve in general.
#'
#' @param x Object of class "\code{\link{epx}}".
#' @param max.cutoff Maximum number of observations selected, equivalently the
#'   maximum shortlist cutoff; default is \code{min(100, length(x$Y))}.
#' @param plot.hc Whether to return a plot of the hit curve; default is
#'   \code{TRUE}.
#' @param ... Further arguments passed to or from other methods.
#' @return Plot of the hit curve (if \code{plot.hc = TRUE}) and a list with the
#'   following vectors:
#' \item{select}{Number of observations in each tied \code{phat} group;
#' \code{select[1]}, \code{select[2]}, \code{...} are the numbers of
#' observations with the largest predicted probability of relevance
#' (\code{max(phat)}), the second largest value in \code{phat}, etc.}
#' \item{p}{Unique \code{phat} values; \code{p[1]}, \code{p[2]}, \code{...} are
#' the largest value in \code{phat}, the second largest value in \code{phat},
#' etc.}
#' \item{nhits}{Number of hits (truly relevant observations) in each tied
#' \code{phat} group.}
#' \item{nhitlast}{Number of hits after \code{max.cutoff} observations
#' selected.}
#' @examples
#' # Example with data(harvest)
#'
#' ## Phalanx-formation using a base classifier with 500 trees (default)
#' \donttest{
#' set.seed(761)
#' model <- epx(x = harvest[, -4], y = harvest[, 4],
#'              classifier.args = list(ntree = 500))
#'
#' ## Hit curve for model with default settings
#' model.hc <- plot(model)
#'
#' ## In the top 100 ranked observations selected, the number that are truly
#' ## relevant is
#' model.hc$nhitlast
#'
#' ## Hit curve with max.cutoff at 150 (Note: Commented off for time.)
#' model.hc.150 <- plot(model, max.cutoff = 150)
#' model.hc.150$nhitlast  # Number of hits in top 150 ranked observations.
#' }
#' @export
plot.epx <- function(x,
                     max.cutoff = min(100, length(x$Y)),
                     plot.hc = TRUE, ...) {

  return(hit.curve(phat = x$ENSEMBLED.FITS,
                   y = x$Y,
                   max.cutoff = max.cutoff,
                   plot.hc = plot.hc))

}


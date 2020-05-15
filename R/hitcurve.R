#' Plot hit curve
#'
#' Plots the hit curve corresponding to \code{phat} and \code{y}.
#'
#' Order the cases by decreasing \code{phat} (predicted probabilities of
#' relevance) values, and plot the expected number and actual number of hits as
#' cases are selected. Cases with tied \code{phat} values are grouped together.
#' See \link{plot.epx} for plotting the hit curve for an "\code{\link{epx}}"
#' object.
#'
#' @param y True binary response vector where 1 denotes the relevant rare class.
#' @param phat Vector of estimated probabilities of relevance.
#' @param max.cutoff Maximum number of observations selected, equivalently the
#'   maximum shortlist cutoff; default is \code{min(100, length(y))}.
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
#' ## Plot hit curve for cross-validated predicted probabilities of relevence
#' set.seed(761)
#' model.cv <- cv.epx(model)
#' preds.cv <- model.cv[-nrow(model.cv), ncol(model.cv)]
#' cv.hc <- hit.curve(phat = as.numeric(preds.cv), y = model$Y)
#' }
#' @export

hit.curve <- function(y,
                      phat,
                      max.cutoff = min(100, length(y)),
                      plot.hc = T, ...) {

  # unique phats, sorted with largest first.
  uniq.phat <- rev(sort(unique(phat)))

  select <- vector("numeric", length(uniq.phat))
  nhits  <- vector("numeric", length(uniq.phat))

  for (i in 1:length(uniq.phat)) {

    cases.sel <- (phat == uniq.phat[i])
    select[i] <- sum(cases.sel)
    nhits[i]  <- sum(y[cases.sel])

    if (sum(select[1:i]) >= max.cutoff || i == length(uniq.phat)) {
      i.max <- i
      break
    }

  }

  uniq.phat <- uniq.phat[1:i.max]
  nhits     <- nhits[1:i.max]
  select    <- select[1:i.max]
  exp.hits  <- uniq.phat * select

  # Plot hitcurve
  if (plot.hc) {

    # Just set up axes, etc.
    x.max <- sum(select)
    y.max <- max(sum(nhits), sum(exp.hits))
    plot(1, 1, type = "n", xlim = c(0, x.max), ylim = c(0, y.max),
         xlab = "Number of relevants selected",
         ylab = "Number of actual hits")

    # Plot the cumulative number of hits.
    cum.select.last <- 0
    cum.nhits.last <- 0
    cum.exp.hits.last <- 0

    for (i in 1:i.max) {

      cum.select <- cum.select.last + select[i]
      cum.nhits <- cum.nhits.last + nhits[i]
      cum.exp.hits <- cum.exp.hits.last + uniq.phat[i] * select[i]

      # Always plot points.
      points(cum.select, cum.nhits, pch = "+", cex = par()$cex * 1) #0.9)
      # points(cum.select, cum.exp.hits, pch = "O",
      #        cex = par()$cex * 1) #0.7)

      # Join points with lines if gap is sufficient.
      if (select[i] / x.max > 0.02) {
        lines(c(cum.select.last + 1, cum.select),
              c(cum.nhits.last + nhits[i] / select[i], cum.nhits),
              lty = 1)
        # lines(c(cum.select.last + 1, cum.select),
        #       c(cum.exp.hits.last + uniq.phat[i], cum.exp.hits),
        #       lty = 2)
      }

      cum.select.last <- cum.select
      cum.nhits.last <- cum.nhits
      cum.exp.hits.last <- cum.exp.hits

    }

  }

  # Calculate number of hits after max.cutoff cases selected
  nhit.select <- 0
  select.cum <- cumsum(select)
  nhits.cum <- cumsum(nhits)

  if (length(select) == 1) {
    nhit.select <- (nhits/select) * max.cutoff
  } else {
    for (i in 1:(length(uniq.phat) - 1)) {

      if(select.cum[i] == max.cutoff) {
        nhit.select <- nhits.cum[i]; break
      }

      if(select.cum[i+1] == max.cutoff) {
        nhit.select <- nhits.cum[i+1]; break
      }

      if(max.cutoff > select.cum[i] & max.cutoff<select.cum[i+1]) {
        nhit.select <- nhits.cum[i] +
          nhits[i+1] / select[i+1] * (max.cutoff-select.cum[i]); break
      }
    }
  }

  return(list(select = select,
              p = uniq.phat,
              nhits = nhits,
              nhitlast = nhit.select))
}

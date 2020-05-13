.getPerformanceMeasure <- function(performance, performance.args, ...) {

  performance <- match.arg(arg = performance,
                           choices = c("AHR", "IE", "TOP1", "RKL"),
                           several.ok = FALSE)

  # AHR ------------------------------------------------------------------------
  if (performance == "AHR") {
    performance.args <- list()
    PM.AHR <- function(y, phat, ties = FALSE, ...) {
      AHR(y = y, phat = phat, ties = ties)
    }
    return(list(AHR, performance.args))
  }

  # IE -------------------------------------------------------------------------
  if (performance == "IE") {

    if (is.null(performance.args$cutoff)) {
      performance.args <- list(cutoff = 300)  # no other possible arguments
    }
    PM.IE <- function(y, phat, ...) {
      IE(y = y, phat = phat,
         cutoff = performance.args$cutoff)  # for now cutoff is at 300 default
    }
    return(list(IE, performance.args))
  }

  # TOP1 -----------------------------------------------------------------------
  if (performance == "TOP1") {
    performance.args <- list()
    PM.TOP1 <- function(y, phat, ...) {
      TOP1(y = y, phat = phat, ...)
    }
    return(list(TOP1, performance.args))
  }

  # RKL ------------------------------------------------------------------------
  if (performance == "RKL") {
    performance.args <- list()
    PM.RKL <- function(y, phat, ...) {
      RKL(y = y, phat = phat, ...)
    }
    return(list(RKL, performance.args))
  }

}

#' Summarising an "\code{epx}" object
#'
#' \code{summary} method for class "\code{\link{epx}}".
#'
#' @param object Object of class "\code{epx}"; a result to a call to
#'   \code{\link{epx}}.
#' @param ... Further arguments passed to or from other methods.
#' @return Prints a summary of the main results of the phalanx-formation
#'   algorithm performed by \code{\link{epx}}.
#' @examples
#' # Example with data(harvest)
#' \donttest{
#' 
#' ## Phalanx-formation using a base classifier with 500 trees (default)
#' 
#' set.seed(761)
#' model <- epx(x = harvest[, -4], y = harvest[, 4],
#'             classifier.args = list(ntree = 500))
#' summary(model)
#'
#' ## The summary corresponds with
#' (model$PHALANXES)[[4]]
#' }
#' @export
summary.epx <- function(object, ...) {

  data <- cbind(object$X, y = object$Y)
  groups <- (object$PHALANXES)[[4]]
  var.names <- colnames(object$X)

  AHR <- object$PHALANXES.FINAL.PERFORMANCE  # AHR is default

  nvars <- ncol(data) - 1
  sort.unique.groups <- 1:max(groups)

  cat("Phalanx-formation algorithm starts with", nvars, "variable(s)", "\n")
  cat("and ends with", length(var.names[groups > 0]),
      "variable(s) grouped into",
      length(sort.unique.groups), "phalanxes:", "\n")
  cat(sort.unique.groups, "\n")

  i <- 0
  for (j in sort.unique.groups) {
    cat("------", "\n")
    i <- i + 1
    namef <- NULL
    namef <- var.names[groups == j]
    cat("Phalanx", j, "contains", length(namef), "variable(s):", "\n")
    cat(namef, sep = ", ")
    cat("\n")
    cat("=>", paste((object$PERFORMANCE.ARGS)$performance.measure),
        "is", AHR[i], "\n")
  }
}



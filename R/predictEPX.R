#' Predict with an "\code{epx}" object
#'
#' Predicted values based on an "\code{\link{epx}}" object; may specify
#' different base classifier arguments than those used for phalanx-formation.
#'
#' @param object Object of class "\code{\link{epx}}".
#' @param newdata An optional data frame specifiying variables with which to
#'   predict; if omitted and \code{classifier.args} are not specified, the
#'   fitted (ensembled) values are used.
#' @param classifier.args Additional arguments for the base classifier; same
#'   base classifier as that used for phalanx-formation (specified in
#'   \code{\link{epx}}).
#' @param ... Further arguments passed to or from other methods.
#' @return Numeric vector of predicted values (double).
#' @examples
#' # Example with data(harvest)
#'
#' ## Phalanx-formation using a base classifier with 500 trees (default)
#' \donttest{
#' set.seed(761)
#' model <- epx(x = harvest[, -4], y = harvest[, 4],
#'              classifier.args = list(ntree = 500))
#'
#' ## Predict training values without additional classifier.args and newdata
#' ## returns the object's ENSEMBLED.FITS
#' all.equal(predict(model), model$ENSEMBLED.FITS)
#'
#' ## Predict training values using 500 trees (default = 500)
#' set.seed(761)
#' preds500 <- predict(model, classifier.args = list(ntree = 500))
#'
#' ## Predict test values by passing dataframe of test predictors to newdata as
#' ## with the predict(model, newdata = . ) function etc.
#' }
#' @export
predict.epx <- function(object, newdata,
                        classifier.args = list(),
                        ...) {

  # get base classifier to use for the emsemble
  FUNS <- .getBaseClassifier(classifier = (object$BASE.CLASSIFIER.ARGS)[[1]],
                             classifier.args)
  BC <- FUNS[[1]]
  BC.predict <- FUNS[[2]]
  classifier.args <- FUNS[[3]]

  ## clarifying what arguments used for the base classifier in predict
  ## vs. what was specified when creating the epx object
  cat("Base classifier:", (object$BASE.CLASSIFIER.ARGS)[[1]], "\n")

  epx.classifier.args <- (object$BASE.CLASSIFIER.ARGS)[[2]]
  cat("Base classifier arguments specified in phalanx-formation:")
  if (length(epx.classifier.args) == 0) {  # no user args in epx
    cat(" none", "\n")
  } else {  # there are user args in epx
    cat("\n")
    print(epx.classifier.args)
  }

  cat("Base classifier arguments specified in prediction:")
  if (length(classifier.args) == 0) {  # no user args from predict
    cat(" none", "\n")
  } else {  # there are user args in predict
    cat("\n")
    print(classifier.args)
  }

  # when newdata is missing values
  if (missing(newdata) && length(classifier.args) == 0) {
    return(object$ENSEMBLED.FITS)
  }

  if (missing(newdata) && length(classifier.args) > 0) {
    # set newdata to be the existing data
    newdata <- object$X
  }

  # have actual newdata, make sure it's okay
  if (ncol(newdata) != ncol(object$X)) {
    stop("Given newdata with incorrect number of columns.")
  }
  #  if (!(colnames(newdata) %in% colnames(object$X)) ||
  #    !(colnames(object$X) %in% colnames(newdata))) {
  if (!all(colnames(newdata) %in% colnames(object$X)) ||
      !all(colnames(object$X) %in% colnames(newdata))) {
    stop("Given newdata with incorrect variable names.")
  }

  # ensemble and fit values ====================================================
  groups <- (object$PHALANXES)[[4]] # to which phalanx each variable belongs
  var.names <- colnames(object$X)
  data <- cbind(object$X, y = object$Y)

  sort.unique.groups <- 1:max(groups)  # of course skip 0th (nothing) phalanx

  # fitting BASE.CLASSIFIER for each phalanx
  predictions <- NULL

  for (j in sort.unique.groups) {
    namef <- NULL
    namef <- var.names[groups == j]
    classifier.j <- BC(.ClassifierFormula(namef), dat = data)
    preds1 <- BC.predict(model = classifier.j,
                         newdata = newdata)

    predictions <- cbind(predictions, preds1)
  }

  colnames(predictions) <- sort.unique.groups
  predicted.values <- as.numeric(apply(predictions, 1, mean))

  return(predicted.values)
}

#' Balanced K-fold cross-validation for an "\code{epx}" object
#'
#' Balanced K-fold cross-validation based on an "\code{\link{epx}}" object.
#' Hence, we have biased cross-validation as we do not re-run the
#' phalanx-formation algorithm for each fold.
#'
#' @param epx Object of class "\code{\link{epx}}".
#' @param folds Optional vector specifying to which fold each observation is to
#'   be assigned. Must be an \eqn{n}-length vector (\eqn{n} being the number of
#'   observations) with integer values only in the range from 1 to \eqn{K}.
#' @param K K-fold cross-validation; default is 10.
#' @param folds.out Indicates whether a vector indicating fold membership for
#'   each of the observations will be output; default is \code{FALSE}.
#' @param classifier.args Arguments for the base classifier specified by
#'   \code{epx}; default is that used in \code{epx} formation.
#' @param performance.args Arguments for the performance measure specified by
#'   \code{epx}; default is that used in \code{epx} formation.
#' @param ... Further arguments passed to or from other methods.
#' @return An \eqn{(n + 1)} by \eqn{(p + 1)} matrix, where \eqn{n} is the number
#'   of observations used to train \code{epx} and \eqn{p} is the number of
#'   (final) phalanxes. The \eqn{p + 1} column of the matrix is the predicted
#'   probabilities of relevance from the ensemble of phalanxes, and the \eqn{n +
#'   1} row is the performance (choice of performance measure determined by the
#'   "\code{\link{epx}}" object) of the corresponding column.
#'
#'   Setting \code{folds.out} as \code{TRUE} changes the output of \code{cv.epx}
#'   into a list of two elements:
#'   \item{EPX.CV}{The \eqn{(n + 1)} by \eqn{(p + 1)} matrix attained by
#'   default when \code{folds.out = FALSE}.}
#'   \item{FOLDS.USED}{A vector of length \eqn{n} with integer values only in
#'   the range from 1 to \code{K} indicating to which fold each observation was
#'   shuffled for cross-validation.}
#' @examples
#' # Example with data(harvest)
#'
#' ## Phalanx-formation using a base classifier with 50 trees (default = 500)
#' \donttest{ 
#' set.seed(761)
#' model <- epx(x = harvest[, -4], y = harvest[, 4],
#'             classifier.args = list(ntree = 50))
#'
#' ## 10-fold balanced cross-validation (different base classifier settings)
#' \dontrun{
#' set.seed(761)
#' cv.100 <- cv.epx(model, classifier.args = list(ntree = 100))
#' tail(cv.100) # see performance (here, AHR) for all phalanxes and the ensemble
#' 
#' 
#' ## Option to output the vector assigning observations to the K folds
#' ## (Commented out for speed.)
#' set.seed(761)
#' cv.folds <- cv.epx(model, folds.out = TRUE)
#' tail(cv.folds[[1]])  # same as first example
#' table(cv.folds[[2]])  # number of observations in each of the 10 folds
#'
#' ## 10 runs of 10-fold balanced cross-validation (using default settings)
#' set.seed(761)
#' cv.ahr <- NULL  # store AHR of each ensemble
#' for (i in 1:10) {
#'   cv.i <- cv.epx(model)
#'   cv.ahr <- c(cv.ahr, cv.i[nrow(cv.i), ncol(cv.i)])
#' }
#' boxplot(cv.ahr)  # to see variation in AHR
#' }
#' }
#' @export
cv.epx <- function(epx,
                   folds = NULL,
                   K = 10,
                   folds.out = FALSE,
                   classifier.args = list(),
                   performance.args = list(),
                   ...) {

  var.names <- colnames(epx$X)
  phalanxes4 <- (epx$PHALANXES)[[4]]
  sort.unique.groups <- 1:max(phalanxes4)
  num.groups <- length(sort.unique.groups)

  # get base classifier and performance measure
  FUNS <- .getBaseClassifier((epx$BASE.CLASSIFIER.ARGS)[[1]], classifier.args)
  BC <- FUNS[[1]]
  BC.PREDICT <- FUNS[[2]]
  classifier.args <- FUNS[[3]]

  performance <- (epx$PERFORMANCE.ARGS)[[1]]
  FUNS.PM <- .getPerformanceMeasure(performance, performance.args)
  PM <- FUNS.PM[[1]]
  performance.args <- FUNS.PM[[2]]

  ## clarifying what arguments used for the base classifier in predict
  ## vs. what was specified when creating the epx object
  cat("Base classifier:", (epx$BASE.CLASSIFIER.ARGS)[[1]], "\n")

  epx.classifier.args <- (epx$BASE.CLASSIFIER.ARGS)[[2]]
  cat("Base classifier arguments specified in phalanx-formation:")
  if (length(epx.classifier.args) == 0) {  # no user args in epx
    cat(" none", "\n")
  } else {  # there are user args in epx
    cat("\n")
    print(epx.classifier.args)
  }

  cat("Base classifier arguments specified in balanced")
  cat(" ", K, "-fold cross-validation:", sep = "")
  if (length(classifier.args) == 0) {  # no user args from cv
    cat(" none", "\n")
  } else {  # there are user args in cv
    cat("\n")
    print(classifier.args)
  }

  # balanced cross-validation starts here ######################################

  data <- cbind(epx$X, y = epx$Y)  # response will for sure be last column
  n <- nrow(data)
  posn <- ncol(data)

  ## setup the folds
  s <- vector(mode = "numeric", length = n)  # fold membership vector
  if (is.null(folds)) {
    nact <- c(1:n)[data[,posn] == 1]  ## index of relevants
    niact <- c(1:n)[data[,posn] == 0]  ## index of irrelevants
    fact <- ceiling(length(nact)/K)  # num. of relevants/fold
    sact <- sample(rep(1:K, fact), length(nact))  # distributing relevants
    fiact <- ceiling(length(niact)/K)  # num. of irrelevants/fold
    siact <- sample(rep(1:K, fiact), length(niact))  # distributing irrelevants
    s[nact] <- sact  # s now holds fold membership for all observations
    s[niact] <- siact
  } else {
    if ( (mean(unique(folds) %in% 1:K) != 1) ||
         (1:K %in% mean(unique(folds)) != 1) ) {
      stop("User-defined folds are not well-defined.")
    }
    s <- folds
  }

  ms <- max(s)  # number of folds

  # matrix FP contains the probability of being relevant for each random forest
  # each column is the predictions for all K folds from a group
  # every n/K rows are all the predictions for a fold
  FP <- matrix(nrow = n, ncol = num.groups)
  RFP <- NULL

  for(i in 1:ms) {

    j.in <- sort(c(1:n)[(s != i)])
    j.out <- sort(c(1:n)[(s == i)])

    data.train <- data[j.in,]
    data.test <- data[j.out,]
    data.train <- as.data.frame(data.train)
    data.test <- as.data.frame(data.test)

    # training-test split is done for a fold
    k <- 0

    # fitting classifier for each clusters relating data
    for(j in sort.unique.groups) {

      k <- k + 1
      namef <- NULL
      namef.ind <- which(phalanxes4 == j)
      namef <- var.names[namef.ind]
      classifier.j <- BC(.ClassifierFormula(namef), dat = data.train)
      FP[j.out, k] <- BC.PREDICT(model = classifier.j, newdata = data.test)

    }
    # training-test split for a fold ends here

  }

  # cross-validation ends here #################################################

  if(num.groups > 1) {
    FPUF <- apply(FP, 1, mean)  # FPUF = ensembled predictions
  } else {  # (num.groups == 1)
    FPUF <- FPFL <- FP
  }

  preds <- cbind(FP, FPUF)
  # last col of preds is for the averaged ensemble filtered

  # labelling preds
  colnames(preds) <- c(paste(sort.unique.groups), "ensemble")

  # computing performance for each classifier/column
  leap <- dim(preds)[[2]]
  avg <- rep(NA, leap)
  for(j in 1:leap){
    avg[j] <- PM(data[,posn], preds[,j], ties = TRUE)
  }

  res <- rbind(preds, performance = avg)

  ## the first (n-1) rows of res contains the probability of being relevant
  ## the last row of res contains the performance for each classifier (phalanx)
  ## last col of preds is for the averaged ensemble

  if (folds.out) {
    return(list(EPX.CV = res,
                FOLDS.USED = s))
  } else {
    return(res)
  }

}


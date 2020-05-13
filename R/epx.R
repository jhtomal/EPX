#' Fitting an Ensemble of Phalanxes
#'
#' \code{epx} forms phalanxes of variables from a given training dataset where
#' the response is a binary response variable with a rare class. Each of these
#' disjoint subsets of variables is then fit with a base classifier together
#' which forms the emsemble. Typically performance will then be assessed with
#' \code{\link{cv.epx}}. Please see Tomal et al. (2015) for more details.
#'
#' @references
#' Tomal, J. H., Welch, W. J., & Zamar, R. H. (2015).
#' Ensembling classification models based on phalanxes of variables with
#' applications in drug discovery.
#' \emph{The Annals of Applied Statistics},
#' \emph{9}(1), 69-93.
#' \url{http://doi.org/10.1214/14-AOAS778}
#'
#' @param x Explanatory variables contained in a data frame.
#' @param y Binary response variable vector (numeric or integer);
#'   1 is the rare class, 0 is non-rare.
#' @param phalanxes.initial Initial variable grouping indices; default is no
#'   grouping. Example: define as vector c(1, 1, 2, 2, 3, ...); note that cannot
#'   indices cannot be skipped i.e. c( 1, 3, 3, 4, 4, 3, 1) is invalid.
#' @param iquant Quantile index for the distribution of the performance measure
#'   (default is average hit rate) of a classifier that does random ranking
#'   (i.e. predictors have no explanatory power at all); default is 0.95.
#' @param nsim Number of simulations done to get the reference distribution of
#'   the performance measure; default is 1000.
#' @param rmin.target Minimum number of variables in each phalanx; default is 1.
#' @param classifier Choice of base classifier from one of:
#'   \code{c("random forest", "logistic regression", "neural network")};
#'   default is "random forest", which uses
#'   \code{\link[randomForest]{randomForest}}.
#' @param classifier.args Arguments relating to the base \code{classifier}
#'   specified in a list as follows: \code{list(argName1 = value1, argName2 =
#'   value2, ...)}. If left empty, the classifier will use whatever its defaults
#'   happen to be. For "random forest", user may specify \code{replace, cutoff,
#'   nodesize, maxnodes}. For "logistic regression" there are no options. For
#'   "neural network", user may specify \code{size, trace}.
#' @param performance Choice of performance measure; default is
#'   \code{\link{AHR}}. Must be one of:
#'   \code{c("\link{AHR}", "\link{IE}", "\link{TOP1}", "\link{RKL}")}.
#' @param performance.args Arguments relating to the \code{performance}
#'   specified in a list as follows:
#'   \code{list(argName1 = value1, argName2 = value2, ...)}.
#'   If left empty, the performance will use whatever its
#'   defaults happens to be. User may only specify \code{cutoff} for
#'   \code{\link{IE}} since all other performance measures have no additional
#'   parameters.
#' @param computing Whether to compute sequentially or in parallel. Input is one
#'   of \code{c("sequential", "parallel")}; default is "sequential".
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The \code{epx} function returns an object of class "\code{epx}". The
#'   function \code{\link{summary.epx}} prints a summary of the results. An
#'   object of class "\code{epx}" is a list containing the following components:
#' \item{PHALANXES}{List of four vectors, each the same length as the number of
#' explanatory variables (columns in \code{X}): \code{phalanxes.initial},
#' \code{phalanxes.filtered}, \code{phalanxes.merged}, \code{phalanxes.final}.
#' Each vector records the phalanx membership of all explanatory variables in
#' the four stages of phalanx-formation. The \eqn{i}th entry of a vector
#' indicates to which phalanx the \eqn{i}th variable belongs. Phalanx "0" does
#' not exist and so membership to phalanx "0" means that the variable does not
#' belong to any phalanx.}
#' \item{PHALANXES.FINAL.PERFORMANCE}{Vector with the performance measure,
#' according to the specified \code{performance}, of each of the final
#' phalanxes. The indices match the phalanx number.}
#' \item{PHALANXES.FINAL.FITS}{A matrix with number of rows equal to the number
#' of observations in the training data and number of columns equal to the
#' number of final phalanxes. The \eqn{i}th column is the predicted
#' probabilities of relevance of each observation produced by fitting the base
#' classifier, specified by \code{classifier}, to the \eqn{i}th phalanx's
#' variables.}
#' \item{ENSEMBLED.FITS}{The predicted probabilities of relevance from the
#' ensemble of phalanxes based on \code{phalanxes.final}.}
#' \item{BASE.CLASSIFIER.ARGS}{(Parsed) record of user-specified arguments for
#' \code{classifier}.}
#' \item{PERFORMANCE.ARGS}{(Parsed) record of user-specified arguments for
#' \code{performance}.}
#' \item{X}{User-provided data frame of explanatory variables.}
#' \item{Y}{User-provided binary response vector.}
#' @examples
#' # Example with data(harvest)
#'
#' ## Phalanx-formation using a base classifier with 500 trees (default)
#' set.seed(761)
#' model <- epx(x = harvest[, -4], y = harvest[, 4],
#'              classifier.args = list(ntree = 500))
#'
#' ## See how phalanx-membership of variables at all four steps
#' ## (0 means not in a phalanx.)
#' (model$PHALANXES)
#'
#' ## Summary of how predictors divided into the final phalanxes (matches above)
#' summary(model)
#'
#' ## Parallel computing
#' ## (Commented out for speed.)
#' # clusters <- parallel::detectCores()
#' # cl <- parallel::makeCluster(clusters)
#' # doSNOW::registerDoSNOW(cl)
#' # set.seed(761)
#' # model.par <- epx(x = harvest[, -4], y = harvest[, 4],
#' #                  computing = "parallel")
#' # parallel::stopCluster(cl)
#' @export
#' @exportClass epx
epx <- function(x, y,
                phalanxes.initial = c(1:ncol(x)),
                iquant = 0.95,
                nsim = 1000,
                rmin.target = 1,
                classifier = "random forest",
                classifier.args = list(),
                performance = "AHR",
                performance.args = list(),
                computing = "sequential",
                ...) {

  p <- ncol(x)

  ## Error handling ============================================================

  if (missing(x) || missing(y)) {
    stop("Arguments 'x', and 'y' are missing.")
  }

  # y is passed to epx_algorithm as a numeric or integer vector with
  # values 0 or 1 because AHR requires a numeric vector and converting a
  # binary factor back to a numeric won't give back 0 and 1s
  # (I keep getting 1s and 2s). I could
  # make it so that the input is a factor and change it to numeric and then
  # correct the coversion results to 0s and 1s, but making the user pass the
  # appropriate 0-1 numeric vector means easier and cleaner code as I can
  # convert it into a factor whenever I want.
  # i.e. y will be a factor when used in random forest due to .ClassifierFormula
  # making y a factor as required.
  if (!is.numeric(y) & !is.integer(y)) {
    stop("'y' must be a numeric or integer vector.")
  }

  if ( !all(y %in% c(0, 1)) || !all(c(0, 1) %in% y) ) {
    stop("'y' must be a binary response vector; 1 is the rare class, 0 is non-rare. ")
  }

  if (class(x) != "data.frame") {
    stop("'x' is not a dataframe.")
  }

  if (is.null(nrow(x))) {
    stop("'x' cannot be empty.")
  }

  if (nrow(x) == 0L) {
    stop("0 (non-NA) cases")
  }

  if (any(!(phalanxes.initial %in% 0:length(phalanxes.initial)))) {
    stop("Initial phalaxes not well-defined.")
  }

  if (length(phalanxes.initial) != ncol(x)) {
    stop("Initial phalanxes are not defined for all variables.")
  }

  if (length(unique(phalanxes.initial)) != max(phalanxes.initial)) {
    stop("Initial phalanxes are not well-defined.")
  }

  if (iquant > 1 || iquant <= 0) {
    stop("Invalid reference distribution quantile.")
  }
  if (nsim <= 0) {
    stop("Invalid number of simulations.")
  }
  if (rmin.target >= ncol(x) || rmin.target < 1) {
    stop("Invalid minimum phalanx size.")
  }

  ## Parallel computing ========================================================
  computing <- match.arg(arg = computing,
                         choices = c("sequential", "parallel"),
                         several.ok = FALSE)

  if (computing == "sequential") { # suppress irritating warning message for seq
    foreach::registerDoSEQ()
  } else if (foreach::getDoParWorkers() == 1) {
    cat("Selected parallel computing, but only 1 execution worker
        in the currently registered doPar backend.")
  }

  ## epx algoritm ==============================================================
  RES <- epxAlgorithm(x = x,
                      y = y,
                      phalanxes.initial = phalanxes.initial,
                      iquant = iquant,
                      nsim = nsim,
                      rmin.target = rmin.target,
                      classifier = classifier,
                      classifier.args = classifier.args,
                      performance = performance,
                      performance.args = performance.args,
                      ...)

  ## Return epx ================================================================
  return(RES)
  }

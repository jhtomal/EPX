#' @import randomForest
#' @import nnet

.getBaseClassifier <- function(classifier, classifier.args, ...) {

  classifier <- match.arg(arg = classifier,
                          choices = c("random forest", "logistic regression",
                                      "neural network", "nnet"),
                          several.ok = FALSE)

  # random forest --------------------------------------------------------------
  if (classifier == "random forest") {
    # make sure no stupid arguments allowed
    if (length(classifier.args) > 0) {
      valid.args <- match.arg(arg = names(classifier.args),
                              choices = c("ntree", "replace", "cutoff",
                                          "nodesize", "maxnodes"),
                              several.ok = TRUE)
      classifier.args <- classifier.args[valid.args]
    }
    classifier.args.def <- classifier.args

    if (is.null(classifier.args$ntree)) {classifier.args.def$ntree <- 500}
    if (is.null(classifier.args$replace)) {classifier.args.def$replace <- TRUE}
    if (is.null(classifier.args$cutoff)) {classifier.args.def$cutoff <- c(1/2, 1/2)}
    if (is.null(classifier.args$nodesize)) {classifier.args.def$nodesize <- 1}
    if (is.null(classifier.args$maxnodes)) {classifier.args.def$maxnodes <- NULL}

    BC <- function(formula, dat) {
      return(randomForest::randomForest(formula,
                                        data = dat,
                                        ntree = classifier.args.def$ntree,
                                        replace = classifier.args.def$replace,
                                        cutoff = classifier.args.def$cutoff,
                                        nodesize = classifier.args.def$nodesize,
                                        maxnodes = classifier.args.def$maxnodes,
                                        keep.forest = TRUE))
    }

    BC.predict <- function(model, newdata) {
      return(predict(model, newdata, type = "prob")[,2])
    }
  }

  # logistic regression --------------------------------------------------------
  if (classifier == "logistic regression") {
    # no options for logistic regression
    classifier.args <- list()

    BC <- function(formula, dat) {
      return(glm(formula,
                 family = "binomial",
                 data = dat))
    }
    BC.predict <- function(model, newdata) {
      return(predict(model, newdata, type = "response"))
    }
  }

  # neural network -------------------------------------------------------------
  if (classifier == "neural network") {
    # ...: size
    if (length(classifier.args) > 0) {
      valid.args <- match.arg(arg = names(classifier.args),
                              choices = c("size", "trace"),
                              several.ok = TRUE)
      classifier.args <- classifier.args[valid.args]
    }
    classifier.args.def <- classifier.args

    if (is.null(classifier.args$size)) {classifier.args.def$size <- 1}
    if (is.null(classifier.args$trace)) {classifier.args.def$trace <- FALSE}
    BC <- function(formula, dat) {
      return(nnet::nnet(formula,
                        size = classifier.args.def$size,
                        data = dat,
                        trace = classifier.args.def$trace))
    }
    BC.predict <- function(model, newdata) {
      return(predict(model, newdata, type = "raw"))
    }
  }

  return(list(BC, BC.predict, classifier.args))

}

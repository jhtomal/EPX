# Helpers used by the package.
#' @import graphics
#' @import stats

################################################################################
# .ClassifierFormula ###########################################################
################################################################################
# Gives formula necessary for many base classifier functions.
#
# Args:
#   explanatory.vars: vector of explanatory variables not as characters,
#                     directly what you would get from colnames(data).
#
# Returns:
#   formula that can be passed into randomForest(), lm(), etc.
.ClassifierFormula <- function(explanatory.vars) {



  return(as.formula(paste("factor(y) ~ ",
                          paste(explanatory.vars, collapse=" + "),
                          sep = "")))
}

################################################################################
# get.nhit.all #################################################################
################################################################################
# Arguments:
#   y: response vector for the entire data
#   phat.predict: predictive hat matrix obtained from get.phat.cv
#   max.cutoff: the number of the compounds selected, e.g., 300
#
# Return values:
#   nhit: number of hits for each model
get.nhit.all <- function(y, phat.predict, max.cutoff) {

  n.model <- nrow(phat.predict)
  nhit <- rep(0,n.model)

  for (i in 1:n.model) {
    nhit[i] <- hit.curve(phat = phat.predict[i,],
                         y = y,
                         max.cutoff,
                         plot.hc = F)$nhitlast
  }

  return(nhit)

}

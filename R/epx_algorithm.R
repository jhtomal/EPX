# Main helper for \link{epx}. For parameter definitions, see \link{epx}
#' @import doSNOW
#' @import rngtools
#' @import foreach
#' @import doRNG
epxAlgorithm <- function(x,
                         y,
                         phalanxes.initial,
                         iquant,
                         nsim,
                         rmin.target,
                         classifier, classifier.args,
                         performance, performance.args,
                         ...) {



  ## Start epx algorithm =======================================================
  # y is passed to this function as a numeric vector with values 0 or 1
  # y will be a factor when used in random forest due to .ClassifierFormula
  # making y a factor as required.

  data <- cbind(x, y)

  n    <- nrow(data)
  posn <- ncol(data)  # num. of columns; posn is where y col. is in data
  p <- ncol(x)

  ## Get base classifier, corresponding prediction fcn, perf fcn ===============
  FUNS <- .getBaseClassifier(classifier, classifier.args)
  BC <- FUNS[[1]]
  BC.PREDICT <- FUNS[[2]]
  classifier.args <- FUNS[[3]]

  FUNS.PM <- .getPerformanceMeasure(performance, performance.args)
  PM <- FUNS.PM[[1]]
  performance.args <- FUNS.PM[[2]]

  ## So that variables in foreach etc. are "globally bound".... ================
  tdata <- NULL
  trfprob <- NULL
  r <- NULL
  j <- NULL

  ## Console output for classifier and performance specifications ==============
  cat("Performance measure:", performance, "\n")
  cat("Performance measure additional arguments:")
  if (length(performance.args) == 0) {
    cat(" none", "\n\n")
  } else {
    cat("\n")
    print(performance.args)
  }

  cat("Base classifier:", classifier, "\n")
  cat("Base classifier arguments specified in phalanx-formation:")
  if (length(classifier.args) == 0) {
    cat(" none", "\n\n")
  } else {
    cat("\n")
    print(classifier.args)
  }

  # ============================================================================
  # Get the preferred quantile (qsim) for the distribution of the performance
  # measure of a random classifier
  qsim <- foreach(i = 1:nsim, .combine = "c") %dorng% {
    pr <- seq(0, 1, length = n)
    pr1 <- sample(pr, n, replace = FALSE)
    PM(y, pr1, ties = FALSE)
  }  # the empirical reference distribution
  qmean <- quantile(qsim, prob = 0.50)
  qsim <- quantile(qsim, prob = iquant)

  cat("Phalanx formation is in progress, please wait ……", "\n")

  # cat("Reference distribution quantiles:", "\n")
  # cat("q", iquant, " = ", qsim, "\n", sep = "")
  # cat("q0.50 =", qmean, "\n\n")

  # STEP 1 INITIAL GROUPING DONE ###############################################
  step1phalanxes <- as.numeric(phalanxes.initial)  # 0 means belongs in no phalanx
  var.names <- colnames(x)  # names of predictors

  # Default performance measure is average hit rate (AHR), so from now on refer
  # to the performance measure as AHR.

  # Get the matrix of the probability of activity (rfprob)
  # and the vector of AHRs (sahr) of all single classifiers #
  rfprob <- foreach(j = 1:max(step1phalanxes),
                    .combine = "cbind",
                    .packages = c("randomForest", "nnet")) %dorng% {
                      inds <- which(step1phalanxes %in% j)
                      namef <- var.names[inds]
                      rf <- BC(formula = .ClassifierFormula(namef), dat = data)
                      pi <- BC.PREDICT(model = rf)
                      # in case RF returns NA for some obs. (very sparse 1s)
                      pi <- replace(pi, list = is.na(pi), 0)
                      PMj <- PM(y, pi, ties = TRUE)
                      c(pi, PMj)
                    }
  sahr <- rfprob[(n+1),] # last row is performances
  rfprob <- rfprob[-(n+1),] # matrix of probabilities
  trfprob <- apply(rfprob, 1, mean)

  # Get the matrix (matsyn) of AHRs when
  #   pairs of variables are merged together,
  # and the matrix (matcom) of AHRs when
  #   pairs of variables are averaged or ensembled
  num.phalanxes <- max(step1phalanxes)
  matsyn <- matcom <- matrix(NA, nrow = num.phalanxes, ncol = num.phalanxes)
  ahrsyn <- ahrcom <- NULL

  # generate sequence of seeds of length the number of computations
  rng <- rngtools::RNGseq((num.phalanxes^2 - num.phalanxes)/2, 761)

  metrics <- foreach(i = 1:(num.phalanxes - 1),
                     .combine = "cbind", .inorder = T) %:%
    foreach(j = (i + 1):num.phalanxes,
            r = rng[num.phalanxes*(i-1) - ((i-1)^2 - (i-1))/2
                    + 1:(num.phalanxes-1)],
            .combine = "cbind", .inorder = T,
            .packages = c("randomForest", "nnet")) %dopar% {
              # set RNG seed
              rngtools::setRNG(r)

              inds <- which(step1phalanxes %in% c(i, j))
              namef <- var.names[inds]
              rf <- BC(formula = .ClassifierFormula(namef), dat = data)

              pij <- BC.PREDICT(model = rf)
              pij <- replace(pij, list = is.na(pij), 0)
              ahrij <- PM(y, pij, ties = TRUE)

              # gotcha: rfprob has fits for each phalanx, so to get pijbar,
              # we avg. across rows for columns == c(i, j) i.e. NOT INDS!!
              pijbar <- apply(rfprob[, c(i, j)], 1, mean)
              ahrijbar <- PM(y, pijbar, ties = TRUE)

              rbind(ahrij, ahrijbar)
            }

  if(is.matrix(metrics)){
    ahrsyn <- metrics[1,]
    ahrcom <- metrics[2,]
  } else {
    ahrsyn <- metrics[1]
    ahrcom <- metrics[2]
  }

  matsyn[lower.tri(matsyn)] <- ahrsyn
  matsyn <- t(matsyn)
  matsyn[lower.tri(matsyn)] <- ahrsyn

  matcom[lower.tri(matcom)] <- ahrcom
  matcom <- t(matcom)
  matcom[lower.tri(matcom)] <- ahrcom

  # STEP 2 COMPUTATION DONE ####################################################

  # groups of variables to keep (CKEEP) after first filtering (step 2)
  smatsyn <- sweep(matsyn, 2, sahr - qmean)
  smatcom <- sweep(matcom, 2, sahr - qmean)

  CKEEP <- NULL
  step2phalanxes <- rep(0, p)
  j <- 1
  for (i in 1:num.phalanxes) {
    if(max(c(sahr[i], smatsyn[i,], smatcom[i,]),
           na.rm = TRUE) >= qsim) {
      CKEEP <- c(CKEEP, i)
      step2phalanxes[step1phalanxes == i] <- j
      j <- j + 1
    }
  }

  if (length(CKEEP) == 0) {
    stop("All initial phalanxes filtered out.
         Consider alternate initial phalanx groupings.")
  }

  # cat("Number of deleted phalanxes after first filtering:",
  #    num.phalanxes - length(CKEEP),
  #    "out of",
  #    num.phalanxes, "\n")

  # STEP 2 FILTERING DONE ######################################################

  # update sahr, rfprob, matsyn and matcom
  sahr <- sahr[CKEEP]
  rfprob <- rfprob[,CKEEP]

  if (length(CKEEP) == 1) {
    matsyn <- matsyn[CKEEP, ]
    matsyn <- matsyn[CKEEP]

    matcom <- matcom[CKEEP, ]
    matcom <- matcom[CKEEP]
  } else {
    matsyn <- matsyn[CKEEP, ]
    matsyn <- matsyn[, CKEEP]

    matcom <- matcom[CKEEP, ]
    matcom <- matcom[, CKEEP]
  }

  num.phalanxes <- max(step2phalanxes)

  # STEP 3 PREPARATION DONE ####################################################

  # STEP 3 STARTS HERE #########################################################
  # hierarchical merging of phalanxes
  step3phalanxes <- step2phalanxes
  rmin <- suppressWarnings( min((matcom / matsyn), na.rm = TRUE) )  # mij

  ## in response to warning: "no non-missing arguments to min; returning Inf"
  if (is.infinite(rmin)) {
    rmin <- rmin.target + 1
  }

  # need to set seeds here because of parallel computing -_-
  # (not sure why)
  seeds <- sample(1:p^2, 3*p, replace = FALSE)
  k <- 1

  # start of while loop (i.e. while AHR(ensemble) < AHR(merged), merge groups)
  while ((rmin < rmin.target) && (num.phalanxes > 1)) {

    set.seed(seeds[k])
    k <- k + 1

    index.rmin <- which((matcom / matsyn) == rmin, arr.ind = TRUE)[1,]
    min.index.rmin <- min(index.rmin)
    max.index.rmin <- max(index.rmin)

    # updating (sahr, rfprob, matsyn, matcom, groups) starts here
    sahr <- sahr[-max.index.rmin]
    rfprob <- rfprob[,-max.index.rmin]
    if (num.phalanxes > 2) {
      matsyn <- matsyn[-max.index.rmin, ] ; matsyn <- matsyn[, -max.index.rmin]
      matcom <- matcom[-max.index.rmin, ] ; matcom <- matcom[, -max.index.rmin]
    } else if(num.phalanxes == 2){
      matsyn <- matsyn[-max.index.rmin, ] ; matsyn <- matsyn[-max.index.rmin]
      matcom <- matcom[-max.index.rmin, ] ; matcom <- matcom[-max.index.rmin]
    }

    # updating phalanx names (i.e. numbers)
    step3phalanxes <-
      replace(step3phalanxes,
              list = which(step3phalanxes == max.index.rmin),
              min.index.rmin)
    merged.px <- min(which(step3phalanxes == min.index.rmin))
    # want to take the min because merged with the max.index.rmin phalanx
    # but we chose to e as the min.index.rmin
    step3phalanxes[step3phalanxes >= max.index.rmin] <-
      step3phalanxes[step3phalanxes >= max.index.rmin] - 1
    num.phalanxes <- max(step3phalanxes)

    merged.px <- step3phalanxes[merged.px]  # the name/number of newly formed px
    # i.e. merged.px is now the index row and column of the matrices that
    # must be recalculated

    # recalculate metric for individual phalanx "min.index.rmin"
    inds <- which(step3phalanxes %in% merged.px)
    namef <- var.names[inds]
    rf <- BC(formula = .ClassifierFormula(namef), dat = data)

    if (num.phalanxes > 1) {
      rfprob[, merged.px] <- BC.PREDICT(model = rf)
      sahr[merged.px] <- PM(data[,posn],
                            rfprob[, merged.px],
                            ties = TRUE)
    } else if (num.phalanxes == 1) {
      rfprob <- BC.PREDICT(model = rf)
      sahr <- PM(data[,posn], rfprob, ties = TRUE)
    }

    # recalculate matcom, matsyn if you have more than 1 phalanx left
    # only for the merged.px colum and row (remember we renamed already!)
    if (num.phalanxes > 1) {

      i <- merged.px
      metrics <-
        foreach(j = 1:num.phalanxes,
                .combine = "cbind", .inorder = T,
                .packages = c("randomForest", "nnet")) %dorng% {
                  inds <- which(step3phalanxes %in% c(i, j))
                  namef <- var.names[inds]
                  rf <- BC(formula = .ClassifierFormula(namef), dat = data)

                  pij <- BC.PREDICT(model = rf)
                  pij <- replace(pij, list = is.na(pij), 0)
                  ahrij <- PM(y, pij, ties = TRUE)

                  pijbar <- apply(rfprob[, c(i,j)], 1, mean)
                  ahrijbar <- PM(y, pijbar, ties = TRUE)

                  rbind(ahrij, ahrijbar)
                }

      metrics[,merged.px] <- NA
      ahrsyn <- metrics[1,]
      ahrcom <- metrics[2,]

      matsyn[merged.px,] <- matsyn[,merged.px] <- ahrsyn
      matcom[merged.px,] <- matcom[,merged.px] <- ahrcom
      # gotcha, rfprob could be a vector, in which case apply doesn't work above

    } else if (num.phalanxes == 1) {

      set.seed(seeds[k])
      k <- k + 1

      inds <- which(step3phalanxes %in% 1)
      namef <- var.names[inds]
      rf <- BC(formula = .ClassifierFormula(namef), dat = data)

      pij <- BC.PREDICT(model = rf)
      ahrij <- PM(y, pij, ties = TRUE)

      pijbar <- mean(rfprob)
      ahrijbar <- PM(y, pijbar, ties = TRUE)

      matsyn <- ahrij
      matcom <- ahrijbar

    }

    # updating (sahr, rfprob, matsyn, matcom, phalanxes) ends here

    if (num.phalanxes > 1) {

      rmin <- suppressWarnings( min((matcom / matsyn), na.rm = TRUE) )  # mij

      ## in response to warning: "no non-missing arguments to min; returning Inf"
      if (is.infinite(rmin)) {
        rmin <- rmin.target + 1
      }

    }
  }

  # end of while loop

  # start of final filtering #

  num.phalanxes <- max(step3phalanxes)
  CKEEP2 <- NULL
  step4phalanxes <- rep(0, p)
  j <- 1
  if (num.phalanxes > 1) {
    smatcom <- sweep(matcom, 2, sahr - qmean)

    for (i in 1:num.phalanxes) {
      if (max(c(sahr[i], smatcom[i,]), na.rm = TRUE) >= qsim) {
        CKEEP2 <- c(CKEEP2, i)
        step4phalanxes[step3phalanxes == i] <- j
        j <- j + 1
      }
    }
  } else {
    cat("No final filtering done since only one phalanx left following merging step.", "\n")
    CKEEP2 <- max(step3phalanxes)  # only 1 phalanx left, no filtering
    step4phalanxes <- step3phalanxes # hence same as step3phalanxes
  }

  # cat("Number of deleted phalanxes after final filtering:",
  #    num.phalanxes - length(CKEEP2),
  #    "out of",
  #    num.phalanxes, "\n\n")

  # end of final filtering

  # WHAT IF ALL GROUPS FILTERED OUT?
  if (length(CKEEP2) == 0) {
    step4phalanxes <- step3phalanxes
    cat("Final filtering step removes all candidate phalanxes. Unfiltered phalanxes returned as final phalanxes.", "\n")
  }

  # getting predicted values etc.
  # gotcha rfprob will be a vector if there's 1 px
  if (length(CKEEP2) > 1) {
    px.predictions <- rfprob[, CKEEP2]
    rownames(px.predictions) <- NULL
    colnames(px.predictions) <- NULL
    fitted.values <- as.numeric(apply(px.predictions, 1, mean))
  } else {
    px.predictions <- as.numeric(rfprob)
    fitted.values <- as.numeric(px.predictions)
  }

  final.perfs <- sahr[CKEEP2]

  # prepare results
  res <- list(PHALANXES = list(phalanxes.initial = step1phalanxes,
                               phalanxes.filtered = step2phalanxes,
                               phalanxes.merged = step3phalanxes,
                               phalanxes.final = step4phalanxes),
              PHALANXES.FINAL.PERFORMANCE = as.numeric(final.perfs),
              PHALANXES.FINAL.FITS = px.predictions,
              ENSEMBLED.FITS = fitted.values,
              BASE.CLASSIFIER.ARGS = list(base.classifier = classifier,
                                          args = classifier.args),
              PERFORMANCE.ARGS = list(performance.measure = performance,
                                      args = performance.args),
              X = x,
              Y = y)

  ## Set the name for the class
  class(res) <- "epx"

  return(res)
  }

# --------------------------------------------------
# Print methods
# --------------------------------------------------
#' Print the output of a \code{"drtmle"} object.
#'
#' @param x A \code{"drtmle"} object
#' @param ... Other arguments (not used)
#' @export
#' @method print drtmle
#
print.drtmle <- function(x, ...) {
  tmp <- list(
    est = cbind(x$drtmle$est),
    cov = x$drtmle$cov
  )
  row.names(tmp$est) <- x$a_0
  colnames(tmp$est) <- ""
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if (length(x$a_0) <= 4) {
    print(tmp)
  } else {
    tmp$cov <- diag(tmp$cov)
    print(tmp)
  }
  invisible(tmp)
}

#' Print the output of a \code{"adaptive_iptw"} object.
#'
#' @param x A \code{"adaptive_iptw"} object.
#' @param ... Other arguments (not used)
#' @export
#' @method print adaptive_iptw
#
print.adaptive_iptw <- function(x, ...) {
  tmp <- list(
    est = cbind(x$iptw_tmle$est),
    cov = x$iptw_tmle$cov
  )
  row.names(tmp$est) <- x$a_0
  colnames(tmp$est) <- ""
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if (length(x$a_0) <= 4) {
    print(tmp)
  } else {
    tmp$cov <- diag(tmp$cov)
  }
  invisible(tmp)
}


#' Print the output of ci.drtmle
#' @export
#' @param x An object of class ci.drtmle
#' @param digits Number of digits to round to
#' @param ... Other options (not currently used)
#' @method print ci.drtmle
#
print.ci.drtmle <- function(x, digits = 3, ...) {
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
  invisible(tmp)
}

#' Print the output of ci.adaptive_iptw
#' @export
#' @param x An object of class ci.adaptive_iptw
#' @param digits Number of digits to round to
#' @param ... Other options (not currently used)
#' @method print ci.adaptive_iptw
#
print.ci.adaptive_iptw <- function(x, digits = 3, ...) {
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
  invisible(tmp)
}

#' Print the output of wald_test.drtmle
#' @export
#' @param x An object of class wald_test.drtmle
#' @param digits Number of digits to round to
#' @param ... Other options (not currently used)
#' @method print wald_test.drtmle
#
print.wald_test.drtmle <- function(x, digits = 3, ...) {
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
  invisible(tmp)
}

#' Print the output of wald_test.adaptive_iptw
#' @export
#' @param x An object of class wald_test.adaptive_iptw
#' @param digits Number of digits to round to
#' @param ... Other options (not currently used)
#' @method print wald_test.adaptive_iptw
#
print.wald_test.adaptive_iptw <- function(x, digits = 3, ...) {
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
  invisible(tmp)
}


# --------------------------------------------------
# Plot methods
# --------------------------------------------------
#' Plot reduced dimension regression fits
#'
#' @param x An object of class \code{"drtmle"}
#' @param nPoints Number of points to plot lines (increase for less bumpy plot,
#'  decrease for faster evaluation).
#' @param a_0 For what value of a_0 should the plot be made for?
#' @param ask Boolean indicating whether R should ask to show each plot
#' @param ... More arguments passed to \code{plot}
#' @export
#' @method plot drtmle
#' @importFrom graphics axis lines par plot
#' @importFrom stats plogis
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # simulate data
#' set.seed(123456)
#' n <- 100
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n, 1, plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
#' # fit drtmle with maxIter = 1 to run fast
#' \donttest{
#' fit1 <- drtmle(
#'   W = W, A = A, Y = Y, a_0 = c(1, 0),
#'   family = binomial(),
#'   stratify = FALSE,
#'   SL_Q = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#'   SL_g = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#'   SL_Qr = "SL.npreg", SL_gr = "SL.npreg",
#'   maxIter = 1, returnModels = TRUE
#' )
#' # plot the reduced-dimension regression fits (not run)
#' 
#' plot(fit1)
#' }
#' #
plot.drtmle <- function(x, nPoints = 500,
                        ask = TRUE,
                        a_0 = x$a_0[1], ...) {
  # ask to see next plot
  par(ask = ask)
  # check if returnModels is null
  if (is.null(x$QrnMod) & is.null(x$grnMod)) {
    stop("Plot function only works if returnModels = TRUE.")
  }

  # which entry in x fits corresponds to this a_0
  listInd <- which(x$a_0 == a_0)

  # get local nuisance fits
  gn <- x$nuisance_drtmle$gnStar[[listInd]]
  Qn <- x$nuisance_drtmle$QnStar[[listInd]]

  # ------------------
  # plot Qrn fit
  # ------------------
  if (!(class(x$QrnMod) == "NULL")) {
    # number of fits (if no CV = 1, if CV > 1)
    nFit <- length(x$QrnMod)
    # xlim = range of gn
    xl <- range(gn)
    # prediction points
    predP <- seq(xl[1], xl[2], length = nPoints)
    # get predictions back for each Qrn fit
    fit_Qrn <- lapply(x$QrnMod, function(y) {
      newDat <- data.frame(gn = predP)
      if ("glm" %in% class(y[[listInd]])) {
        predict(y[[listInd]], newdata = newDat, type = "response")
      } else if (class(y[[listInd]]) == "SuperLearner") {
        pred <- predict(y[[listInd]], newdata = newDat)
        # get sl prediction if meta learning did not fail
        if (!all(y[[listInd]]$coef == 0)) {
          pred$pred
          # otherwise get discrete super learner
        } else {
          pred$library.predict[, which.min(y$cvRisk)]
        }
      } else {
        predict(y[[listInd]]$fit, newdata = newDat, type = "response")
      }
    })
    # get ylimits
    yl <- range(unlist(fit_Qrn))
    # set up empty plot
    plot(
      0,
      type = "n", xlim = xl, ylim = yl,
      xaxt = "n", yaxt = "n", bty = "n",
      xlab = expression(g[n](W)),
      ylab = expression("E[Y-" * Q[n](W) * " | " * g[n](W) * "]"), ...
    )
    # add axes
    axis(side = 1)
    axis(side = 2)
    # add lines
    invisible(lapply(fit_Qrn, lines, x = predP, lwd = 2, col = "gray50"))
  }
  # ------------------
  # plot grn fit
  # ------------------
  if (!(class(x$grnMod) == "NULL")) {
    # only plot if univariate reduction
    reduction <- as.list(x$call)$reduction
    if (is.null(reduction)) reduction <- "univariate"
    if (reduction == "univariate") {
      # xlim = range of gn
      xl <- range(Qn)
      # prediction points
      predP <- seq(xl[1], xl[2], length = nPoints)
      ## get fitted values of g_{n,r,1}
      fit_grn1 <- lapply(x$grnMod, function(y) {
        newDat <- data.frame(Qn = predP)
        if ("glm" %in% class(y[[listInd]]$fm1)) {
          predict(y[[listInd]]$fm1, newdata = newDat, type = "response")
        } else if (class(y[[listInd]]$fm1) == "SuperLearner") {
          pred <- predict(y[[listInd]]$fm1, newdata = newDat)
          # get sl prediction if meta learning did not fail
          if (!all(y[[listInd]]$fm1$coef == 0)) {
            pred$pred
            # otherwise get discrete super learner
          } else {
            pred$library.predict[, which.min(y$cvRisk)]
          }
        } else {
          predict(y[[listInd]]$fm1$fit, newdata = newDat, type = "response")
        }
      })
      # get ylimits
      yl <- range(unlist(fit_grn1))
      # set up empty plot
      plot(
        0,
        type = "n", xlim = xl, ylim = yl,
        xaxt = "n", yaxt = "n", bty = "n",
        xlab = expression(Q[n](W)),
        ylab = expression("E[{" * A - g[n](W) * "} / " * g[n](W) * " | " *
          Q[n](W) * "]", ...)
      )
      # add axes
      axis(side = 1)
      axis(side = 2)
      # add lines
      invisible(lapply(fit_grn1, lines, x = predP, lwd = 2, col = "gray50"))

      ## get fitted values of g_{n,r,2}
      fit_grn2 <- lapply(x$grnMod, function(y) {
        newDat <- data.frame(Qn = predP)
        if ("glm" %in% class(y[[listInd]]$fm2)) {
          predict(y[[listInd]]$fm2, newdata = newDat, type = "response")
        } else if (class(y[[listInd]]$fm2) == "SuperLearner") {
          pred <- predict(y[[listInd]]$fm2, newdata = newDat)
          # get sl prediction if meta learning did not fail
          if (!all(y[[listInd]]$fm2$coef == 0)) {
            pred$pred
            # otherwise get discrete super learner
          } else {
            pred$library.predict[, which.min(y$cvRisk)]
          }
        } else {
          predict(y[[listInd]]$fm2$fit, newdata = newDat, type = "response")
        }
      })
      # get ylimits
      yl <- range(unlist(fit_grn2))
      # set up empty plot
      plot(
        0,
        type = "n", xlim = xl, ylim = yl,
        xaxt = "n", yaxt = "n", bty = "n",
        xlab = expression(Q[n](W)),
        ylab = expression("E[A | " * Q[n](W) * "]"),
        ...
      )
      # add axes
      axis(side = 1)
      axis(side = 2)
      # add lines
      invisible(lapply(fit_grn2, lines, x = predP, lwd = 2, col = "gray50"))
    }
  }
}

#' Helper function to reorder lists according to cvFolds
#'
#' @param a_list Structured list of nuisance parameters
#' @param a_0 Treatment levels
#' @param validRows List of rows of data in validation data for
#' each split.
#' @param grn_ind Structure of grn call is slightly different
#' @param n_SL Number of super learners. If >1, then predictions
#' are averaged
#' @param n Sample size
#' @param for_se_cv Is this being used to average over 
#' cross-validated standard errors? Affects index of \code{a_list}. 
reorder_list <- function(a_list,
                         a_0,
                         validRows,
                         n_SL = 1,
                         grn_ind = FALSE,
                         n,
                         for_se_cv = FALSE) {
  n_cvFolds <- length(validRows) / n_SL
  est_index <- ifelse(for_se_cv, 3, 1)

  reduced_outList <- vector(mode = "list", length = length(a_0))

  for (i in seq_along(reduced_outList)) {
    if (!grn_ind) {
      reduced_outList[[i]] <- rep(0, n)
    } else {
      reduced_outList[[i]] <- data.frame(grn1 = rep(0, n), grn2 = rep(0, n))
    }
  }

  # re-order predictions
  for (v in seq_len(n_SL)) {
    outListValid <- unlist(a_list[(n_cvFolds * (v - 1) + 1):(v * n_cvFolds)],
      recursive = FALSE, use.names = FALSE
    )
    # this is in 0/1 format
    outListUnOrd <- do.call(Map, c(c, outListValid[seq(est_index, length(outListValid), length(a_list[[1]]))]))
    outList <- vector(mode = "list", length = length(a_0))
    if (!grn_ind) {
      for (i in seq_along(a_0)) {
        outList[[i]] <- rep(NA, n)
        # works because validRows are the same across repeated SLs
        outList[[i]][unlist(validRows)[1:n]] <- outListUnOrd[[i]]
      }
    } else {
      for (i in seq_along(a_0)) {
        outList[[i]] <- data.frame(grn1 = rep(NA, n), grn2 = rep(NA, n))
        outList[[i]][unlist(validRows), "grn1"] <- unlist(outListUnOrd[[i]][seq(1, 2 * n_cvFolds, by = 2)],
          use.names = FALSE
        )
        outList[[i]][unlist(validRows), "grn2"] <- unlist(outListUnOrd[[i]][seq(2, 2 * n_cvFolds, by = 2)],
          use.names = FALSE
        )
      }
    }
    reduced_outList <- mapply(x = reduced_outList, y = outList, function(x, y) {
      x + y
    }, SIMPLIFY = FALSE)
  }
  out <- lapply(reduced_outList, function(x) {
    x / n_SL
  })
  return(out)
}

#' Help function to extract models from fitted object
#' @param a_list Structured list of nuisance parameters
extract_models <- function(a_list) {
  outListValid <- unlist(a_list, recursive = FALSE, use.names = FALSE)
  outListValid[seq(2, length(outListValid), 3)]
}

#' Make list of rows in each validation fold.
#' @param cvFolds Numeric number of cv folds
#' @param n Number of observations
#' @param ... Other arguments
make_validRows <- function(cvFolds, n, ...) {
  if (length(cvFolds) > 1) {
    stopifnot(length(cvFolds) == n)
    # comes in as vector of fold assignments
    # split up into a list of id's
    validRows <- sapply(sort(unique(cvFolds)), function(f) {
      which(cvFolds == f)
    }, simplify = FALSE)
  } else if (cvFolds != 1) {
    # split data up
    validRows <- split(sample(seq_len(n)), rep(seq_len(cvFolds), length = n))
  } else {
    # no cross-validation
    validRows <- list(seq_len(n))
  }
  return(validRows)
}

#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function() {
  computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights,
                          errorsInLibrary = NULL, ...) {
    cvRisk <- apply(Z, 2, function(x) {
      mean(obsWeights * (x -
        Y)^2)
    })
    names(cvRisk) <- libraryNames
    compute <- function(x, y, wt = rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      fit <- tryCatch(
        {
          quadprog::solve.QP(
            Dmat = D, dvec = d, Amat = A,
            bvec = bvec, meq = 1
          )
        },
        error = function(e) {
          out <- list()
          class(out) <- "error"
          out
        }
      )
      invisible(fit)
    }
    modZ <- Z
    naCols <- which(apply(Z, 2, function(z) {
      all(z == 0)
    }))
    anyNACols <- length(naCols) > 0
    if (anyNACols) {
      warning(paste0(
        paste0(libraryNames[naCols], collapse = ", "),
        " have NAs.", "Removing from super learner."
      ))
    }
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    # if (anyDupCols) {
    #   warning(paste0(
    #     paste0(libraryNames[dupCols], collapse = ", "),
    #     " are duplicates of previous learners.", " Removing from super learner."
    #   ))
    # }
    if (anyDupCols | anyNACols) {
      rmCols <- unique(c(naCols, dupCols))
      modZ <- Z[, -rmCols, drop = FALSE]
    }
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    if (class(fit) != "error") {
      coef <- fit$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (class(fit) != "error") {
      if (anyDupCols | anyNACols) {
        ind <- c(seq_along(coef), rmCols - 0.5)
        coef <- c(coef, rep(0, length(rmCols)))
        coef <- coef[order(ind)]
      }
      coef[coef < 1e-04] <- 0
      coef <- coef / sum(coef)
    }
    if (!sum(coef) > 0) {
      warning("All algorithms have zero weight", call. = FALSE)
    }
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  computePred <- function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(
    require = "quadprog", computeCoef = computeCoef,
    computePred = computePred
  )
  invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
tmp_method.CC_nloglik <- function() {
  computePred <- function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      stop("All metalearner coefficients are zero, cannot compute prediction.")
    }
    stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
      matrix(coef[coef != 0]))
  }
  computeCoef <- function(Z, Y, libraryNames, obsWeights, control,
                          verbose, ...) {
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    modZ <- Z
    if (anyDupCols) {
      # warning(paste0(
      #   paste0(libraryNames[dupCols], collapse = ", "),
      #   " are duplicates of previous learners.", " Removing from super learner."
      # ))
      modZ <- modZ[, -dupCols, drop = FALSE]
    }
    modlogitZ <- trimLogit(modZ, control$trimLogit)
    logitZ <- trimLogit(Z, control$trimLogit)
    cvRisk <- apply(logitZ, 2, function(x) {
      -sum(2 * obsWeights *
        ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x,
          log.p = TRUE,
          lower.tail = FALSE
        )))
    })
    names(cvRisk) <- libraryNames
    obj_and_grad <- function(y, x, w = NULL) {
      y <- y
      x <- x
      function(beta) {
        xB <- x %*% cbind(beta)
        loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
          y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w)) {
          loglik <- loglik * w
        }
        obj <- -2 * sum(loglik)
        p <- stats::plogis(xB)
        grad <- if (is.null(w)) {
          2 * crossprod(x, cbind(p - y))
        } else {
          2 * crossprod(x, w * cbind(p - y))
        }
        list(objective = obj, gradient = grad)
      }
    }
    lower_bounds <- rep(0, ncol(modZ))
    upper_bounds <- rep(1, ncol(modZ))
    if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] <- 0
    }
    r <- tryCatch(
      {
        nloptr::nloptr(
          x0 = rep(1 / ncol(modZ), ncol(modZ)),
          eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds,
          ub = upper_bounds, eval_g_eq = function(beta) {
            (sum(beta) -
              1)
          }, eval_jac_g_eq = function(beta) rep(1, length(beta)),
          opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08)
        )
      },
      error = function(e) {
        out <- list()
        class(out) <- "error"
        out
      }
    )
    if (r$status < 1 || r$status > 4) {
      warning(r$message)
    }
    if (class(r) != "error") {
      coef <- r$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (anyDupCols) {
      ind <- c(seq_along(coef), dupCols - 0.5)
      coef <- c(coef, rep(0, length(dupCols)))
      coef <- coef[order(ind)]
    }
    coef[coef < 1e-04] <- 0
    coef <- coef / sum(coef)
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}

#' Helper function for averaging lists of estimates
#' generated in the main \code{for} loop of \code{drtmle}
#' 
#' @param est_cov_list A list with named entries \code{est} and \code{cov}
average_est_cov_list <- function(est_cov_list){
  length_list <- length(est_cov_list)
  all_ests <- lapply(est_cov_list, "[[", "est")
  all_cov <- lapply(est_cov_list, "[[", "cov")
  avg_est <- Reduce("+", all_ests) / length_list
  avg_cov <- Reduce("+", all_cov) / length_list
  return(list(est = avg_est, cov = avg_cov))
}

#' Helper function to average convergence results and drtmle
#' influence function estimates over multiple fits
#' @param ic_list List of influence function estimates
average_ic_list <- function(ic_list){
  length_list <- length(ic_list)
  all_mean_eif <- unlist(lapply(ic_list, "[[", "mean_eif"), use.names = FALSE)
  all_mean_missQ <- unlist(lapply(ic_list, "[[", "mean_missQ"), use.names = FALSE)
  all_mean_missg <- unlist(lapply(ic_list, "[[", "mean_missg"), use.names = FALSE)
  all_ic <- lapply(ic_list, "[[", "ic")
  avg_ic <- Reduce("+", all_ic) / length_list
  return(list(mean_eif = all_mean_eif,
              mean_missQ = all_mean_missQ,
              mean_missg = all_mean_missg,
              ic  = avg_ic))
}

#' Helper function to properly format partially cross-validated predictions 
#' from a fitted super learner.
#' 
#' @param fit_sl A fitted \code{SuperLearner} object with 
#' \code{control$saveCVFitLibrary = TRUE}
#' @param a_0 Treatment level to set. If \code{NULL}, assume this function
#' is being used to get partially cross-validated propensity score predictions.
#' @param W A \code{data.frame} of named covariates.
#' @param include A boolean vector indicating which observations were actually
#' used to fit the regression. 
#' @param easy A boolean indicating whether the predictions can be 
#' computed the "easy" way, i.e., based just on the Z matrix from SuperLearner.
#' This is possible for propensity score models when no missing data AND no 
#' stratification. 
partial_cv_preds <- function(fit_sl, a_0, W = NULL,
                             include = NULL, easy = FALSE){
  n_algo <- length(fit_sl$cvRisk)
  n_folds <- length(fit_sl$validRows)

  if(!easy){ 
    n <- length(W[,1])
  }else{ # if used in easy scenario, then fit_sl will have been 
         # fit using all observations and no W will enter, so we'll check
         # the Z matrix that is used to generate predictions below
    n <- length(fit_sl$Z[,1])
  }
  alpha_hat <- matrix(fit_sl$coef, nrow = n_algo)
  if(!easy){
    # rslt_list will eventually hold cross-validated predictions for 
    # all observations that were actually used to fit the regression
    rslt_list <- vector(mode = "list", length = n_folds)
    for(v in seq_len(n_folds)){
      # who's in the validation fold of the included folks
      foldv_ids <- fit_sl$validRows[[v]]
      # these are the people who were not used to fit these models
      foldv_models <- fit_sl$cvFitLibrary[[v]]
      rslt <- matrix(NA, nrow = length(foldv_ids), ncol = n_algo)
      for(k in seq_len(n_algo)){
        # predict under a_0
        if(!is.null(a_0)){
          rslt[ , k] <- predict(foldv_models[[k]], newdata = data.frame(A = a_0, W)[include,][foldv_ids,])
        }else{
          rslt[ , k] <- predict(foldv_models[[k]], newdata = W[include,][foldv_ids,])
        }
      }
      rslt_list[[v]] <- rslt
    }
    # combine using weights from full super learner
    pred_list <- vector(mode = "list", length = n_folds)
    for(v in seq_len(n_folds)){
      pred_list[[v]] <- rslt_list[[v]] %*% alpha_hat
    }
    reorder_preds <- rep(NA, n)
    # fill in observations in regression with cross-validated prediction
    reorder_preds[include][unlist(fit_sl$validRows)] <- unlist(pred_list, use.names = FALSE)
    # all others fill in with prediction from super learner
    if(any(!include)){
      if(!is.null(a_0)){
        reorder_preds[!include] <- predict(fit_sl, newdata = data.frame(A = a_0, W)[!include,])[[1]]
      }else{
        reorder_preds[!include] <- predict(fit_sl, newdata = W[!include,])[[1]]
      }
    }
  }else{ # when a_0 is NULL
    # in this case, we're operating on a propensity score model 
    # in which case, fit_sl$Z already has cross-validated predictions
    reorder_preds <- as.numeric(fit_sl$Z %*% alpha_hat)
  }
  return(reorder_preds)
}
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
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # simulate data
#' set.seed(123456)
#' n <- 100
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#' # fit drtmle with maxIter = 1 to run fast
#' fit1 <- drtmle(W = W, A = A, Y = Y, a_0 = c(1,0),
#'                family=binomial(),
#'                stratify=FALSE,
#'                SL_Q=c("SL.glm","SL.mean","SL.glm.interaction"),
#'                SL_g=c("SL.glm","SL.mean","SL.glm.interaction"),
#'                SL_Qr="SL.npreg", SL_gr="SL.npreg",
#'                maxIter = 1, returnModels = TRUE)
#' # plot the reduced-dimension regression fits (not run)
#' \dontrun{plot(fit1)}
#
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

  # ------------------
  # plot grn fit
  # ------------------
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

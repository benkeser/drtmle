#' Compute confidence intervals for drtmle and adaptive_iptw@
#' @param ... Arguments to be passed to method
#' @export
ci <- function(...) {
  UseMethod("ci")
}

#' Confidence intervals for drtmle objects
#'
#' @param object An object of class \code{"drtmle"}
#' @param est A vector indicating for which estimators to return a
#'  confidence interval. Possible estimators include the TMLE with doubly robust
#'  inference (\code{"drtmle"}, recommended), the AIPTW with additional
#'  correction for misspecification (\code{"aiptw_c"}, not recommended), the
#'  standard TMLE (\code{"tmle"}, recommended only for comparison to "drtmle"),
#'  the standard AIPTW (\code{"aiptw"}, recommended only for comparison to
#'  "drtmle"), and G-computation (\code{"gcomp"}, not recommended).
#' @param level The nominal coverage probability of the desired confidence
#'  interval (should be between 0 and 1). Default computes 95\% confidence
#'  intervals.
#' @param contrast Specifies the parameter for which to return confidence
#'  intervals. If \code{contrast=NULL}, then confidence intervals for the
#'  marginal means are computed. If instead, \code{contrast} is a numeric vector
#'  of ones, negative ones, and zeros to define linear combinations of the
#'  various means (e.g., to estimate an average treatment effect, see example).
#'  Finally, \code{contrast} can be a list with named functions \code{f},
#'  \code{f_inv}, \code{h}, and \code{fh_grad}. The first two functions should
#'  take as input argument \code{eff}. Respectively, these specify which
#'  transformation of the effect measure to compute the confidence interval for
#'  and the inverse transformation to put the confidence interval back on the
#'  original scale. The function \code{h} defines the contrast to be estimated
#'  and should take as input \code{est}, a vector of the same length as
#'  \code{object$a_0}, and output the desired contrast. The function
#'  \code{fh_grad} is the gradient of the function \code{h}. See examples and
#'  vignette for more information.
#' @param ... Other options (not currently used).
#' @importFrom stats qnorm
#' @export
#' @method ci drtmle
#' @return An object of class \code{"ci.drtmle"} with point estimates and
#'  confidence intervals of the specified level.
#'
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # simulate data
#' set.seed(123456)
#' n <- 100
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n, 1, plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
#' \donttest{
#' # fit drtmle with maxIter = 1 to run fast
#' fit1 <- drtmle(
#'   W = W, A = A, Y = Y, a_0 = c(1, 0),
#'   family = binomial(),
#'   stratify = FALSE,
#'   SL_Q = c("SL.glm", "SL.mean"),
#'   SL_g = c("SL.glm", "SL.mean"),
#'   SL_Qr = "SL.npreg",
#'   SL_gr = "SL.npreg", maxIter = 1
#' )
#'
#' # get confidence intervals for each mean
#' ci_mean <- ci(fit1)
#'
#' # get confidence intervals for ATE
#' ci_ATE <- ci(fit1, contrast = c(1, -1))
#'
#' # get confidence intervals for risk ratio by
#' # computing CI on log scale and back-transforming
#' myContrast <- list(
#'   f = function(eff) {
#'     log(eff)
#'   },
#'   f_inv = function(eff) {
#'     exp(eff)
#'   },
#'   h = function(est) {
#'     est[1] / est[2]
#'   },
#'   fh_grad = function(est) {
#'     c(1 / est[1], -1 / est[2])
#'   }
#' )
#' ci_RR <- ci(fit1, contrast = myContrast)
#' }
ci.drtmle <- function(object, est = c("drtmle"), level = 0.95,
                      contrast = NULL, ...) {
  if (class(object) != "drtmle") {
    stop("ci only works with drtmle objects")
  }
  out <- vector(mode = "list", length = length(est))
  names(out) <- est
  # if no contrast then return an CI for each
  # covariate-adjusted mean
  if (is.null(contrast)) {
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 3)
      for (j in seq_along(object$a_0)) {
        out[[i]][j, ] <-
          rep(object[[est[i]]]$est[j], 3) +
          stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
            rep(sqrt(object[[est[i]]]$cov[j, j]), 3)
      }
      row.names(out[[i]]) <- object$a_0
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  } else if (is.numeric(contrast)) {
    # check that contrast came in correctly
    if (length(contrast) != length(object$a_0)) {
      stop("length of contrast vector not equal to length of a_0")
    }
    if (!all(contrast %in% c(-1, 1, 0))) {
      stop("contrast should only be -1, 1, or 0")
    }
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
      g <- matrix(contrast, nrow = length(contrast))
      p <- matrix(object[[est[i]]]$est, nrow = length(contrast))
      v <- object[[est[i]]]$cov
      thisC <- t(g) %*% p
      thisSe <- sqrt(t(g) %*% v %*% g)
      out[[i]][1, ] <- rep(thisC, 3) +
        stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
          rep(thisSe, 3)
      # E[Y(a_0[1])] - E[Y(a_0[2])]
      indMinus <- which(contrast == -1)
      indPlus <- which(contrast == 1)
      plusTerms <- minusTerms <- ""
      if (length(indMinus) > 0) {
        minusTerms <- paste0("E[Y(", object$a_0[indMinus], ")]")
      }
      if (length(indPlus) > 0) {
        plusTerms <- paste0("E[Y(", object$a_0[indPlus], ")]")
      }
      thisName <- paste0(
        paste0(plusTerms, collapse = "+"),
        ifelse(length(indMinus) > 0, "-", ""),
        paste0(minusTerms, collapse = "-")
      )

      row.names(out[[i]]) <- thisName
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  } else if (is.list(contrast)) {
    if (!all(c("f", "f_inv", "h", "fh_grad") %in% names(contrast))) {
      stop("some function missing in contrast. see ?ci.drtmle for help.")
    }
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
      thisC <- do.call(contrast$h, args = list(est = object[[est[i]]]$est))
      f_thisC <- do.call(contrast$f, args = list(eff = thisC))

      grad <- matrix(do.call(
        contrast$fh_grad,
        args = list(est = object[[est[i]]]$est)
      ), nrow = length(object$a_0))
      v <- object[[est[i]]]$cov
      thisSe <- sqrt(t(grad) %*% v %*% grad)
      transformCI <- rep(f_thisC, 3) +
        stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
          rep(thisSe, 3)
      out[[i]][1, ] <- do.call(contrast$f_inv, args = list(eff = transformCI))
      row.names(out[[i]]) <- c("user contrast")
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  }
  class(out) <- "ci.drtmle"
  return(out)
}

#' Confidence intervals for adaptive_iptw objects
#'
#' Estimate confidence intervals for objects of class \code{"adaptive_iptw"}
#'
#' @param object An object of class \code{"adaptive_iptw"}
#' @param est A vector indicating for which estimators to return a
#' confidence interval. Possible estimators include the TMLE IPTW
#' (\code{"iptw_tmle"}, recommended), the one-step IPTW
#' (\code{"iptw_os"}, not recommended), the standard IPTW
#' (\code{"iptw"}, recommended only for comparison to the other two estimators).
#' @param level The nominal coverage probability of the desired confidence
#'  interval (should be between 0 and 1). Default computes 95\% confidence
#'  intervals.
#' @param contrast Specifies the parameter for which to return confidence
#'  intervals. If \code{contrast=NULL}, then confidence intervals for the
#'  marginal means are computed. If instead, \code{contrast} is a numeric vector
#'  of ones, negative ones, and zeros to define linear combinations of the
#'  various means (e.g., to estimate an average treatment effect, see example).
#'  Finally, \code{contrast} can be a list with named functions \code{f},
#'  \code{f_inv}, \code{h}, and \code{fh_grad}. The first two functions should
#'  take as input argument \code{eff}. Respectively, these specify which
#'  transformation of the effect measure to compute the confidence interval for
#'  and the inverse transformation to put the confidence interval back on the
#'  original scale. The function \code{h} defines the contrast to be estimated
#'  and should take as input \code{est}, a vector of the same length as
#'  \code{object$a_0}, and output the desired contrast. The function
#'  \code{fh_grad} is the gradient of the function \code{h}. See examples and
#'  vignette for more information.
#' @param ... Other options (not currently used).
#' @export
#' @method ci adaptive_iptw
#' @return An object of class \code{"ci.adaptive_iptw"} with point estimates and
#'  confidence intervals of the specified level.
#'
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # fit adaptive_iptw
#' set.seed(123456)
#' n <- 200
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n, 1, plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
#'
#' fit1 <- adaptive_iptw(
#'   W = W, A = A, Y = Y, a_0 = c(1, 0),
#'   SL_g = c("SL.glm", "SL.mean", "SL.step"),
#'   SL_Qr = "SL.glm"
#' )
#'
#' # get confidence intervals for each mean
#' ci_mean <- ci(fit1)
#'
#' # get confidence intervals for ATE
#' ci_ATE <- ci(fit1, contrast = c(1, -1))
#'
#' # get confidence intervals for risk ratio
#' # by inputting own contrast function
#' # this computes CI on log scale and back transforms
#' myContrast <- list(
#'   f = function(eff) {
#'     log(eff)
#'   },
#'   f_inv = function(eff) {
#'     exp(eff)
#'   },
#'   h = function(est) {
#'     est[1] / est[2]
#'   },
#'   fh_grad = function(est) {
#'     c(1 / est[1], -1 / est[2])
#'   }
#' )
#' ci_RR <- ci(fit1, contrast = myContrast)
ci.adaptive_iptw <- function(object, est = c("iptw_tmle"), level = 0.95,
                             contrast = NULL, ...) {
  if (any(est == "iptw")) {
    stop("Theory does not support inference for naive IPTW with super learner.")
  }
  if (class(object) != "adaptive_iptw") {
    stop("ci only works with adaptive_iptw objects")
  }
  out <- vector(mode = "list", length = length(est))
  names(out) <- est
  # if no contrast then return an CI for each
  # covariate-adjusted mean
  if (is.null(contrast)) {
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 3)
      for (j in seq_along(object$a_0)) {
        out[[i]][j, ] <-
          rep(object[[est[i]]]$est[j], 3) +
          stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
            rep(sqrt(object[[est[i]]]$cov[j, j]), 3)
      }
      row.names(out[[i]]) <- object$a_0
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  } else if (is.numeric(contrast)) {
    # check that contrast came in correctly
    if (length(contrast) != length(object$a_0)) {
      stop("length of contrast vector not equal to length of a_0")
    }
    if (!all(contrast %in% c(-1, 1, 0))) {
      stop("contrast should only be -1, 1, or 0")
    }
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
      g <- matrix(contrast, nrow = length(contrast))
      p <- matrix(object[[est[i]]]$est, nrow = length(contrast))
      v <- object[[est[i]]]$cov
      thisC <- t(g) %*% p
      thisSe <- sqrt(t(g) %*% v %*% g)
      out[[i]][1, ] <- rep(thisC, 3) +
        stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) * rep(thisSe, 3)
      indMinus <- which(contrast == -1)
      indPlus <- which(contrast == 1)
      plusTerms <- minusTerms <- ""
      if (length(indMinus) > 0) {
        minusTerms <- paste0("E[Y(", object$a_0[indMinus], ")]")
      }
      if (length(indPlus) > 0) {
        plusTerms <- paste0("E[Y(", object$a_0[indPlus], ")]")
      }
      thisName <- paste0(
        paste0(plusTerms, collapse = "+"),
        ifelse(length(indMinus) > 0, "-", ""),
        paste0(minusTerms, collapse = "-")
      )

      row.names(out[[i]]) <- thisName
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  } else if (is.list(contrast)) {
    if (!all(c("f", "f_inv", "h", "fh_grad") %in% names(contrast))) {
      stop("some function missing in contrast.")
    }
    for (i in seq_along(est)) {
      out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
      thisC <- do.call(contrast$h, args = list(est = object[[est[i]]]$est))
      f_thisC <- do.call(contrast$f, args = list(eff = thisC))

      grad <- matrix(do.call(
        contrast$fh_grad,
        args = list(est = object[[est[i]]]$est)
      ), nrow = length(object$a_0))
      v <- object[[est[i]]]$cov
      thisSe <- sqrt(t(grad) %*% v %*% grad)
      transformCI <- rep(f_thisC, 3) +
        stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
          rep(thisSe, 3)
      out[[i]][1, ] <- do.call(contrast$f_inv, args = list(eff = transformCI))
      row.names(out[[i]]) <- c("user contrast")
      colnames(out[[i]]) <- c("est", "cil", "ciu")
    }
  }
  class(out) <- "ci.adaptive_iptw"
  return(out)
}

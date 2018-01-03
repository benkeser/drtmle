#' fluctuateQ1
#'
#' Function called internally by drtmle to perform the first fluctuation
#' of the initial estimator of Q (i.e., solves the original EIF estimating eqn)
#'
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param a_0 A list of fixed treatment values
#'
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#'

fluctuateQ1 <- function(Y, A, W, DeltaA, DeltaY, Qn, gn, a_0) {
  QnStar <- mapply(a = a_0, Q = Qn, g = gn, FUN = function(x, a, Q, g) {
    l <- min(Y, na.rm = TRUE)
    u <- max(Y, na.rm = TRUE)
    Yscale <- (Y - l) / (u - l)
    off <- SuperLearner::trimLogit((Q - l) / (u - l))
    H1 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / g
    suppressWarnings(
      fm <- stats::glm(
        Yscale ~ -1 + offset(off) + H1, start = 0,
        data = data.frame(Y = Y, off = off, H1 = H1), family = "binomial"
      )
    )
    Qnstar <- stats::predict(
      fm, type = "response",
      newdata = data.frame(off = off, H1 = 1 / g)
    ) * (u - l) + l
    list(est = Qnstar, eps = fm$coef)
  }, SIMPLIFY = FALSE)
  QnStar
}

#' fluctuateG
#'
#' Function called internally by drtmle to perform the fluctuation
#' of the initial estimator of g (i.e., solves the new estimating eqn that results
#' from misspecification of Q)
#'
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param Qrn A list of reduced-dimension regression estimates evaluated on observed data
#' @param coefTol A tolerance level on the magnitude of the coefficient that flags the
#' result as potentially the result of numeric instability.
#' @param tolg The lower bound on propensity score estimates
#' @param a_0 A list of fixed treatment values
#'
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#'

fluctuateG <- function(Y, A, W, DeltaY, DeltaA, a_0, gn, Qrn, tolg, coefTol=1e5) {
  gnStar <- mapply(a = a_0, g = gn, Qr = Qrn, FUN = function(x, a, g, Qr) {
    H1 <- Qr / g
    off <- SuperLearner::trimLogit(g, tolg)
    thisA <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1)
    suppressWarnings(
      fm <- stats::glm(
        thisA ~ -1 + offset(off) + H1, start = 0,
        data = data.frame(thisA = thisA, off = off, H1 = H1),
        family = "binomial"
      )
    )
    if (!fm$converged | abs(fm$coef) > coefTol) {
      suppressWarnings(
        fm <- stats::glm(
          thisA ~ -1 + offset(off) + H1,
          data = data.frame(thisA = thisA, off = off, H1 = H1),
          family = "binomial"
        )
      )
      if (!fm$converged | abs(fm$coef) > coefTol) {
        # warning("No sane fluctuation found for G this iteration. Check mean of IC.")
      }
    }
    pred <- stats::predict(fm, type = "response")
    pred[pred < tolg] <- tolg
    list(est = pred, eps = fm$coef)
  }, SIMPLIFY = FALSE)
  gnStar
}


#' fluctuateQ2
#'
#' Function called internally by drtmle to perform the second fluctuation
#' of the initial estimator of Q (i.e., solves the new estimating eqn that results
#' from misspecification of g)
#'
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param grn A list of reduced-dimension regression estimates evaluated on observed data
#' @param coefTol A tolerance level on the magnitude of the coefficient that flags the
#' result as potentially the result of numeric instability.
#' @param reduction A character indicating what reduced dimension regression was used.
#' @param a_0 A list of fixed treatment values
#'
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#'


fluctuateQ2 <- function(Y, A, W, DeltaY, DeltaA,
                        Qn, gn, grn, a_0, reduction, coefTol=1e5) {
  QnStar <- mapply(a = a_0, Q = Qn, g = gn, gr = grn, FUN = function(a, Q, g, gr) {
    l <- min(Y, na.rm = TRUE)
    u <- max(Y, na.rm = TRUE)
    Yscale <- (Y - l) / (u - l)
    off <- SuperLearner::trimLogit((Q - l) / (u - l))
    if (reduction == "univariate") {
      H2 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / gr$grn2 * gr$grn1
    } else if (reduction == "bivariate") {
      H2 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / gr$grn2 *
        (gr$grn2 - g) / g
    }
    suppressWarnings(
      fm <- stats::glm(
        Yscale ~ -1 + offset(off) + H2, start = 0,
        data = data.frame(Y = Y, off = off, H2 = H2), family = "binomial"
      )
    )
    if (!fm$converged | abs(max(fm$coef)) > coefTol) {
      # if it doesn't converge, try with no starting values
      suppressWarnings(
        fm <- stats::glm(
          Yscale ~ -1 + offset(off) + H2,
          data = data.frame(Y = Y, off = off, H2 = H2),
          family = "binomial"
        )
      )
      if (!fm$converged | abs(max(fm$coef)) > coefTol) {
        # warning("No sane fluctuation found. Proceeding using current estimates.")
        if (reduction == "univariate") {
          return(list(est = Q, eps = rep(0, 2)))
        } else if (reduction == "bivariate") {
          return(list(est = Q, eps = rep(0, 2)))
        }
      }
    }

    if (reduction == "univariate") {
      Qnstar <- stats::predict(
        fm, type = "response", newdata = data.frame(
          off = off, H2 = 1 / gr$grn2 * gr$grn1
        )
      ) * (u - l) + l

      return(list(est = Qnstar, eps = fm$coef))
    } else if (reduction == "bivariate") {
      Qnstar <- stats::predict(
        fm, type = "response", newdata = data.frame(
          off = off,
          H2 = 1 / gr$grn2 * (gr$grn2 - g) / g
        )
      ) * (u - l) + l
      return(list(est = Qnstar, eps = fm$coef))
    }
  }, SIMPLIFY = FALSE)
  QnStar
}

#' fluctuateQ
#'
#' Function called internally by drtmle to perform simultaneous fluctuation
#' of the initial estimator of Q (i.e., solves both EIF estimating eqn and
#' the new estimating eqn that results from misspecification of g)
#'
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param grn A list of reduced-dimension regression estimates evaluated on observed data
#' @param coefTol A tolerance level on the magnitude of the coefficient that flags the
#' result as potentially the result of numeric instability.
#' @param reduction A character indicating what reduced dimension regression was used.
#' @param a_0 A list of fixed treatment values
#'
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#'
fluctuateQ <- function(Y, A, W, DeltaY, DeltaA,
                       Qn, gn, grn, a_0, reduction, coefTol=1e5) {
  QnStar <- mapply(a = a_0, Q = Qn, g = gn, gr = grn, FUN = function(a, Q, g, gr) {
    l <- min(Y, na.rm = TRUE)
    u <- max(Y, na.rm = TRUE)
    Yscale <- (Y - l) / (u - l)
    off <- SuperLearner::trimLogit((Q - l) / (u - l))

    H1 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / g
    if (reduction == "univariate") {
      H2 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / gr$grn2 * gr$grn1
    } else if (reduction == "bivariate") {
      H2 <- as.numeric(A == a & DeltaA == 1 & DeltaY == 1) / gr$grn2 * (gr$grn2 - g) / g
    }

    suppressWarnings(
      fm <- stats::glm(
        Yscale ~ -1 + offset(off) + H1 + H2, start = c(0, 0),
        data = data.frame(Y = Y, off = off, H1 = H1, H2 = H2),
        family = "binomial"
      )
    )
    if (!fm$converged | abs(max(fm$coef)) > coefTol) {
      # if it doesn't converge, try with no starting values
      suppressWarnings(
        fm <- stats::glm(
          Yscale ~ -1 + offset(off) + H1 + H2,
          data = data.frame(Y = Y, off = off, H1 = H1, H2 = H2),
          family = "binomial"
        )
      )
      if (!fm$converged | abs(max(fm$coef)) > coefTol) {
        # warning("No sane fluctuation found. Proceeding using current estimates.")
        if (reduction == "univariate") {
          return(list(est = Q, eps = rep(0, 2)))
        } else if (reduction == "bivariate") {
          return(list(est = Q, eps = rep(0, 2)))
        }
      }
    }

    if (reduction == "univariate") {
      Qnstar <- stats::predict(
        fm, type = "response", newdata = data.frame(
          off = off, H1 = 1 / g,
          H2 = 1 / gr$grn2 * gr$grn1
        )
      ) * (u - l) + l
      return(list(est = Qnstar, eps = fm$coef))
    } else if (reduction == "bivariate") {
      Qnstar <- stats::predict(
        fm, type = "response", newdata = data.frame(
          off = off, H1 = 1 / g,
          H2 = 1 / gr$grn2 * (gr$grn2 - g) / g
        )
      ) * (u - l) + l
      return(list(est = Qnstar, eps = fm$coef))
    }
  }, SIMPLIFY = FALSE)
  QnStar
}

#' Compute asymptotically linear IPTW estimators with super learning
#' for the propensity score
#'
#' @param W A \code{data.frame} of named covariates
#' @param A A \code{numeric} vector of binary treatment assignment (assumed to
#'  be equal to 0 or 1)
#' @param Y A \code{numeric} numeric of continuous or binary outcomes.
#' @param DeltaY A \code{numeric} indicator of missing outcome (assumed to be
#'  equal to 0 if missing 1 if observed)
#' @param DeltaA A \code{numeric} indicator of missing treatment (assumed to be
#'  equal to 0 if missing 1 if observed)
#' @param a_0 A vector of \code{numeric} treatment values at which to return
#'  marginal mean estimates.
#' @param stratify A \code{logical} indicating whether to estimate the missing
#'  outcome regression separately for observations with different levels of
#'  \code{A} (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param family A \code{family} object equal to either \code{binomial()} or
#'  \code{gaussian()}, to be passed to the \code{SuperLearner} or \code{glm}
#'  function.
#' @param SL_g A vector of characters describing the super learner library to be
#'  used for each of the propensity score regressions (\code{DeltaA}, \code{A},
#'  and \code{DeltaY}). To use the same library for each of the regressions (or
#'  if there is no missing data in \code{A} nor \code{Y}), a single library may
#'  be input. See \code{link{SuperLearner::SuperLearner}} for details on how
#'  super learner libraries can be specified.
#' @param SL_Qr A vector of characters or a list describing the Super Learner
#'  library to be used for the reduced-dimension outcome regression.
#' @param glm_g A list of characters describing the formulas to be used
#'  for each of the propensity score regressions (\code{DeltaA}, \code{A}, and
#'  \code{DeltaY}). To use the same formula for each of the regressions (or if
#'  there is no missing data in \code{A} nor \code{Y}), a single character
#'  formula may be input.
#' @param glm_Qr A character describing a formula to be used in the call to
#'  \code{glm} for reduced-dimension outcome regression. Ignored if
#'  \code{SL_Qr!=NULL}. The formula should use the variable name \code{'gn'}.
#' @param maxIter A numeric that sets the maximum number of iterations the TMLE
#'  can perform in its fluctuation step.
#' @param tolIC A numeric that defines the stopping criteria based on the
#'  empirical mean of the influence function.
#' @param tolg A numeric indicating the minimum value for estimates of the
#'  propensity score.
#' @param verbose A logical indicating whether to print status updates.
#' @param returnModels A logical indicating whether to return model fits for the
#'  propensity score and reduced-dimension regressions.
#' @param cvFolds A numeric equal to the number of folds to be used in
#'  cross-validated fitting of nuisance parameters. If \code{cvFolds = 1}, no
#'  cross-validation is used.
#' @param parallel A logical indicating whether to use parallelization based on
#'  \code{future} to estimate nuisance parameters in parallel. Only useful if
#'  \code{cvFolds > 1}. By default, a \code{multiprocess} evaluation scheme is
#'  invoked, using forked R processes (if supported on the OS) and background R
#'  sessions otherwise. Users may also register their own backends using the
#'  \code{future.batchtools} package.
#' @param future_hpc A character string identifying a high-performance computing
#'  backend to be used with parallelization. This should match exactly one of
#'  the options available from the \code{future.batchtools} package.
#' @param gn An optional list of propensity score estimates. If specified, the
#'  function will ignore the nuisance parameter estimation specified by
#'  \code{SL_g} and \code{glm_g}. The entries in the list should correspond to
#'  the propensity for the observed values of \code{W}, with order determined by
#'  the input to \code{a_0} (e.g., if \code{a_0 = c(0,1)} then \code{gn[[1]]}
#'  should be propensity of \code{A} = 0 and \code{gn[[2]]} should be propensity
#'  of \code{A} = 1).
#' @param ... Other options (not currently used).
#' @return An object of class \code{"adaptive_iptw"}.
#' \describe{
#'  \item{\code{iptw_tmle}}{A \code{list} of point estimates and
#'        covariance matrix for the IPTW estimator based on a targeted
#'        propensity score. }
#'  \item{\code{iptw_tmle_nuisance}}{A \code{list} of the final TMLE estimates
#'        of the propensity score (\code{$gnStar}) and reduced-dimension
#'        regression (\code{$QrnStar}) evaluated at the observed data values.}
#'  \item{\code{iptw_os}}{A \code{list} of point estimates and covariance matrix
#'        for the one-step correct IPTW estimator.}
#'  \item{\code{iptw_os_nuisance}}{A \code{list} of the initial estimates of the
#'        propensity score and reduced-dimension regression evaluated at the
#'        observed data values.}
#'  \item{\code{iptw}}{A \code{list} of point estimates for the standard IPTW
#'        estimator. No estimate of the covariance matrix is provided because
#'        theory does not support asymptotic Normality of the IPTW estimator if
#'        super learning is used to estimate the propensity score.}
#'  \item{\code{gnMod}}{The fitted object for the propensity score. Returns
#'        \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{QrnMod}}{The fitted object for the reduced-dimension regression
#'        that guards against misspecification of the outcome regression.
#'        Returns \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{a_0}}{The treatment levels that were requested for computation
#'        of covariate-adjusted means.}
#'  \item{\code{call}}{The call to \code{adaptive_iptw}.}
#' }
#'
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom doFuture registerDoFuture
#' @importFrom stats cov
#'
#' @export
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
#' # fit iptw with maxIter = 1 to run fast
#' \donttest{
#' fit1 <- adaptive_iptw(
#'   W = W, A = A, Y = Y, a_0 = c(1, 0),
#'   SL_g = c("SL.glm", "SL.mean", "SL.step"),
#'   SL_Qr = "SL.npreg", maxIter = 1
#' )
#' }
adaptive_iptw <- function(W, A, Y,
                          DeltaY = as.numeric(!is.na(Y)),
                          DeltaA = as.numeric(!is.na(A)),
                          stratify = FALSE,
                          family = if (all(Y %in% c(0, 1))) {
                            stats::binomial()
                          } else {
                            stats::gaussian()
                          },
                          a_0 = unique(A[!is.na(A)]),
                          SL_g = NULL, glm_g = NULL,
                          SL_Qr = NULL, glm_Qr = NULL,
                          returnModels = TRUE,
                          verbose = FALSE,
                          maxIter = 2,
                          tolIC = 1 / length(Y),
                          tolg = 1e-2,
                          cvFolds = 1,
                          parallel = FALSE,
                          future_hpc = NULL,
                          gn = NULL,
                          ...) {
  call <- match.call()
  # if cvFolds non-null split data into cvFolds pieces
  n <- length(Y)
  if (cvFolds != 1) {
    validRows <- split(sample(seq_len(n)), rep(1:cvFolds, length = n))
  } else {
    validRows <- list(seq_len(n))
  }
  # use futures with foreach if parallel mode
  if (!parallel) {
    future::plan(future::transparent)
  } else {
    doFuture::registerDoFuture()
    if (all(c("sequential", "uniprocess") %in% class(future::plan())) &
      is.null(future_hpc)) {
      future::plan(future::multiprocess)
    } else if (!is.null(future_hpc)) {
      if (future_hpc == "batchtools_torque") {
        future::plan(future.batchtools::batchtools_torque)
      } else if (future_hpc == "batchtools_slurm") {
        future::plan(future.batchtools::batchtools_slurm)
      } else if (future_hpc == "batchtools_sge") {
        future::plan(future.batchtools::batchtools_sge)
      } else if (future_hpc == "batchtools_lsf") {
        future::plan(future.batchtools::batchtools_lsf)
      } else if (future_hpc == "batchtools_openlava") {
        future::plan(future.batchtools::batchtools_openlava)
      } else {
        stop("The currently specified HPC backend is not (yet) available.")
      }
    }
  }
  # -------------------------------
  # estimate propensity score
  # -------------------------------
  if (is.null(gn)) {
    gnOut <- future.apply::future_lapply(
      X = validRows, FUN = estimateG,
      A = A, W = W,
      DeltaA = DeltaA, DeltaY = DeltaY,
      tolg = tolg, verbose = verbose,
      returnModels = returnModels,
      SL_g = SL_g, glm_g = glm_g,
      a_0 = a_0, stratify = stratify
    )
    # re-order predictions
    gnValid <- unlist(gnOut, recursive = FALSE, use.names = FALSE)
    gnUnOrd <- do.call(Map, c(c, gnValid[seq(1, length(gnValid), 2)]))
    gn <- vector(mode = "list", length = length(a_0))
    for (i in seq_along(a_0)) {
      gn[[i]] <- rep(NA, n)
      gn[[i]][unlist(validRows)] <- gnUnOrd[[i]]
    }
    # obtain list of propensity score fits
    gnMod <- gnValid[seq(2, length(gnValid), 2)]
  } else {
    gnMod <- NULL
  }
  # compute iptw estimator
  psi_n <- mapply(a = split(a_0, seq_along(a_0)), g = gn, function(a, g) {
    modA <- A
    modA[is.na(A)] <- -999
    modY <- Y
    modY[is.na(Y)] <- -999
    mean(as.numeric(modA == a & DeltaA == 1 & DeltaY == 1) / g * modY)
  })

  # estimate influence function
  Dno <- eval_Diptw(
    A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY, gn = gn,
    psi_n = psi_n, a_0 = a_0
  )

  # -------------------------------------
  # estimate reduced dimension Q
  # -------------------------------------
  # note that NULL is input to estimateQrn -- internally the function
  # assign Qn = 0 for all a_0 because estimateQrn estimates the regression
  # of Y - Qn on gn (which is needed for drtmle), while here we just need
  # the regression of Y on gn.
  QrnOut <- future.apply::future_lapply(
    X = validRows, FUN = estimateQrn,
    Y = Y, A = A, W = W,
    DeltaA = DeltaA, DeltaY = DeltaY,
    Qn = NULL, gn = gn, glm_Qr = glm_Qr,
    family = family, SL_Qr = SL_Qr, a_0 = a_0,
    returnModels = returnModels
  )

  # re-order predictions
  QrnValid <- unlist(QrnOut, recursive = FALSE, use.names = FALSE)
  QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1, length(QrnValid), 2)]))
  Qrn <- vector(mode = "list", length = length(a_0))
  for (i in seq_along(a_0)) {
    Qrn[[i]] <- rep(NA, n)
    Qrn[[i]][unlist(validRows)] <- QrnUnOrd[[i]]
  }
  # obtain list of propensity score fits
  QrnMod <- QrnValid[seq(2, length(QrnValid), 2)]

  Dngo <- eval_Diptw_g(
    A = A, DeltaA = DeltaA, DeltaY = DeltaY,
    Qrn = Qrn, gn = gn, a_0 = a_0
  )
  PnDgn <- lapply(Dngo, mean)

  # one-step iptw estimator
  psi.o <- mapply(
    a = psi_n, b = PnDgn, SIMPLIFY = FALSE,
    FUN = function(a, b) {
      a - b
    }
  )

  # targeted g estimator
  gnStar <- gn
  QrnStar <- Qrn
  PnDgnStar <- Inf
  ct <- 0
  # fluctuate
  while (max(abs(unlist(PnDgnStar))) > tolIC & ct < maxIter) {
    ct <- ct + 1

    # fluctuate gnStar
    gnStarOut <- fluctuateG(
      Y = Y, A = A, W = W, DeltaA = DeltaA,
      DeltaY = DeltaY, a_0 = a_0, tolg = tolg,
      gn = gnStar, Qrn = QrnStar
    )
    gnStar <- lapply(gnStarOut, function(x) {
      unlist(x$est)
    })
    # re-estimate reduced dimension regression
    QrnStarOut <- future.apply::future_lapply(
      X = validRows, FUN = estimateQrn,
      Y = Y, A = A, W = W,
      DeltaA = DeltaA, DeltaY = DeltaY,
      Qn = NULL, gn = gnStar,
      glm_Qr = glm_Qr, family = family,
      SL_Qr = SL_Qr, a_0 = a_0,
      returnModels = returnModels
    )
    # re-order predictions
    QrnValid <- unlist(QrnStarOut, recursive = FALSE, use.names = FALSE)
    QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1, length(QrnValid), 2)]))
    QrnStar <- vector(mode = "list", length = length(a_0))
    for (i in seq_along(a_0)) {
      QrnStar[[i]] <- rep(NA, n)
      QrnStar[[i]][unlist(validRows)] <- QrnUnOrd[[i]]
    }
    # obtain list of propensity score fits
    QrnMod <- QrnValid[seq(2, length(QrnValid), 2)]

    # compute influence function for fluctuated estimators
    DngoStar <- eval_Diptw_g(
      A = A, DeltaA = DeltaA, DeltaY = DeltaY,
      Qrn = QrnStar, gn = gnStar, a_0 = a_0
    )
    PnDgnStar <- future.apply::future_lapply(DngoStar, mean)
    if (verbose) {
      cat("Mean of IC       =", round(unlist(PnDgnStar), 10), "\n")
    }
  }

  # compute final tmle-iptw estimate
  # compute iptw estimator
  psi_nStar <- mapply(
    a = split(a_0, seq_along(a_0)), g = gnStar,
    function(a, g) {
      modA <- A
      modA[is.na(A)] <- -999
      modY <- Y
      modY[is.na(Y)] <- -999
      mean(as.numeric(modA == a & DeltaA == 1 &
        DeltaY == 1) / g * modY)
    }
  )

  # compute variance estimators
  # original influence function
  DnoStar <- eval_Diptw(
    A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY,
    gn = gnStar, psi_n = psi_nStar, a_0 = a_0
  )

  # covariance for tmle iptw
  DnoStarMat <- matrix(
    unlist(DnoStar) - unlist(DngoStar),
    nrow = n,
    ncol = length(a_0)
  )
  cov.t <- stats::cov(DnoStarMat) / n
  # covariate for one-step iptw
  DnoMat <- matrix(unlist(Dno) - unlist(Dngo), nrow = n, ncol = length(a_0))
  cov.os <- stats::cov(DnoMat) / n

  # output
  out <- list(
    iptw_tmle = list(est = unlist(psi_nStar), cov = cov.t),
    iptw_tmle_nuisance = list(gn = gnStar, QrnStar = QrnStar),
    iptw_os = list(est = unlist(psi.o), cov = cov.os),
    iptw_os_nuisance = list(gn = gn, Qrn = Qrn),
    iptw = list(est = unlist(psi_n)),
    gnMod = NULL, QrnMod = NULL, a_0 = a_0, call = call
  )

  if (returnModels) {
    out$gnMod <- gnMod
    out$QrnMod <- QrnMod
  }
  class(out) <- "adaptive_iptw"
  return(out)
}

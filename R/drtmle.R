#' TMLE estimate of the average treatment effect with doubly-robust inference
#'
#' @param W A \code{data.frame} of named covariates.
#' @param A A \code{numeric} vector of discrete-valued treatment assignment.
#' @param Y A \code{numeric} continuous or binary outcomes.
#' @param DeltaY A \code{numeric} vector of missing outcome indicator (assumed
#'  to be equal to 0 if missing 1 if observed).
#' @param DeltaA A \code{numeric} vector of missing treatment indicator (assumed
#'  to be equal to 0 if missing 1 if observed).
#' @param a_0 A \code{numeric} vector of fixed treatment values at which to
#'  return marginal mean estimates.
#' @param family A \code{family} object equal to either \code{binomial()} or
#'  \code{gaussian()}, to be passed to the \code{SuperLearner} or \code{glm}
#'  function.
#' @param stratify A \code{boolean} indicating whether to estimate the outcome
#'  regression separately for different values of \code{A} (if \code{TRUE}) or
#'  to pool across \code{A} (if \code{FALSE}).
#' @param SL_Q A vector of characters or a list describing the Super Learner
#'  library to be used for the outcome regression. See
#'  \code{\link[SuperLearner]{SuperLearner}} for details.
#' @param SL_g A vector of characters describing the super learner library to be
#'  used for each of the propensity score regressions (\code{DeltaA}, \code{A},
#'  and \code{DeltaY}). To use the same library for each of the regressions (or
#'  if there is no missing data in \code{A} nor \code{Y}), a single library may
#'  be input. See \code{\link[SuperLearner]{SuperLearner}} for details on how
#'  super learner libraries can be specified.
#' @param SL_Qr A vector of characters or a list describing the Super Learner
#'  library to be used for the reduced-dimension outcome regression.
#' @param SL_gr A vector of characters or a list describing the Super Learner
#'  library to be used for the reduced-dimension propensity score.
#' @param n_SL Number of repeated Super Learners to run (default 1) for the
#'  each nuisance parameter. Repeat Super Learners more times to obtain more stable
#'  inference.
#' @param glm_Q A character describing a formula to be used in the call to
#'  \code{glm} for the outcome regression. Ignored if \code{SL_Q!=NULL}.
#' @param glm_g A list of characters describing the formulas to be used
#'  for each of the propensity score regressions (\code{DeltaA}, \code{A}, and
#'  \code{DeltaY}). To use the same formula for each of the regressions (or if
#'  there are no missing data in \code{A} nor \code{Y}), a single character
#'  formula may be input. In general the formulas can reference any variable in 
#'  \code{colnames(W)}, unless \code{adapt_g = TRUE} in which case the formulas
#'  should reference variables \code{QaW} where \code{a} takes values in \code{a_0}.
#' @param glm_Qr A character describing a formula to be used in the call to
#'  \code{glm} for reduced-dimension outcome regression. Ignored if
#'  \code{SL_Qr!=NULL}. The formula should use the variable name \code{'gn'}.
#' @param glm_gr A character describing a formula to be used in the call to
#'  \code{glm} for the reduced-dimension propensity score. Ignored if
#'  \code{SL_gr!=NULL}. The formula should use the variable name \code{'Qn'} and
#'  \code{'gn'} if \code{reduction='bivariate'} and \code{'Qn'} otherwise.
#' @param adapt_g A boolean indicating whether the propensity score should be 
#'  outcome adaptive. If \code{TRUE} then the propensity score is estimated as the
#'  regression of \code{A} onto covariates \code{QaW} for \code{a} in each value
#'  contained in \code{a_0}. See vignette for more details. 
#' @param guard A character vector indicating what pattern of misspecifications
#'  to guard against. If \code{guard} contains \code{"Q"}, then the TMLE guards
#'  against misspecification of the outcome regression by estimating the
#'  reduced-dimension outcome regression specified by \code{glm_Qr} or
#'  \code{SL_Qr}. If \code{guard} contains \code{"g"} then the TMLE
#'  (additionally) guards against misspecification of the propensity score by
#'  estimating the reduced-dimension propensity score specified by \code{glm_gr}
#'  or \code{SL_gr}. If \code{guard} is set to \code{NULL}, then only standard TMLE
#'  and one-step estimators are computed.
#' @param reduction A character equal to \code{"univariate"} for a univariate
#'  misspecification correction (default) or \code{"bivariate"} for the
#'  bivariate version.
#' @param returnModels A boolean indicating whether to return model fits for the
#'  outcome regression, propensity score, and reduced-dimension regressions.
#' @param maxIter A numeric that sets the maximum number of iterations the TMLE
#'  can perform in its fluctuation step.
#' @param tolIC A numeric that defines the stopping criteria based on the
#'  empirical mean of the influence function.
#' @param tolg A numeric indicating the minimum value for estimates of the
#'  propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param Qsteps A numeric equal to 1 or 2 indicating whether the fluctuation
#'  submodel for the outcome regression should be fit using a single
#'  minimization (\code{Qsteps = 1}) or a backfitting-type minimization
#'  (\code{Qsteps=2}). The latter was found to be more stable in simulations and
#'  is the default.
#' @param cvFolds A numeric equal to the number of folds to be used in
#'  cross-validated fitting of nuisance parameters. If \code{cvFolds = 1}, no
#'  cross-validation is used. Alternatively, \code{cvFolds} may be entered as a
#'  vector of fold assignments for observations, in which case its length should
#'  be the same length as \code{Y}.
#' @param parallel A boolean indicating whether to use parallelization based on
#'  \code{future} when estimating nuisance parameters. Only useful if
#'  \code{cvFolds > 1}. By default, a \code{multiprocess} evaluation scheme is
#'  invoked, using forked R processes (if supported on the OS) and background R
#'  sessions otherwise. Users may also register their own backends using the
#'  \code{future.batchtools} package.
#' @param future_hpc A character string identifying a high-performance computing
#'  backend to be used with parallelization. This should match exactly one of
#'  the options available from the \code{future.batchtools} package.
#' @param Qn An optional list of outcome regression estimates. If specified, the
#'  function will ignore the nuisance parameter estimation specified by
#'  \code{SL_Q} and \code{glm_Q}. The entries in the list should correspond to
#'  the outcome regression evaluated at \code{A} and the observed values of
#'  \code{W}, with order determined by the input to \code{a_0} (e.g., if
#'  \code{a_0 = c(0, 1)} then \code{Qn[[1]]} should be outcome regression at
#'  \code{A} = 0 and \code{Qn[[2]]} should be outcome regression at
#'  \code{A} = 1).
#' @param gn An optional list of propensity score estimates. If specified, the
#'  function will ignore the nuisance parameter estimation specified by
#'  \code{SL_g} and \code{glm_g}. The entries in the list should correspond to
#'  the propensity for the observed values of \code{W}, with order determined by
#'  the input to \code{a_0} (e.g., if \code{a_0 = c(0,1)} then \code{gn[[1]]}
#'  should be propensity of \code{A} = 0 and \code{gn[[2]]} should be propensity
#'  of \code{A} = 1).
#' @param use_future Boolean indicating whether to use \code{future_lapply} or
#' instead to just use lapply. The latter can be easier to run down errors.
#' @param ... Other options (not currently used).
#'
#' @return An object of class \code{"drtmle"}.
#' \describe{
#'  \item{\code{drtmle}}{A \code{list} of doubly-robust point estimates and
#'        a doubly-robust covariance matrix}
#'  \item{\code{nuisance_drtmle}}{A \code{list} of the final TMLE estimates of
#'        the outcome regression (\code{$QnStar}), propensity score
#'        (\code{$gnStar}), and reduced-dimension regressions (\code{$QrnStar},
#'        \code{$grnStar}) evaluated at the observed data values.}
#'  \item{\code{ic_drtmle}}{A \code{list} of the empirical mean of the efficient
#'        influence function (\code{$eif}) and the extra pieces of the influence
#'        function resulting from misspecification. All should be smaller than
#'        \code{tolIC} (unless \code{maxIter} was reached first). Also includes
#'        a matrix of the influence function values at the estimated nuisance
#'        parameters evaluated at the observed data.}
#'  \item{\code{aiptw_c}}{A \code{list} of doubly-robust point estimates and
#'        a non-doubly-robust covariance matrix. Theory does not guarantee
#'        performance of inference for these estimators, but simulation studies
#'        showed they often perform adequately.}
#'  \item{\code{nuisance_aiptw}}{A \code{list} of the initial estimates of the
#'        outcome regression, propensity score, and reduced-dimension
#'        regressions evaluated at the observed data values.}
#'  \item{\code{tmle}}{A \code{list} of doubly-robust point estimates and
#'        non-doubly-robust covariance for the standard TMLE estimator.}
#'  \item{\code{aiptw}}{A \code{list} of doubly-robust point estimates and
#'        non-doubly-robust covariance matrix for the standard AIPTW estimator.}
#'  \item{\code{gcomp}}{A \code{list} of non-doubly-robust point estimates and
#'        non-doubly-robust covariance matrix for the standard G-computation
#'        estimator. If super learner is used there is no guarantee of correct
#'        inference for this estimator.}
#'  \item{\code{QnMod}}{The fitted object for the outcome regression. Returns
#'        \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{gnMod}}{The fitted object for the propensity score. Returns
#'        \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{QrnMod}}{The fitted object for the reduced-dimension regression
#'        that guards against misspecification of the outcome regression.
#'        Returns \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{grnMod}}{The fitted object for the reduced-dimension regression
#'        that guards against misspecification of the propensity score. Returns
#'        \code{NULL} if \code{returnModels = FALSE}.}
#'  \item{\code{a_0}}{The treatment levels that were requested for computation
#'        of covariate-adjusted means.}
#' }
#'
#' @importFrom future plan sequential multiprocess
#' @importFrom future.apply future_lapply
#' @importFrom doFuture registerDoFuture
#' @importFrom future.batchtools batchtools_slurm batchtools_lsf batchtools_sge
#'  batchtools_torque batchtools_openlava
#' @importFrom stats cov
#'
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
#' # A quick example of drtmle:
#' # We note that more flexible super learner libraries
#' # are available, and that we recommend the user use more flexible
#' # libraries for SL_Qr and SL_gr for general use.
#' fit1 <- drtmle(
#'   W = W, A = A, Y = Y, a_0 = c(1, 0),
#'   family = binomial(),
#'   stratify = FALSE,
#'   SL_Q = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#'   SL_g = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#'   SL_Qr = "SL.glm",
#'   SL_gr = "SL.glm", maxIter = 1
#' )
drtmle <- function(Y, A, W,
                   DeltaA = as.numeric(!is.na(A)),
                   DeltaY = as.numeric(!is.na(Y)),
                   a_0 = unique(A[!is.na(A)]),
                   family = if (all(Y %in% c(0, 1))) {
                     stats::binomial()
                   } else {
                     stats::gaussian()
                   },
                   stratify = TRUE,
                   SL_Q = NULL,
                   SL_g = NULL,
                   SL_Qr = NULL,
                   SL_gr = NULL,
                   n_SL = 1,
                   glm_Q = NULL, glm_g = NULL,
                   glm_Qr = NULL, glm_gr = NULL,
                   adapt_g = FALSE, 
                   guard = c("Q", "g"),
                   reduction = "univariate",
                   returnModels = FALSE,
                   cvFolds = 1,
                   maxIter = 3,
                   tolIC = 1 / length(Y),
                   tolg = 1e-2,
                   verbose = FALSE,
                   Qsteps = 2,
                   parallel = FALSE,
                   future_hpc = NULL,
                   Qn = NULL,
                   gn = NULL,
                   use_future = TRUE,
                   ...) {
  call <- match.call()
  n <- length(Y)
  # save user input Qn and gn, because these values
  # get overwritten later.
  Qn_user <- !is.null(Qn)
  gn_user <- !is.null(gn)

  # check if any additional targeting is performed
  extra_targeting <- !is.null(guard)
  # if using adaptive propensity, extra targeting not needed
  if(adapt_g){
    extra_targeting <- FALSE
    guard <- NULL
  }
  # if cvFolds non-null split data into cvFolds pieces
  validRows <- make_validRows(cvFolds, n = length(Y), n_SL = n_SL)
  if (n_SL > 1) {
    validRows <- rep(validRows, n_SL)
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
  # estimate outcome regression
  # -------------------------------
  if (is.null(Qn)) {
    if (use_future) {
      QnOut <- future.apply::future_lapply(
        X = validRows, FUN = estimateQ,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        verbose = verbose,
        returnModels = returnModels,
        SL_Q = SL_Q, a_0 = a_0,
        stratify = stratify,
        glm_Q = glm_Q,
        family = family
      )
    } else {
      QnOut <- lapply(
        X = validRows, FUN = estimateQ,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        verbose = verbose,
        returnModels = returnModels,
        SL_Q = SL_Q, a_0 = a_0,
        stratify = stratify,
        glm_Q = glm_Q,
        family = family
      )
    }
    # re-order predictions
    Qn <- reorder_list(QnOut,
      a_0 = a_0, validRows = validRows,
      n_SL = n_SL, n = n
    )

    # obtain list of outcome regression fits
    QnMod <- extract_models(QnOut)
  }

  # -------------------------------
  # estimate propensity score
  # -------------------------------
  if (is.null(gn)) {
    if (use_future) {
      gnOut <- future.apply::future_lapply(
        X = validRows, FUN = estimateG, A = A,
        W = W, DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, verbose = verbose,
        stratify = stratify,
        returnModels = returnModels, SL_g = SL_g,
        glm_g = glm_g, a_0 = a_0, 
        Qn = Qn, adapt_g = adapt_g
      )
    } else {
      gnOut <- lapply(
        X = validRows, FUN = estimateG, A = A,
        W = W, DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, verbose = verbose,
        stratify = stratify,
        returnModels = returnModels, SL_g = SL_g,
        glm_g = glm_g, a_0 = a_0,
        Qn = Qn, adapt_g = adapt_g
      )
    }
    # re-order predictions
    gn <- reorder_list(gnOut,
      a_0 = a_0, validRows = validRows,
      n_SL = n_SL, n = n
    )
    # obtain list of propensity score fits
    gnMod <- extract_models(gnOut)
  } else {
    # truncate too-small predictions
    gn <- lapply(gn, function(g) {
      g[g < tolg] <- tolg
      g
    })
  }
  # naive g-computation estimate
  psi_n <- lapply(Qn, mean)

  # estimate influence function
  Dno <- eval_Dstar(
    A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
    Qn = Qn, gn = gn, psi_n = psi_n, a_0 = a_0
  )

  # estimate bias correction
  PnDn <- lapply(Dno, mean)

  # additional bias terms
  Dngo <- rep(0, n)
  DnQo <- rep(0, n)
  PnDQn <- PnDgn <- 0

  if ("Q" %in% guard) {
    if (use_future) {
      QrnOut <- future.apply::future_lapply(
        X = validRows, FUN = estimateQrn,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        Qn = Qn, gn = gn, glm_Qr = glm_Qr,
        family = stats::gaussian(), SL_Qr = SL_Qr,
        a_0 = a_0, returnModels = returnModels
      )
    } else {
      QrnOut <- lapply(
        X = validRows, FUN = estimateQrn,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        Qn = Qn, gn = gn, glm_Qr = glm_Qr,
        family = stats::gaussian(), SL_Qr = SL_Qr,
        a_0 = a_0, returnModels = returnModels
      )
    }
    # re-order predictions
    Qrn <- reorder_list(QrnOut,
      a_0 = a_0, validRows = validRows,
      n_SL = n_SL, n = n
    )

    # obtain list of reduced dimension regression fits
    QrnMod <- extract_models(QrnOut)

    Dngo <- eval_Dstar_g(
      A = A, DeltaY = DeltaY, DeltaA = DeltaA, Qrn = Qrn,
      gn = gn, a_0 = a_0
    )
    PnDgn <- lapply(Dngo, mean)
  } else {
    Qrn <- NULL
  }
  if ("g" %in% guard) {
    if (use_future) {
      grnOut <- future.apply::future_lapply(
        X = validRows, FUN = estimategrn,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, Qn = Qn, gn = gn,
        glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
        reduction = reduction,
        returnModels = returnModels
      )
    } else {
      grnOut <- lapply(
        X = validRows, FUN = estimategrn,
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, Qn = Qn, gn = gn,
        glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
        reduction = reduction,
        returnModels = returnModels
      )
    }
    # re-order predictions
    grn <- reorder_list(grnOut,
      a_0 = a_0, validRows = validRows,
      grn_ind = TRUE,
      n_SL = n_SL, n = n
    )

    # obtain list of outcome regression fits
    grnMod <- extract_models(grnOut)

    # evaluate extra piece of influence function
    DnQo <- eval_Dstar_Q(
      A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
      Qn = Qn, grn = grn, gn = gn, a_0 = a_0,
      reduction = reduction
    )
    PnDQn <- lapply(DnQo, mean)
  } else {
    grn <- NULL
  }

  # one step estimates
  psi_o1 <- mapply(
    a = psi_n, b = PnDn, SIMPLIFY = FALSE,
    FUN = function(a, b) {
      a + b
    }
  )
  psi_o <- mapply(
    a = psi_n, b = PnDn, c = PnDQn, d = PnDgn, SIMPLIFY = FALSE,
    FUN = function(a, b, c, d) {
      a + b - c - d
    }
  )

  # covariance for one step
  Dno1Mat <- matrix(unlist(Dno), nrow = n, ncol = length(a_0))
  DnoMat <- matrix(
    unlist(Dno) - unlist(DnQo) - unlist(Dngo),
    nrow = n, ncol = length(a_0)
  )

  cov_o1 <- stats::cov(Dno1Mat) / n
  cov_o <- stats::cov(DnoMat) / n

  # initialize fluctuations
  QnStar <- Qn
  if ("g" %in% guard) {
    grnStar <- grn
    PnDQnStar <- Inf
  } else {
    grnStar <- NULL
    PnDQnStar <- 0
  }
  if ("Q" %in% guard) {
    QrnStar <- Qrn
    PnDgnStar <- Inf
  } else {
    QrnStar <- NULL
    PnDgnStar <- 0
  }
  gnStar <- gn
  PnDnoStar <- Inf
  ct <- 0

  if (extra_targeting) {
    # fluctuate
    while (max(abs(c(unlist(PnDQnStar), unlist(PnDgnStar), unlist(PnDnoStar)))) >
      tolIC & ct < maxIter) {
      ct <- ct + 1

      # re-estimate Qrn
      if ("Q" %in% guard) {
        # fluctuate gnStar
        gnStarOut <- fluctuateG(
          Y = Y, A = A, W = W, DeltaA = DeltaA,
          DeltaY = DeltaY, a_0 = a_0, tolg = tolg,
          gn = gnStar, Qrn = QrnStar
        )
        gnStar <- lapply(gnStarOut, function(x) {
          unlist(x$est)
        })
      }

      # fluctuate QnStar
      if ("g" %in% guard) {
        if (use_future) {
          grnStarOut <- future.apply::future_lapply(
            X = validRows, FUN = estimategrn,
            Y = Y, A = A, W = W,
            DeltaA = DeltaA, DeltaY = DeltaY,
            tolg = tolg, Qn = QnStar,
            gn = gnStar, glm_gr = glm_gr,
            SL_gr = SL_gr, a_0 = a_0,
            reduction = reduction,
            returnModels = returnModels
          )
        } else {
          grnStarOut <- lapply(
            X = validRows, FUN = estimategrn,
            Y = Y, A = A, W = W,
            DeltaA = DeltaA, DeltaY = DeltaY,
            tolg = tolg, Qn = QnStar,
            gn = gnStar, glm_gr = glm_gr,
            SL_gr = SL_gr, a_0 = a_0,
            reduction = reduction,
            returnModels = returnModels
          )
        }
        # re-order predictions
        grnStar <- reorder_list(grnStarOut,
          a_0 = a_0, validRows = validRows,
          grn_ind = TRUE, n_SL = n_SL, n = n
        )

        # obtain list of outcome regression fits
        grnMod <- extract_models(grnStarOut)

        if (Qsteps == 1) {
          QnStarOut <- fluctuateQ(
            Y = Y, A = A, W = W, DeltaA = DeltaA,
            DeltaY = DeltaY, a_0 = a_0, Qn = QnStar,
            gn = gnStar, grn = grnStar,
            reduction = reduction
          )
          QnStar <- lapply(QnStarOut, function(x) {
            unlist(x$est)
          })
        } else if (Qsteps == 2) {
          # do the extra targeting
          QnStarOut2 <- fluctuateQ2(
            Y = Y, A = A, W = W, DeltaA = DeltaA,
            DeltaY = DeltaY, a_0 = a_0, Qn = QnStar,
            gn = gnStar, grn = grnStar,
            reduction = reduction
          )
          QnStar <- lapply(QnStarOut2, function(x) {
            unlist(x[[1]])
          })

          # do the usual targeting
          QnStarOut1 <- fluctuateQ1(
            Y = Y, A = A, W = W, DeltaA = DeltaA,
            DeltaY = DeltaY, a_0 = a_0, Qn = QnStar,
            gn = gnStar
          )
          QnStar <- lapply(QnStarOut1, function(x) {
            unlist(x[[1]])
          })
        }
      } else {
        QnStarOut <- fluctuateQ1(
          Y = Y, A = A, W = W, DeltaA = DeltaA,
          DeltaY = DeltaY, a_0 = a_0, Qn = QnStar,
          gn = gnStar
        )
        QnStar <- lapply(QnStarOut, function(x) {
          unlist(x[[1]])
        })
      }

      if ("Q" %in% guard) {
        if (use_future) {
          QrnStarOut <- future.apply::future_lapply(
            X = validRows, FUN = estimateQrn,
            Y = Y, A = A, W = W,
            DeltaA = DeltaA, DeltaY = DeltaY,
            Qn = QnStar, gn = gnStar,
            glm_Qr = glm_Qr,
            family = stats::gaussian(),
            SL_Qr = SL_Qr, a_0 = a_0,
            returnModels = returnModels
          )
        } else {
          QrnStarOut <- lapply(
            X = validRows, FUN = estimateQrn,
            Y = Y, A = A, W = W,
            DeltaA = DeltaA, DeltaY = DeltaY,
            Qn = QnStar, gn = gnStar,
            glm_Qr = glm_Qr,
            family = stats::gaussian(),
            SL_Qr = SL_Qr, a_0 = a_0,
            returnModels = returnModels
          )
        }
        # re-order predictions
        QrnStar <- reorder_list(QrnStarOut,
          a_0 = a_0, validRows = validRows,
          n_SL = n_SL, n = n
        )
        # obtain list of reduced dimension regression fits
        QrnMod <- extract_models(QrnStarOut)
      }

      # tmle estimates
      psi_t <- lapply(QnStar, mean)

      # calculate influence functions
      DnoStar <- eval_Dstar(
        A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
        Qn = QnStar, gn = gnStar, psi_n = psi_t, a_0 = a_0
      )
      PnDnoStar <- lapply(DnoStar, mean)

      if ("g" %in% guard) {
        DnQoStar <- eval_Dstar_Q(
          A = A, Y = Y, DeltaY = DeltaY,
          DeltaA = DeltaA, Qn = QnStar, grn = grnStar, gn = gn,
          a_0 = a_0, reduction = reduction
        )
        PnDQnStar <- lapply(DnQoStar, mean)
      }
      if ("Q" %in% guard) {
        DngoStar <- eval_Dstar_g(
          A = A, DeltaY = DeltaY, DeltaA = DeltaA,
          Qrn = QrnStar, gn = gnStar, a_0 = a_0
        )
        PnDgnStar <- lapply(DngoStar, mean)
      }
      if (verbose) {
        cat("Mean of IC       =", round(c(
          unlist(PnDnoStar),
          unlist(PnDQnStar),
          unlist(PnDgnStar)
        ), 10), "\n")
      }
    }
  }

  # standard tmle fluctuations
  QnStar1Out <- fluctuateQ1(
    Y = Y, A = A, W = W, DeltaA = DeltaA,
    DeltaY = DeltaY, Qn = Qn, gn = gn, a_0 = a_0
  )
  QnStar1 <- lapply(QnStar1Out, function(x) {
    unlist(x[[1]])
  })

  # tmle estimates
  psi_t1 <- lapply(QnStar1, mean)
  if (extra_targeting) {
    psi_t <- lapply(QnStar, mean)
  } else {
    psi_t <- psi_t1
  }


  # covariance for tmle
  Dno1Star <- eval_Dstar(
    A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY,
    Qn = QnStar1, gn = gn, psi_n = psi_t1, a_0 = a_0
  )
  Dno1StarMat <- matrix(unlist(Dno1Star), nrow = n, ncol = length(a_0))
  cov_t1 <- stats::cov(Dno1StarMat) / n

  # covariance for drtmle
  DnoStar <- eval_Dstar(
    A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY,
    Qn = QnStar, gn = gnStar, psi_n = psi_t, a_0 = a_0
  )
  PnDnoStar <- lapply(DnoStar, mean)

  DnQoStar <- rep(list(rep(0, n)), length(a_0))
  DngoStar <- rep(list(rep(0, n)), length(a_0))

  if ("g" %in% guard) {
    DnQoStar <- eval_Dstar_Q(
      A = A, Y = Y, DeltaY = DeltaY,
      DeltaA = DeltaA, Qn = QnStar, grn = grnStar, gn = gn,
      a_0 = a_0, reduction = reduction
    )
    PnDQnStar <- lapply(DnQoStar, mean)
  }
  if ("Q" %in% guard) {
    DngoStar <- eval_Dstar_g(
      A = A, DeltaY = DeltaY, DeltaA = DeltaA,
      Qrn = QrnStar, gn = gnStar, a_0 = a_0
    )
    PnDgnStar <- lapply(DngoStar, mean)
  }

  DnoStarMat <- matrix(
    unlist(DnoStar) - unlist(DnQoStar) - unlist(DngoStar),
    nrow = n, ncol = length(a_0)
  )
  cov_t <- stats::cov(DnoStarMat) / n

  out <- list(
    drtmle = list(est = unlist(psi_t), cov = cov_t),
    nuisance_drtmle = list(
      QnStar = QnStar, gnStar = gnStar,
      QrnStar = QrnStar, grnStar = grn,
      meanIC = unlist(c(
        PnDnoStar, PnDQnStar,
        PnDgnStar
      ))
    ),
    ic_drtmle = list(
      mean_eif = PnDnoStar, mean_missQ = PnDgnStar,
      mean_missg = PnDQnStar, ic = DnoStarMat
    ),
    aiptw_c = list(est = unlist(psi_o), cov = cov_o),
    nuisance_aiptw_c = list(Qn = Qn, gn = gn, Qrn = Qrn, grn = grn),
    tmle = list(est = unlist(psi_t1), cov = cov_t1),
    aiptw = list(est = unlist(psi_o1), cov = cov_o1),
    gcomp = list(est = unlist(psi_n), cov = cov_o1),
    QnMod = NULL, gnMod = NULL, QrnMod = NULL, grnMod = NULL,
    a_0 = a_0, call = call
  )

  # tack on models if requested
  if (returnModels) {
    if (!Qn_user) {
      out$QnMod <- QnMod
    }
    if (!gn_user) {
      out$gnMod <- gnMod
    }
    if ("Q" %in% guard) {
      out$QrnMod <- QrnMod
    }
    if ("g" %in% guard) {
      out$grnMod <- grnMod
    }
  }
  class(out) <- "drtmle"
  return(out)
}

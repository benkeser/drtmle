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
#' @param avg_over If multiple Super Learners are run, on which scale should the
#' results be aggregated. Options include: \code{"SL"} = 
#' repeated nuisance parameter estimates are averaged before subsequently 
#' generating a single vector of point estimates based on the averaged models;
#' \code{"drtmle"} = repeated vectors of point estimates are generated and 
#' averaged. Both can be specified, recognizing that this adds considerable
#' computational expense. In this case, the final estimates are the average 
#' of \code{n_SL} point estimates where each is built by averaging \code{n_SL} 
#' fits. If \code{NULL}, no averaging is performed (in which case \code{n_SL} 
#' should be set equal to 1). 
#' @param se_cv Should cross-validated nuisance parameter estimates be used 
#' for computing standard errors? 
#' Options are \code{"none"} = no cross-validation is performed; \code{"partial"} = 
#' only applicable if Super Learner is used for nuisance parameter estimates; 
#' \code{"full"} = full cross-validation is performed. See vignette for further 
#' details. Ignored if \code{cvFolds > 1}, since then
#' cross-validated nuisance parameter estimates are used by default and it is 
#' assumed that you want full cross-validated standard errors. 
#' @param se_cvFolds If cross-validated nuisance parameter estimates are used
#' to compute standard errors, how many folds should be used in this computation. 
#' If \code{se_cv = "partial"}, then this option sets the number of folds used
#' by the \code{SuperLearner} fitting procedure. 
#' @param targeted_se A boolean indicating whether the targeted nuisance 
#' parameters should be used in standard error computation or the initial 
#' estimators. If \code{se_cv} is not set to \code{"none"}, this option is 
#' ignored and standard errors are computed based on non-targeted, cross-validated 
#' nuisance parameter fits. 
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
#' @param returnNuisance A boolean indicating whether to return the estimated 
#'  nuisance regressions evaluated on the observed data. Defaults to \code{TRUE}. 
#'  If \code{n_SL} is large and \code{"drtmle"} is in \code{avg_over}, then 
#'  consider setting to \code{FALSE} in order to reduce size of resultant object.
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
#' @importFrom future.apply future_lapply
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
                   stratify = FALSE,
                   SL_Q = NULL,
                   SL_g = NULL,
                   SL_Qr = NULL,
                   SL_gr = NULL,
                   n_SL = 1,
                   avg_over = "drtmle",
                   se_cv = "none", 
                   se_cvFolds = ifelse(se_cv == "partial", 10, 1),
                   targeted_se = se_cv != "partial", 
                   glm_Q = NULL, glm_g = NULL,
                   glm_Qr = NULL, glm_gr = NULL,
                   adapt_g = FALSE, 
                   guard = c("Q", "g"),
                   reduction = "univariate",
                   returnModels = FALSE,
                   returnNuisance = TRUE, 
                   cvFolds = 1,
                   maxIter = 3,
                   tolIC = 1 / length(Y),
                   tolg = 1e-2,
                   verbose = FALSE,
                   Qsteps = 2,
                   Qn = NULL,
                   gn = NULL,
                   use_future = FALSE,
                   ...) {
  call <- match.call()
  n <- length(Y)
  # save user input Qn and gn, because these values
  # get overwritten later.
  Qn_user <- !is.null(Qn)
  if(Qn_user){
    Qn_se <- Qn_user
  }
  gn_user <- !is.null(gn)
  if(gn_user){
    gn_se <- gn_user
  }

  # check if any additional targeting is performed
  extra_targeting <- !is.null(guard)
  # if using adaptive propensity, extra targeting not needed
  if(adapt_g){
    extra_targeting <- FALSE
    guard <- NULL
  }

  if(n_SL == 1 | !("drtmle" %in% avg_over)){ 
    n_loops <- 1
    n_SL_avg <- n_SL
  }else{
    n_loops <- n_SL
    if("SL" %in% avg_over){
      n_SL_avg <- n_SL
    }else{
      n_SL_avg <- 1
    }
  }

  # if cvtmle requested, then cross-validated standard errors
  # will be computed by default with a normal run through the code
  # so we force se_cv to be "none"
  if(cvFolds > 1){
    se_cv <- "none"
  }

  # make empty holder objects
  nuisance_drtmle_list <- vector(mode = "list", length = n_loops) 
  nuisance_aiptw_c_list <- vector(mode = "list", length = n_loops)
  ic_drtmle_list <- vector(mode = "list", length = n_loops) 
  QnMod_list <- vector(mode = "list", length = n_loops)
  gnMod_list <- vector(mode = "list", length = n_loops)
  QrnMod_list <- vector(mode = "list", length = n_loops)
  grnMod_list <- vector(mode = "list", length = n_loops)
  drtmle_list <- vector(mode = "list", length = n_loops)
  aiptw_c_list <- vector(mode = "list", length = n_loops)
  tmle_list <- vector(mode = "list", length = n_loops)
  aiptw_list <- vector(mode = "list", length = n_loops)
  gcomp_list <- vector(mode = "list", length = n_loops)
  validRows_list <- vector(mode = "list", length = n_loops)

  for(i in seq_len(n_loops)){
    # if cvFolds non-null split data into cvFolds pieces
    validRows_list[[i]] <- make_validRows(cvFolds, n = length(Y))
    if (n_SL > 1 & ("SL" %in% avg_over)) {
      validRows_list[[i]] <- rep(validRows_list[[i]], n_SL)
    }

    # -------------------------------
    # estimate outcome regression
    # -------------------------------
    if (is.null(Qn)) {
      QnOut <- estimateQ_loop(
        validRows = validRows_list[[i]],
        Y = Y, A = A, W = W, 
        DeltaA = DeltaA, DeltaY = DeltaY,
        verbose = verbose,
        returnModels = returnModels,
        SL_Q = SL_Q, a_0 = a_0,
        stratify = stratify,
        glm_Q = glm_Q,
        family = family,
        use_future = use_future,
        se_cv = se_cv, se_cvFolds = se_cvFolds
      )
      
      # re-order predictions 
      Qn <- reorder_list(QnOut,
        a_0 = a_0, validRows = validRows_list[[i]],
        n_SL = n_SL_avg, n = n
      )

      # obtain list of outcome regression fits
      QnMod <- extract_models(QnOut)

      if(se_cv == "full"){
        validRows_se <- make_validRows(se_cvFolds, n = length(Y))
        QnOut_se <- estimateQ_loop(
          validRows = validRows_se,
          Y = Y, A = A, W = W, 
          DeltaA = DeltaA, DeltaY = DeltaY,
          verbose = verbose,
          returnModels = FALSE, # don't keep models here
          SL_Q = SL_Q, a_0 = a_0,
          stratify = stratify,
          glm_Q = glm_Q,
          family = family,
          use_future = use_future,
          se_cv = se_cv, se_cvFolds = se_cvFolds
        )
        # re-order predictions
        Qn_se <- reorder_list(QnOut_se,
          a_0 = a_0, validRows = validRows_se,
          n_SL = n_SL_avg, n = n
        )
      }else if(se_cv == "partial"){
        Qn_se <- reorder_list(QnOut,
          a_0 = a_0, validRows = validRows_list[[i]],
          n_SL = n_SL_avg, n = n, for_se_cv = TRUE
        )
      }else{
        Qn_se <- Qn
      }
    }

    # -------------------------------
    # estimate propensity score
    # -------------------------------
    if (is.null(gn)) {
      gnOut <- estimateG_loop(
        validRows = validRows_list[[i]], A = A,
        W = W, DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, verbose = verbose,
        stratify = stratify,
        returnModels = returnModels, SL_g = SL_g,
        glm_g = glm_g, a_0 = a_0, 
        Qn = Qn, adapt_g = adapt_g,
        use_future = use_future,
        se_cv = se_cv, se_cvFolds = se_cvFolds
      )

      # re-order predictions
      gn <- reorder_list(gnOut,
        a_0 = a_0, validRows = validRows_list[[i]],
        n_SL = n_SL_avg, n = n
      )
      # obtain list of propensity score fits
      gnMod <- extract_models(gnOut)

      if(se_cv == "full"){
        # validRows_se defined above
        gnOut_se <- estimateG_loop(
          validRows = validRows_se, A = A,
          W = W, DeltaA = DeltaA, DeltaY = DeltaY,
          tolg = tolg, verbose = verbose,
          stratify = stratify,
          returnModels = returnModels, SL_g = SL_g,
          glm_g = glm_g, a_0 = a_0, 
          Qn = Qn, adapt_g = adapt_g,
          use_future = use_future, se_cv = se_cv,
          se_cvFolds = se_cvFolds
        )
        gn_se <- reorder_list(gnOut_se,
          a_0 = a_0, validRows = validRows_se,
          n_SL = n_SL_avg, n = n
        )
      }else if(se_cv == "partial"){
        gn_se <- reorder_list(gnOut,
          a_0 = a_0, validRows = validRows_list[[i]],
          n_SL = n_SL_avg, n = n, for_se_cv = TRUE
        )
      }else{
        gn_se <- gn
      }
    } else {
      # truncate too-small predictions
      gn <- lapply(gn, function(g) {
        g[g < tolg] <- tolg
        g
      })
      gn_se <- gn
    }
    # naive g-computation estimate
    psi_n <- lapply(Qn, mean)
    psi_n_se <- lapply(Qn_se, mean)
    # estimate influence function
    Dno <- eval_Dstar(
      A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
      Qn = Qn, gn = gn, psi_n = psi_n, a_0 = a_0
    )
    # and the one used for variance computation (possibly same as Dno)
    Dno_se <- eval_Dstar(
      A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
      Qn = Qn_se, gn = gn_se, psi_n = psi_n_se, a_0 = a_0
    )

    # estimate bias correction
    PnDn <- lapply(Dno, mean)

    # additional bias terms
    Dngo <- rep(0, n)
    Dngo_se <- Dngo
    DnQo <- rep(0, n)
    DnQo_se <- DnQo
    PnDQn <- PnDgn <- 0

    if ("Q" %in% guard) {
      QrnOut <- estimateQrn_loop(
        validRows = validRows_list[[i]],
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        Qn = Qn, gn = gn, glm_Qr = glm_Qr,
        family = stats::gaussian(), SL_Qr = SL_Qr,
        a_0 = a_0, returnModels = returnModels,
        use_future = use_future
      )
      # re-order predictions
      Qrn <- reorder_list(QrnOut,
        a_0 = a_0, validRows = validRows_list[[i]],
        n_SL = n_SL_avg, n = n
      )

      # obtain list of reduced dimension regression fits
      QrnMod <- extract_models(QrnOut)

      # extra piece of influence function
      Dngo <- eval_Dstar_g(
        A = A, DeltaY = DeltaY, DeltaA = DeltaA, Qrn = Qrn,
        gn = gn, a_0 = a_0
      )
      
      PnDgn <- lapply(Dngo, mean)

      if(se_cv == "full"){
        QrnOut_se <- estimateQrn_loop(
          validRows = validRows_se, 
          Y = Y, A = A, W = W,
          DeltaA = DeltaA, DeltaY = DeltaY,
          Qn = Qn_se, gn = gn_se, glm_Qr = glm_Qr,
          family = stats::gaussian(), SL_Qr = SL_Qr,
          a_0 = a_0, returnModels = FALSE,
          use_future = use_future
        )

        Qrn_se <- reorder_list(QrnOut_se,
          a_0 = a_0, validRows = validRows_se,
          n_SL = n_SL_avg, n = n
        )
        
        Dngo_se <- eval_Dstar_g(
          A = A, DeltaY = DeltaY, DeltaA = DeltaA, Qrn = Qrn_se,
          gn = gn_se, a_0 = a_0
        )
      }else{
        Qrn_se <- Qrn
        Dngo_se <- Dngo
      }
    } else {
      Qrn <- NULL
      Qrn_se <- NULL
    }


    if ("g" %in% guard) {
      grnOut <- estimategrn_loop(
        validRows = validRows_list[[i]], 
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, Qn = Qn, gn = gn,
        glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
        reduction = reduction,
        returnModels = returnModels,
        use_future = use_future
      )

      # re-order predictions
      grn <- reorder_list(grnOut,
        a_0 = a_0, validRows = validRows_list[[i]],
        grn_ind = TRUE,
        n_SL = n_SL_avg, n = n
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

      if(se_cv == "full"){
        grnOut_se <- estimategrn_loop(
          validRows = validRows_se, 
          Y = Y, A = A, W = W,
          DeltaA = DeltaA, DeltaY = DeltaY,
          tolg = tolg, Qn = Qn_se, gn = gn_se,
          glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
          reduction = reduction,
          returnModels = FALSE,
          use_future = use_future
        )

        grn_se <- reorder_list(grnOut,
          a_0 = a_0, validRows = validRows_se,
          grn_ind = TRUE,
          n_SL = n_SL_avg, n = n
        )
        
        DnQo_se <- eval_Dstar_Q(
          A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA,
          Qn = Qn_se, grn = grn_se, gn = gn_se, a_0 = a_0,
          reduction = reduction
        )
      }else{
        grn_se <- grn
        DnQo_se <- DnQo
      }
    } else {
      grn <- NULL
      grn_se <- NULL
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
    Dno1Mat <- matrix(unlist(Dno_se), nrow = n, ncol = length(a_0))
    DnoMat <- matrix(
      unlist(Dno_se) - unlist(DnQo_se) - unlist(Dngo_se),
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
          grnStarOut <- estimategrn_loop(
            validRows = validRows_list[[i]], 
            Y = Y, A = A, W = W,
            DeltaA = DeltaA, DeltaY = DeltaY,
            tolg = tolg, Qn = QnStar, gn = gnStar,
            glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
            reduction = reduction,
            returnModels = returnModels,
            use_future = use_future
          )
          
          # re-order predictions
          grnStar <- reorder_list(grnStarOut,
            a_0 = a_0, validRows = validRows_list[[i]],
            grn_ind = TRUE, n_SL = n_SL_avg, n = n
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
              X = validRows_list[[i]], FUN = estimateQrn,
              Y = Y, A = A, W = W,
              DeltaA = DeltaA, DeltaY = DeltaY,
              Qn = QnStar, gn = gnStar,
              glm_Qr = glm_Qr,
              family = stats::gaussian(),
              SL_Qr = SL_Qr, a_0 = a_0,
              returnModels = returnModels,
              future.seed = TRUE
            )
          } else {
            QrnStarOut <- lapply(
              X = validRows_list[[i]], FUN = estimateQrn,
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
            a_0 = a_0, validRows = validRows_list[[i]],
            n_SL = n_SL_avg, n = n
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
    if(targeted_se){
      Dno1Star <- eval_Dstar(
        A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY,
        Qn = QnStar1, gn = gn, psi_n = psi_t1, a_0 = a_0
      )
      Dno1StarMat <- matrix(unlist(Dno1Star), nrow = n, ncol = length(a_0))
      cov_t1 <- stats::cov(Dno1StarMat) / n
    }else{
      cov_t1 <- cov_o1
    }

    if(targeted_se){
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
    }else{
      cov_t <- cov_o
    }
    # add results to relevant lists to store relevant output 
    # from this iteration of the for loop
    nuisance_drtmle_list[[i]] <- list(
      QnStar = QnStar, gnStar = gnStar,
      QrnStar = QrnStar, grnStar = grn,
      meanIC = unlist(c(
        PnDnoStar, PnDQnStar,
        PnDgnStar
      ))
    )
    nuisance_aiptw_c_list[[i]] = list(Qn = Qn, gn = gn, 
                                      Qrn = Qrn, grn = grn)

    if(targeted_se){
      ic_drtmle_list[[i]] = list(
        mean_eif = PnDnoStar, mean_missQ = PnDgnStar,
        mean_missg = PnDQnStar, ic = DnoStarMat
      )
    }else{
      ic_drtmle_list[[i]] = list(
        mean_eif = PnDnoStar, mean_missQ = PnDgnStar,
        mean_missg = PnDQnStar, ic = DnoMat
      )
    }
    if(returnModels){
      if(!Qn_user){
        QnMod_list[[i]] <- QnMod
      }
      if(!gn_user){
        gnMod_list[[i]] <- gnMod
      }
      if ("Q" %in% guard) {
        QrnMod_list[[i]] <- QrnMod
      }
      if ("g" %in% guard) {
        grnMod_list[[i]] <- grnMod
      }
    }
    drtmle_list[[i]] <- list(est = unlist(psi_t), cov = cov_t)
    tmle_list[[i]] <- list(est = unlist(psi_t1), cov = cov_t1)
    aiptw_list[[i]] <- list(est = unlist(psi_o1), cov = cov_o1)
    gcomp_list[[i]] <- list(est = unlist(psi_n), cov = cov_o1)
    aiptw_c_list[[i]] <- list(est = unlist(psi_o), cov = cov_o)
  }


  #   QnMod, gnMod, QrnMod, grnMod
  #   drtmle -- eventually avg over to get final point estimates
  #   aiptw_c -- eventually avg over 
  #   tmle -- eventually avg over
  #   aiptw -- eventually avg over
  #   gcomp -- eventually avg over
  #   ic_drtmle
  #   average over these guys!!!
  #   validRows -- should this be added to the output?

  out <- list(
    drtmle = average_est_cov_list(drtmle_list),
    nuisance_drtmle = NULL,
    ic_drtmle = average_ic_list(ic_drtmle_list),
    aiptw_c = average_est_cov_list(aiptw_c_list),
    nuisance_aiptw_c = NULL,
    tmle = average_est_cov_list(tmle_list), 
    aiptw = average_est_cov_list(aiptw_list),
    gcomp = average_est_cov_list(gcomp_list),
    QnMod = NULL, gnMod = NULL, QrnMod = NULL, grnMod = NULL,
    a_0 = a_0, call = call
  )
  # tack on nuisance if requested
  if(returnNuisance){
    # for compatibility with previous versions
    if(length(nuisance_drtmle_list) == 1){
      out$nuisance_drtmle <- nuisance_drtmle_list[[1]]
      out$nuisance_aiptw_c <- nuisance_aiptw_c_list[[1]] 
    }else{
      out$nuisance_drtmle <- nuisance_drtmle_list[[1]]
      out$nuisance_aiptw_c <- nuisance_aiptw_c_list[[1]] 
    }
    
  }

  # tack on models if requested
  if (returnModels) {
    if (!Qn_user) {
      if(n_loops == 1){
        out$QnMod <- QnMod_list[[1]]
      }else{
        out$QnMod <- QnMod_list
      }
    }
    if (!gn_user) {
      if(n_loops == 1){
        out$gnMod <- gnMod_list[[1]]
      }else{
        out$gnMod <- gnMod_list
      }
    }
    if ("Q" %in% guard) {
      if(n_loops == 1){
        out$QrnMod <- QrnMod_list[[1]]
      }else{
        out$QrnMod <- QrnMod_list
      }
    }
    if ("g" %in% guard) {
      if(n_loops == 1){
        out$grnMod <- grnMod_list[[1]]
      }else{
        out$grnMod <- grnMod_list
      }
    }
  }
  class(out) <- "drtmle"
  return(out)
}

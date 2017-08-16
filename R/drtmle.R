globalVariables(c("v", "%dopar%"))

#' TMLE estimate of the average treatment effect with doubly-robust inference
#' 
#' @param W A \code{data.frame} of named covariates.
#' @param A A vector of discrete-valued treatment assignment (assumed to be 
#' equal to 0 or 1).
#' @param Y numeric of continuous or binary outcomes.
#' @param DeltaY A vector of missing outcome indicator (assumed to be equal to 0 
#' if missing 1 if observed).
#' @param DeltaA A vector of missing treatment indicator (assumed to be equal to 
#' 0 if missing 1 if observed).
#' @param a_0 A vector of fixed treatment values at which to return marginal 
#' mean estimates.
#' @param family A \code{family} object equal to either \code{binomial()} or 
#' \code{gaussian()}, to be passed to the \code{SuperLearner} or \code{glm} 
#' function.
#' @param stratify A \code{boolean} indicating whether to estimate the outcome 
#' regression separately for different values of \code{A} (if \code{TRUE}) or to 
#' pool across \code{A} (if \code{FALSE}).
#' @param SL_Q A vector of characters or a list describing the Super Learner 
#' library to be used for the outcome regression. See 
#' \code{link{SuperLearner::SuperLearner}} for details.
#' @param SL_g A vector of characters describing the super learner library to be 
#' used for each of the propensity score regressions (\code{DeltaA}, \code{A}, 
#' and \code{DeltaY}). To use the same library for each of the regressions (or 
#' if there is no missing data in \code{A} nor \code{Y}), a single library may 
#' be input. See \code{link{SuperLearner::SuperLearner}} for details on how 
#' super learner libraries can be specified.
#' @param SL_Qr A vector of characters or a list describing the Super Learner 
#' library to be used for the reduced-dimension outcome regression. 
#' @param SL_gr A vector of characters or a list describing the Super Learner 
#' library to be used for the second reduced-dimension propensity score. 
#' @param glm_Q A character describing a formula to be used in the call to 
#' \code{glm} for the outcome regression. Ignored if \code{SL_Q!=NULL}.
#' @param glm_g A list of characters describing the formulas to be used
#' for each of the propensity score regressions (\code{DeltaA}, \code{A}, and 
#' \code{DeltaY}). To use the same formula for each of the regressions (or if 
#' there is no missing data in \code{A} nor \code{Y}), a single character 
#' formula may be input.
#' @param glm_Qr A character describing a formula to be used in the call to 
#' \code{glm} for reduced-dimension outcome regression. Ignored if 
#' \code{SL_Qr!=NULL}. The formula should use the variable name \code{'gn'}.
#' @param glm_gr A character describing a formula to be used in the call to 
#' \code{glm} for the reduced-dimension propensity score. Ignored if 
#' \code{SL_gr!=NULL}. The formula should use the variable name \code{'Qn'} and 
#' \code{'gn'} if \code{reduction='bivariate'} and \code{'Qn'} otherwise.
#' @param guard A character vector indicating what pattern of misspecifications 
#' to guard against. If \code{guard} contains \code{"Q"}, then the TMLE guards 
#' against misspecification of the outcome regression by estimating the 
#' reduced-dimension outcome regression specified by \code{glm_Qr} or 
#' \code{SL_Qr}. If \code{guard} contains \code{"g"} then the TMLE 
#' (additionally) guards against misspecification of the propensity score by 
#' estimating the reduced-dimension propensity score specified by \code{glm_gr} 
#' or \code{SL_gr}. 
#' @param reduction A character equal to \code{"univariate"} for a univariate 
#' misspecification correction (default) or \code{"bivariate"}
#' for the bivariate version. 
#' @param returnModels A boolean indicating whether to return model fits for the 
#' outcome regression, propensity score,
#' and reduced-dimension regressions.
#' @param maxIter A numeric that sets the maximum number of iterations the TMLE 
#' can perform in its fluctuation step.
#' @param tolIC A numeric that defines the stopping criteria based on the 
#' empirical mean of the influence function. 
#' @param tolg A numeric indicating the minimum value for estimates of the 
#' propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param Qsteps A numeric equal to 1 or 2 indicating whether the fluctuation 
#' submodel for the outcome regression
#' should be fit using a single minimization (\code{Qsteps = 1}) or a 
#' backfitting-type minimization (\code{Qsteps=2}). The latter was found to be 
#' more stable in simulations and is the default. 
#' @param cvFolds A numeric equal to the number of folds to be used in 
#' cross-validated fitting of nuisance parameters. If \code{cvFolds = 1}, no 
#' cross-validation is used.
#' @param parallel A boolean indicating whether to use \code{foreach}
#' to estimate nuisance parameters in parallel. Only useful if there is a 
#' registered parallel backend and \code{cvFolds > 1}.
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
#'        \code{tolIC} (unless \code{maxIter} was reached first).}
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
#' 
#' @importFrom plyr llply laply
#' @importFrom foreach foreach
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
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#' # fit drtmle with maxIter = 1 to run fast
#' fit1 <- drtmle(W = W, A = A, Y = Y, a_0 = c(1,0),
#'                family=binomial(),
#'                stratify=FALSE,
#'                SL_Q=c("SL.glm","SL.mean","SL.glm.interaction"),
#'                SL_g=c("SL.glm","SL.mean","SL.glm.interaction"),
#'                SL_Qr="SL.glm",
#'                SL_gr="SL.glm", maxIter = 1)

drtmle <- function(Y, A, W, 
                   DeltaA = as.numeric(!is.na(A)),
                   DeltaY = as.numeric(!is.na(Y)),
                   a_0 = unique(A[!is.na(A)]),
                   family = if(all(Y %in% c(0,1))){
                    stats::binomial()
                   }else{ stats::gaussian() },
                   stratify = TRUE,
                   SL_Q = NULL, SL_g = NULL,
                   SL_Qr = NULL, SL_gr = NULL,
                   glm_Q = NULL, glm_g = NULL,
                   glm_Qr = NULL, glm_gr = NULL,
                   guard = c("Q","g"),
                   reduction = "univariate",
                   returnModels = FALSE,
                   cvFolds = 1, 
                   maxIter = 3,
                   tolIC = 1/length(Y), 
                   tolg = 1e-2,
                   verbose = FALSE,
                   Qsteps = 2,
                   parallel = FALSE,
                   ...){
  call <- match.call()
  # if cvFolds non-null split data into cvFolds pieces
  n <- length(Y)
  if(cvFolds!=1){
    validRows <- split(sample(1:n), rep(1:cvFolds, length = n))       
  }else{
    validRows <- list(1:n)
  }
  #-------------------------------
  # estimate propensity score
  #-------------------------------
  if(!parallel){
    gnOut <- lapply(X = validRows, FUN = estimateG,
                    A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
                    tolg = tolg, verbose = verbose, stratify = stratify,
                    returnModels = returnModels, SL_g = SL_g,
                    glm_g = glm_g, a_0 = a_0)
  }else{
    gnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% 
    {
      estimateG(A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                tolg = tolg, verbose = verbose, stratify = stratify,
                returnModels = returnModels, SL_g = SL_g,
                glm_g = glm_g, a_0 = a_0, validRows = validRows[[v]])
    }
  }

  # # re-order predictions
  gnValid <- unlist(gnOut, recursive = FALSE, use.names = FALSE)
  gnUnOrd <- do.call(Map, c(c, gnValid[seq(1, length(gnValid), 2)]))
  gn <- vector(mode = "list", length = length(a_0))
  for(i in 1:length(a_0)){
    gn[[i]] <- rep(NA, n)
    gn[[i]][unlist(validRows)] <- gnUnOrd[[i]]
  }
  # obtain list of propensity score fits
  gnMod <- gnValid[seq(2, length(gnValid), 2)]

  #-------------------------------
  # estimate outcome regression
  #-------------------------------
  if(!parallel){
    QnOut <- lapply(X = validRows, FUN = estimateQ,
                    Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
                    verbose = verbose, returnModels = returnModels, 
                    SL_Q = SL_Q, a_0 = a_0, stratify = stratify, 
                    glm_Q = glm_Q, family = family)
  }else{
    QnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% 
    {
      estimateQ(Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
                verbose = verbose, returnModels = returnModels, 
                SL_Q = SL_Q, a_0 = a_0, glm_Q = glm_Q, family = family, 
                stratify = stratify, validRows = validRows[[v]])
    }
  }
  # re-order predictions
  QnValid <- unlist(QnOut, recursive = FALSE, use.names = FALSE)
  QnUnOrd <- do.call(Map, c(c, QnValid[seq(1,length(QnValid),2)]))
  Qn <- vector(mode = "list", length = length(a_0))
  for(i in 1:length(a_0)){
    Qn[[i]] <- rep(NA, n)
    Qn[[i]][unlist(validRows)] <- QnUnOrd[[i]]
  }
  # obtain list of outcome regression fits
  QnMod <- QnValid[seq(2,length(QnValid),2)]
  
  # naive g-computation estimate
  psi_n <- lapply(Qn, mean)
  
  # estimate influence function
  Dno <- eval_Dstar(A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA, 
                    Qn = Qn, gn = gn, psi_n = psi_n, a_0 = a_0)
  
  # estimate bias correction
  PnDn <- lapply(Dno, mean)
  
  # additional bias terms
  Dngo <- rep(0, n)
  DnQo <- rep(0, n)
  PnDQn <- PnDgn <- 0
  
  if("Q" %in% guard){
    if(!parallel){
      QrnOut <- lapply(X = validRows, FUN = estimateQrn, 
                       Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                       Qn = Qn, gn = gn, glm_Qr = glm_Qr, 
                       family = stats::gaussian(), SL_Qr = SL_Qr, 
                       a_0 = a_0, returnModels = returnModels)
    }else{
      QrnOut <- foreach::foreach(v = 1:cvFolds, 
                                 .packages = "SuperLearner") %dopar% {
        estimateQrn(Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                    Qn = Qn, gn = gn, glm_Qr = glm_Qr, 
                    family = stats::gaussian(), SL_Qr = SL_Qr, a_0 = a_0, 
                    returnModels = returnModels, validRows = validRows[[v]])
      }
    }
    # re-order predictions
    QrnValid <- unlist(QrnOut, recursive = FALSE, use.names = FALSE)
    QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid),2)]))
    Qrn <- vector(mode = "list", length = length(a_0))
    for(i in 1:length(a_0)){
      Qrn[[i]] <- rep(NA, n)
      Qrn[[i]][unlist(validRows)] <- QrnUnOrd[[i]]
    }
    # obtain list of reduced dimension regression fits
    QrnMod <- QrnValid[seq(2,length(QrnValid),2)]

    Dngo <- eval_Dstar_g(A = A, DeltaY = DeltaY, DeltaA = DeltaA, Qrn = Qrn, 
                         gn = gn, a_0 = a_0)
    PnDgn <- lapply(Dngo, mean)
  }
  if("g" %in% guard){
    if(!parallel){
      grnOut <- lapply(X = validRows, FUN = estimategrn,
                       Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                       tolg = tolg, Qn = Qn, gn = gn, 
                       glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0, 
                       reduction = reduction, returnModels = returnModels)      
    }else{
      grnOut <- foreach::foreach(v = 1:cvFolds, 
                                 .packages = "SuperLearner") %dopar% {
        estimategrn(Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                    tolg = tolg, Qn = Qn, gn = gn, glm_gr = glm_gr, 
                    SL_gr = SL_gr, reduction = reduction, a_0 = a_0,
                    returnModels = returnModels, validRows = validRows[[v]])
      }
    }
    # re-order predictions
    grnValid <- unlist(grnOut, recursive = FALSE, use.names = FALSE)
    grnUnOrd <- do.call(Map, c(rbind, grnValid[seq(1,length(grnValid),2)]))
    grn <- vector(mode = "list", length = length(a_0))
    for(i in 1:length(a_0)) {
      grn[[i]] <- data.frame(grn1 = rep(NA, n), grn2=rep(NA, n))
      grn[[i]][unlist(validRows),] <- cbind(grnUnOrd[[i]])
    }
    # obtain list of outcome regression fits
    grnMod <- grnValid[seq(2, length(grnValid), 2)]

    # evaluate extra piece of influence function
    DnQo <- eval_Dstar_Q(A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA, 
                         Qn = Qn, grn = grn, gn = gn, a_0 = a_0, 
                         reduction = reduction)
    PnDQn <- lapply(DnQo, mean)
  }
  
  # one step estimates
  psi_o1 <- mapply(a = psi_n, b = PnDn, SIMPLIFY = FALSE, 
                   FUN = function(a, b){ a + b })
  psi_o <- mapply(a = psi_n, b = PnDn, c = PnDQn, d = PnDgn, SIMPLIFY = FALSE,
                  FUN = function(a, b, c, d){ a + b - c - d }) 
  
  # covariance for one step
  Dno1Mat <- matrix(unlist(Dno), nrow = n, ncol = length(a_0))
  DnoMat <- matrix(unlist(Dno) - unlist(DnQo) - unlist(Dngo), 
                   nrow = n, ncol = length(a_0))
  
  cov_o1 <- stats::cov(Dno1Mat)/n 
  cov_o <- stats::cov(DnoMat)/n
  
  # initialize fluctuations
  QnStar <- Qn
  if("g" %in% guard) grnStar <- grn
  if("Q" %in% guard) QrnStar <- Qrn 
  gnStar <- gn; 
  PnDQnStar <- PnDgnStar <- PnDnoStar <- Inf
  eps <- list(Inf)
  ct <- 0

  # fluctuate
  while(max(abs(c(unlist(PnDQnStar), unlist(PnDgnStar), unlist(PnDnoStar)))) > 
        tolIC & ct < maxIter){
    ct <- ct + 1
    
    # re-estimate Qrn
    if("Q" %in% guard){
      # fluctuate gnStar
      gnStarOut <- fluctuateG(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                              DeltaY = DeltaY, a_0 = a_0, tolg = tolg, 
                              gn = gnStar, Qrn = QrnStar)
      gnStar <- plyr::llply(gnStarOut, function(x){ unlist(x$est) })
      epsg <- plyr::laply(gnStarOut, function(x){ x$eps })
    }else{
      epsg <- NA
    }
    
    # fluctuate QnStar
    if("g" %in% guard){
      if(!parallel){
        grnStarOut <- lapply(X = validRows, FUN = estimategrn,
                         Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                         tolg = tolg, Qn = QnStar, gn = gnStar, 
                         glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0, 
                         reduction = reduction, returnModels = returnModels)
      }else{
        grnStarOut <- foreach::foreach(v = 1:cvFolds, 
                                       .packages = "SuperLearner") %dopar% {
          estimategrn(Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                      tolg = tolg, Qn = QnStar, gn = gnStar, glm_gr = glm_gr, 
                      SL_gr = SL_gr, a_0=a_0, reduction = reduction, 
                      returnModels = returnModels, validRows = validRows[[v]])
        }
      }
      # re-order predictions
      grnValid <- unlist(grnStarOut, recursive = FALSE, use.names = FALSE)
      grnUnOrd <- do.call(Map, c(rbind, grnValid[seq(1,length(grnValid),2)]))
      grnStar <- vector(mode = "list", length = length(a_0))
      for(i in 1:length(a_0)){
        grnStar[[i]] <- data.frame(grn1 = rep(NA, n), grn2 = rep(NA, n))
        grnStar[[i]][unlist(validRows), ] <- cbind(grnUnOrd[[i]])
      }
      # obtain list of outcome regression fits
      grnMod <- grnValid[seq(2, length(grnValid), 2)]

      if(Qsteps==1){
        QnStarOut <- fluctuateQ(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                                DeltaY = DeltaY, a_0 = a_0, Qn = QnStar, 
                                gn = gnStar, grn = grnStar, 
                                reduction = reduction)
        QnStar <- plyr::llply(QnStarOut, function(x){ unlist(x$est) })
        epsQ <- plyr::laply(QnStarOut, function(x){ x$eps })
      }else if(Qsteps == 2){
        # do the extra targeting
        QnStarOut2 <- fluctuateQ2(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                                  DeltaY = DeltaY, a_0 = a_0, Qn = QnStar,   
                                  gn = gnStar, grn = grnStar, 
                                  reduction = reduction)
        QnStar <- plyr::llply(QnStarOut2, function(x){ unlist(x[[1]]) })
        
        # do the usual targeting
        QnStarOut1 <- fluctuateQ1(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                                  DeltaY = DeltaY, a_0 = a_0, Qn = QnStar, 
                                  gn = gnStar)
        QnStar <- plyr::llply(QnStarOut1, function(x){ unlist(x[[1]]) })
        
        # for later retrieval of fluct coefficients
        epsQ <- mapply(q1 = QnStarOut1, q2 = QnStarOut2, function(q1, q2){
          c(q1$eps,q2$eps)
        })
      }
    }else{
      QnStarOut <- fluctuateQ1(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                               DeltaY = DeltaY, a_0 = a_0, Qn = QnStar, 
                               gn = gnStar)
      QnStar <- plyr::llply(QnStarOut, function(x){ unlist(x[[1]]) })
      epsQ <- plyr::laply(QnStarOut, function(x){ x$eps })
    }
    
    if("Q" %in% guard){
      if(!parallel){
        QrnStarOut <- lapply(X = validRows, FUN = estimateQrn, 
                         Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                         Qn = QnStar, gn = gnStar, glm_Qr = glm_Qr, 
                         family = stats::gaussian(), SL_Qr = SL_Qr, a_0 = a_0,
                         returnModels = returnModels)
      }else{
        QrnStarOut <- foreach::foreach(v = 1:cvFolds, 
                                       .packages = "SuperLearner") %dopar% {
          estimateQrn(Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                      Qn = QnStar, gn = gnStar, glm_Qr = glm_Qr, 
                      family = stats::gaussian(), SL_Qr = SL_Qr, a_0 = a_0,
                      returnModels = returnModels, validRows = validRows[[v]])
        }
      }
      # re-order predictions
      QrnValid <- unlist(QrnStarOut, recursive = FALSE, use.names = FALSE)
      QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid), 2)]))
      QrnStar <- vector(mode = "list", length = length(a_0))
      for(i in 1:length(a_0)){
        QrnStar[[i]] <- rep(NA, n)
        QrnStar[[i]][unlist(validRows)] <- QrnUnOrd[[i]]
      }      
      # obtain list of reduced dimension regression fits
      QrnMod <- QrnValid[seq(2,length(QrnValid),2)]
    }
    
    # get fluctuation parameters
    eps <- c(epsQ, epsg)
    
    # tmle estimates
    psi_t <- lapply(QnStar, mean)
    
    # calculate influence functions
    DnoStar <- eval_Dstar(A = A, Y = Y, DeltaY = DeltaY, DeltaA = DeltaA, 
                    Qn = QnStar, gn = gnStar, psi_n = psi_t, a_0 = a_0)
    PnDnoStar <- lapply(DnoStar, mean)
    
    if("g" %in% guard){
      DnQoStar <- eval_Dstar_Q(A = A, Y = Y, DeltaY = DeltaY, 
                         DeltaA = DeltaA, Qn = QnStar, grn = grnStar, gn = gn,
                         a_0 = a_0, reduction = reduction)
      PnDQnStar <- lapply(DnQoStar, mean)
    }
    if("Q" %in% guard){
      DngoStar <- eval_Dstar_g(A = A, DeltaY = DeltaY, DeltaA = DeltaA, 
                               Qrn = QrnStar, gn = gnStar, a_0 = a_0)
      PnDgnStar <- lapply(DngoStar, mean)
    }
    if(verbose){
      cat("TMLE Iteration", ct, "=", round(unlist(eps), 5), "\n")
      cat("Mean of IC       =", round(c(unlist(PnDnoStar), 
                                        unlist(PnDQnStar), 
                                        unlist(PnDgnStar)), 10),"\n")
    }
  }
  
  # standard tmle fluctuations
  QnStar1Out <- fluctuateQ1(Y = Y, A = A, W = W, DeltaA = DeltaA, 
                            DeltaY = DeltaY, Qn = Qn, gn = gn, a_0 = a_0)
  QnStar1 <- plyr::llply(QnStar1Out, function(x){ unlist(x[[1]]) })
  
  # tmle estimates
  psi_t <- lapply(QnStar, mean)
  psi_t1 <- lapply(QnStar1, mean)
  
  # covariance for tmle
  Dno1Star <- eval_Dstar(A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY, 
                         Qn = QnStar1, gn = gn, psi_n = psi_t1, a_0 = a_0)
  Dno1StarMat <- matrix(unlist(Dno1Star), nrow = n, ncol = length(a_0))
  cov_t1 <- stats::cov(Dno1StarMat)

  # covariance for drtmle
  DnoStar <- eval_Dstar(A = A, Y = Y, DeltaA = DeltaA, DeltaY = DeltaY, 
                        Qn = QnStar, gn = gnStar, psi_n = psi_t, a_0 = a_0)
  PnDnoStar <- lapply(DnoStar, mean)
  
  DnQoStar <- rep(list(rep(0, n)), length(a_0))
  DngoStar <- rep(list(rep(0, n)), length(a_0))
  
  if("g" %in% guard){
    DnQoStar <- eval_Dstar_Q(A = A, Y = Y, DeltaY = DeltaY, 
                       DeltaA = DeltaA, Qn = QnStar, grn = grnStar, gn = gn, 
                       a_0 = a_0, reduction = reduction)
    PnDQnStar <- lapply(DnQoStar, mean)
  }
  if("Q" %in% guard){
    DngoStar <- eval_Dstar_g(A = A, DeltaY = DeltaY, DeltaA = DeltaA, 
                             Qrn = QrnStar, gn = gnStar, a_0 = a_0)
    PnDgnStar <- lapply(DngoStar, mean)
  }
  
  DnoStarMat <- matrix(unlist(DnoStar) - unlist(DnQoStar) - unlist(DngoStar), 
                       nrow = n, ncol = length(a_0))
  cov_t <- stats::cov(DnoStarMat)/n 

  out <- list(drtmle = list(est = unlist(psi_t), cov = cov_t),
              nuisance_drtmle = list(QnStar = QnStar, gnStar = gnStar,
                                     QrnStar = QrnStar, grnStar = grn,
                                     meanIC = unlist(c(PnDnoStar, PnDQnStar, 
                                                       PnDgnStar))),
              ic_drtmle = list(eif = PnDnoStar, missQ = PnDgnStar, 
                               missg = PnDQnStar),
              aiptw_c = list(est = unlist(psi_o),cov=cov_o),
              nuisance_aiptw_c = list(Qn = Qn, gn = gn, Qrn = Qrn, grn = grn),
              tmle = list(est = unlist(psi_t1),cov = cov_t1),
              aiptw = list(est = unlist(psi_o1), cov = cov_o1),
              gcomp=list(est = unlist(psi_n), cov = cov_o1),
              QnMod = NULL, gnMod = NULL, QrnMod = NULL, grnMod = NULL,
              a_0 = a_0, call = call)

  # tack on models if requested
  if(returnModels){
    out$QnMod <- QnMod
    out$gnMod <- gnMod
    out$QrnMod <- QrnMod
    out$grnMod <- grnMod
  }
  class(out) <- "drtmle"
  return(out)
}



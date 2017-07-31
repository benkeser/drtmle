#' TMLE estimate of the average treatment effect with doubly-robust inference
#' 
#' @param W A \code{data.frame} of named covariates
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes. 
#' @param a_0 A vector of treatment levels at which to compute the adjusted mean outcome. 
#' @param family A \code{family} object equal to either \code{binomial()} or \code{gaussian()}, 
#' to be passed to the \code{SuperLearner} or \code{glm} function.
#' @param stratify A \code{boolean} indicating whether to estimate the outcome regression separately
#' for observations with \code{A} equal to 0/1 (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param SL_Q A vector of characters or a list describing the Super Learner library to be used 
#' for the outcome regression. See \code{link{SuperLearner::SuperLearner}} for details.
#' @param SL_g A vector of characters or a list describing the Super Learner library to be used 
#' for the propensity score. See \code{link{SuperLearner::SuperLearner}} for details.
#' @param SL_Qr A vector of characters or a list describing the Super Learner library to be used 
#' for the first reduced-dimension regression. 
#' @param SL_gr A vector of characters or a list describing the Super Learner library to be used 
#' for the second reduced-dimension regression.
#' @param glm_Q A character describing a formula to be used in the call to \code{glm} for the outcome regression. Ignored
#' if \code{SL_Q!=NULL}.
#' @param glm_g A character describing a formula to be used in the call to \code{glm} for the propensity score. Ignored
#' if \code{SL_g!=NULL}.
#' @param glm_Qr A character describing a formula to be used in the call to \code{glm} for the first reduced-dimension regression. Ignored
#' if \code{SL_Qr!=NULL}. The formula should use the variable name \code{'gn'}.
#' @param glm_gr A character describing a formula to be used in the call to \code{glm} for the second reduced-dimension regression. Ignored
#' if \code{SL_gr!=NULL}. The formula should use the variable name \code{'Qn'} and \code{'gn'} if 
#' \code{reduction='bivariate'} and \code{'Qn'} otherwise.
#' @param guard A character vector indicating what pattern of misspecifications to guard against. If \code{guard} contains \code{"Q"}, 
#' then the TMLE guards against misspecification of the outcome regression by estimating the reduced dimension regression
#' specified by \code{glm_Qr} or \code{SL_Qr}. If \code{guard} contains \code{"g"} then the TMLE 
#' (additionally) guards against misspecification of the propensity score by estimating the reduced dimension regression
#' specified by \code{glm_gr} or \code{SL_gr}. If \code{NULL}, the usual TMLE is computed.
#' @param reduction A character equal to \code{"univariate"} for a univariate misspecification correction or \code{"bivariate"}
#' for the bivariate version. 
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' @param maxIter A numeric that sets the maximum number of iterations the TMLE can perform in its fluctuation step.
#' @param tolIC A numeric that defines the stopping criteria based on the empirical mean
#' of the scores of the fluctuation submodels submodels. Setting to \code{"default"}
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param Qsteps A numeric equal to 1/2 indicating whether the fluctuation submodel for the outcome regression
#' should be fit using a single minimization (\code{Qsteps = 1}) or a backfitting-type minimization (\code{Qsteps=2}). 
#' The latter was found to be more stable in simulations. 
#' @param cvFolds A numeric equal to the number of folds to be used in cross-validated fitting of 
#' nuisance parameters. If \code{NULL}, no cross-validation is used. 
#' @param ... Other options (not currently used)
#' 
#' @return An object of class \code{"drtmle"}.
#' \describe{
#'  \item{\code{drtmle}}{A \code{list} of doubly-robust point estimates and 
#'        a doubly-robust covariance matrix}
#'  \item{\code{nuisance_drtmle}}{A \code{list} of the final TMLE estimates of the
#'        outcome regression (\code{$QnStar}), propensity score (\code{$gnStar}), 
#'        and reduced-dimension regressions (\code{$QrnStar}, \code{$grnStar}) evaluated
#'        at the observed data values.}
#'  \item{\code{ic_drtmle}}{A \code{list} of the empirical mean of the efficient 
#'        influence function (\code{$eif}) and the extra pieces of the influence
#'        function resulting from misspecification. All should be smaller than 
#'        \code{tolIC} (unless \code{maxIter} was reached first).}
#'  \item{\code{aiptw_c}}{A \code{list} of doubly-robust point estimates and 
#'        a non-doubly-robust covariance matrix. Theory does not guarantee performance of 
#'        inference for these estimators, but simulation studies showed they often
#'        perform adequately.}
#'  \item{\code{nuisance_aiptw}}{A \code{list} of the initial estimates of the
#'        outcome regression, propensity score, and reduced-dimension regressions 
#'        evaluated at the observed data values.}
#'  \item{\code{tmle}}{A \code{list} of doubly-robust point estimates and non-doubly-robust
#'        covariance for the standard TMLE estimator.}
#'  \item{\code{aiptw}}{A \code{list} of doubly-robust point estimates and non-doubly-robust
#'        covariance matrix for the standard AIPTW estimator.}
#'  \item{\code{gcomp}}{A \code{list} of non-doubly-robust point estimates and non-doubly-robust
#'        covariance matrix for the standard G-computation estimator. If super learner is used
#'        there is no guarantee of correct inference for this estimator.}
#'  \item{\code{QnMod}}{The fitted object for the outcome regression. Returns \code{NULL}
#'        if \code{returnModels = FALSE}.}
#'  \item{\code{gnMod}}{The fitted object for the propensity score. Returns \code{NULL} if
#'        \code{returnModels = FALSE}.}
#'  \item{\code{QrnMod}}{The fitted object for the reduced-dimension regression that guards 
#'        against misspecification of the outcome regression. Returns \code{NULL} if 
#'        \code{returnModels = FALSE}.}
#'  \item{\code{grnMod}}{The fitted object for the reduced-dimension regression that guards
#'        against misspecification of the propensity score. Returns \code{NULL} if 
#'        \code{returnModels = FALSE}.}
#'  \item{\code{a_0}}{The treatment levels that were requested for computation of 
#'        covariate-adjusted means.}
#' }
#' 
#' 
#' @importFrom plyr llply laply
#' 
#' 
#' @export 
#' 
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # simulate data
#' set.seed(123456)
#' n <- 200
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#' # fit drtmle
#' fit1 <- drtmle(W = W, A = A, Y = Y, a_0 = c(1,0),
#'                family=binomial(),
#'                stratify=FALSE,
#'                SL_Q=c("SL.glm","SL.mean","SL.step"),
#'                SL_g=c("SL.glm","SL.mean","SL.step"),
#'                SL_Qr="SL.npreg",
#'                SL_gr="SL.npreg")

drtmle <- function(Y, A, W, a_0 = unique(A),
                   family = stats::binomial(),
                   stratify = TRUE,
                   SL_Q = NULL,
                   SL_g = NULL,
                   SL_Qr = NULL,
                   SL_gr = NULL,
                   glm_Q = NULL,
                   glm_g = NULL,
                   glm_Qr = NULL,
                   glm_gr = NULL,
                   guard = c("Q","g"),
                   reduction = "univariate",
                   returnModels = FALSE,
                   cvFolds = NULL, 
                   maxIter = 3,
                   tolIC = 1/length(Y), 
                   tolg = 1e-2,
                   verbose = FALSE,
                   Qsteps = 2,
                   ...){
  # if cvFolds non-null split data into cvFolds pieces
  n <- length(Y)
  if(!is.null(cvFolds)){
    validRows <- split(sample(1:n), rep(1:cvFolds, length = n))       
    ordVR <- order(unlist(validRows, use.names = FALSE))
  }else{
    validRows <- list(NULL)
    ordVR <- 1:n
  }
  #-------------------------------
  # estimate propensity score
  #-------------------------------
  gnOut <- lapply(X = validRows, FUN = estimateG,
                  A=A, W=W, tolg=tolg, verbose=verbose, 
                  returnModels=returnModels,SL_g=SL_g,
                  glm_g=glm_g, a_0=a_0)
  # re-order predictions
  gnValid <- unlist(gnOut, recursive = FALSE, use.names = FALSE)
  gnUnOrd <- do.call(Map, c(c, gnValid[seq(1,length(gnValid),2)]))
  gn <- lapply(gnUnOrd, function(x){ x[ordVR] })
  # obtain list of propensity score fits
  gnMod <- gnValid[seq(2,length(gnValid),2)]
  # TO DO: Add reasonable names to gnMod?

  #-------------------------------
  # estimate outcome regression
  #-------------------------------
  QnOut <- lapply(X = validRows, FUN = estimateQ,
                  Y=Y, A=A, W=W, verbose=verbose, 
                  returnModels=returnModels, SL_Q=SL_Q, a_0=a_0,
                  stratify=stratify, glm_Q=glm_Q,family=family)
  # re-order predictions
  QnValid <- unlist(QnOut, recursive = FALSE, use.names = FALSE)
  QnUnOrd <- do.call(Map, c(c, QnValid[seq(1,length(QnValid),2)]))
  Qn <- lapply(QnUnOrd, function(x){ x[ordVR] })
  # obtain list of outcome regression fits
  QnMod <- QnValid[seq(2,length(QnValid),2)]
  # TO DO: Add reasonable names to QnMod?
  

  # naive g-computation estimate
  psi.n <- lapply(Qn, mean)
  
  # estimate influence function
  Dno <- mapply(a=split(a_0,1:length(a_0)),Q=Qn,g=gn,p=psi.n,FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  
  # estimate bias correction
  PnDn <- lapply(Dno, mean)
  
  # additional bias terms
  PnDQn <- PnDgn <- 0
  Dngo <- rep(0, length(Y))
  DnQo <- rep(0, length(Y))
  
  if("Q" %in% guard){
      QrnOut <- lapply(X = validRows, FUN = estimateQrn, 
                       Y=Y, A=A, W=W, Qn=Qn, gn=gn, glm_Qr=glm_Qr, 
                       SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels)
      # re-order predictions
      QrnValid <- unlist(QrnOut, recursive = FALSE, use.names = FALSE)
      QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid),2)]))
      Qrn <- lapply(QrnUnOrd, function(x){ x[ordVR] })
      # obtain list of reduced dimension regression fits
      QrnMod <- QrnValid[seq(2,length(QrnValid),2)]
      # TO DO: Add reasonable names to QrnMod?

      Dngo <- mapply(a=split(a_0,1:length(a_0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
        Qr/g * (as.numeric(A==a) - g)
      },SIMPLIFY=FALSE)
      PnDgn <- lapply(Dngo, mean)
  }
  if("g" %in% guard){
      grnOut <- lapply(X=validRows, FUN = estimategrn,
                       Y=Y,A=A, W=W, tolg=tolg, Qn=Qn, gn=gn, 
                       glm_gr=glm_gr, SL_gr=SL_gr, a_0=a_0, 
                       reduction=reduction,returnModels = returnModels)      
      # re-order predictions
      grnValid <- unlist(grnOut, recursive = FALSE, use.names = FALSE)
      grnUnOrd <- do.call(Map, c(rbind, grnValid[seq(1,length(grnValid),2)]))
      grn <- lapply(grnUnOrd, function(x){ x[ordVR,] })
      # obtain list of outcome regression fits
      grnMod <- grnValid[seq(2,length(grnValid),2)]
      # TO DO: Add reasonable names to grnMod?

      if(reduction=="univariate"){
        DnQo <- mapply(a=split(a_0,1:length(a_0)),Q=Qn,gr=grn,FUN=function(a,Q,gr){
          as.numeric(A==a)/gr$grn2 * gr$grn1 * (Y-Q)
        },SIMPLIFY=FALSE)
      }else if(reduction=="bivariate"){
        DnQo <- mapply(a=split(a_0,1:length(a_0)),Q=Qn,g=gn, gr=grn,FUN=function(a,Q,gr,g){
          as.numeric(A==a)/gr$grn2 * (gr$grn2 - g)/g * (Y-Q)
        },SIMPLIFY=FALSE)
      }
      PnDQn <- lapply(DnQo, mean)
  }
  
  # one step estimates
  psi.o1 <- mapply(a=psi.n,b=PnDn,SIMPLIFY=FALSE, 
                   FUN=function(a,b){a+b})
  psi.o <- mapply(a=psi.n,b=PnDn,c=PnDQn,d=PnDgn,SIMPLIFY=FALSE,
                  FUN=function(a,b,c,d){a+b-c-d}) 
  
  # covariance for one step
  Dno1Mat <- matrix(unlist(Dno),ncol=length(Y), nrow=length(a_0), byrow=TRUE)
  DnoMat <- matrix(unlist(Dno) - unlist(DnQo) - unlist(Dngo), ncol=length(Y), nrow=length(a_0), byrow=TRUE)
  
  cov.o1 <- (Dno1Mat-mean(Dno1Mat))%*%t(Dno1Mat-mean(Dno1Mat)) /(length(Y)^2)
  cov.o <- (DnoMat-mean(DnoMat))%*%t(DnoMat-mean(DnoMat))/(length(Y)^2)
  
  # initialize fluctuations
  QnStar <- Qn
  if("g" %in% guard) grnStar <- grn
  if("Q" %in% guard) QrnStar <- Qrn 
  gnStar <- gn; 
  PnDQnStar <- PnDgnStar <- PnDnoStar <- Inf
  eps <- list(Inf)
  ct <- 0

  # fluctuate
  while(max(abs(c(unlist(PnDQnStar),unlist(PnDgnStar),unlist(PnDnoStar)))) > tolIC & ct < maxIter){
    ct <- ct + 1
    
    # re-estimate Qrn
    if("Q" %in% guard){
      # fluctuate gnStar
      gnStarOut <- fluctuateG(Y=Y, A=A, W=W, a_0=a_0, tolg=tolg, gn=gnStar, Qrn=QrnStar)
      gnStar <- plyr::llply(gnStarOut, function(x){unlist(x$est)})
      epsg <- plyr::laply(gnStarOut, function(x){x$eps})
    }else{
      epsg <- NA
    }
    
    # fluctuate QnStar
    if("g" %in% guard){
      grnStarOut <- lapply(X=validRows, FUN = estimategrn,
                       Y=Y,A=A, W=W, tolg=tolg, Qn=QnStar, gn=gnStar, 
                       glm_gr=glm_gr, SL_gr=SL_gr, a_0=a_0, 
                       reduction=reduction,returnModels = returnModels)
      # re-order predictions
      grnValid <- unlist(grnStarOut, recursive = FALSE, use.names = FALSE)
      grnUnOrd <- do.call(Map, c(rbind, grnValid[seq(1,length(grnValid),2)]))
      grn <- lapply(grnUnOrd, function(x){ x[ordVR,] })
      # obtain list of outcome regression fits
      grnMod <- grnValid[seq(2,length(grnValid),2)]
      # TO DO: Add reasonable names to grnMod?

      if(Qsteps==1){
        QnStarOut <- fluctuateQ(Y=Y, A=A, W=W, a_0=a_0, Qn=QnStar, 
                                gn=gnStar, grn=grnStar, reduction=reduction)
        QnStar <- plyr::llply(QnStarOut, function(x){unlist(x$est)})
        epsQ <- plyr::laply(QnStarOut, function(x){x$eps})
      }else if(Qsteps==2){
        # do the extra targeting
        QnStarOut2 <- fluctuateQ2(Y=Y, A=A, W=W, a_0=a_0, Qn=QnStar, 
                                  gn=gnStar, grn=grnStar, reduction=reduction)
        QnStar <- plyr::llply(QnStarOut2, function(x){unlist(x[[1]])})
        
        # do the usual targeting
        QnStarOut1 <- fluctuateQ1(Y=Y, A=A, W=W, a_0=a_0, Qn=QnStar, gn=gnStar)
        QnStar <- plyr::llply(QnStarOut1, function(x){unlist(x[[1]])})
        
        # for later retrieval of fluct coefficients
        epsQ <- mapply(q1=QnStarOut1, q2=QnStarOut2, function(q1,q2){c(q1$eps,q2$eps)})
      }
    }else{
      QnStarOut <- fluctuateQ1(Y=Y, A=A, W=W, a_0=a_0, Qn=QnStar, gn=gnStar)
      QnStar <- plyr::llply(QnStarOut, function(x){unlist(x[[1]])})
      epsQ <- plyr::laply(QnStarOut, function(x){x$eps})
    }
    
    if("Q" %in% guard){
      QrnStarOut <- lapply(X = validRows, FUN = estimateQrn, 
                       Y=Y, A=A, W=W, Qn=QnStar, gn=gnStar, glm_Qr=glm_Qr, 
                       SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels)
      # re-order predictions
      QrnValid <- unlist(QrnStarOut, recursive = FALSE, use.names = FALSE)
      QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid),2)]))
      QrnStar <- lapply(QrnUnOrd, function(x){ x[ordVR] })
      # obtain list of reduced dimension regression fits
      QrnMod <- QrnValid[seq(2,length(QrnValid),2)]
      # TO DO: Add reasonable names to QrnMod?
    }
    
    # get fluctuation parameters
    eps <- c(epsQ, epsg)
    
    # tmle estimates
    psi.t <- lapply(QnStar, mean)
    
    # calculate influence functions
    DnoStar <- mapply(a=split(a_0,1:length(a_0)),
                      Q=QnStar,g=gnStar,p=psi.t,
                      FUN=function(a,Q,g,p){
      as.numeric(A==a)/g * (Y - Q) + Q - p
    },SIMPLIFY=FALSE)
    PnDnoStar <- lapply(DnoStar, mean)
    
    if("g" %in% guard){
      if(reduction=="univariate"){
        DnQoStar <- mapply(a=split(a_0,1:length(a_0)),
                           Q=QnStar,gr=grnStar,
                           FUN=function(a,Q,gr){
          as.numeric(A==a)/gr$grn2 * gr$grn1 * (Y-Q)
        }, SIMPLIFY=FALSE)
      }else if(reduction=="bivariate"){
        DnQoStar <- mapply(a=split(a_0,1:length(a_0)),
                           Q=QnStar,g=gnStar,gr=grnStar,
                           FUN=function(a,Q,gr,g){
          as.numeric(A==a)/gr$grn2 * (gr$grn2 - g)/g * (Y-Q)
        },SIMPLIFY=FALSE)
      }
      PnDQnStar <- lapply(DnQoStar, mean)
    }
    if("Q" %in% guard){
      DngoStar <- mapply(a=split(a_0,1:length(a_0)),
                         Qr=QrnStar,
                         g=gnStar,FUN=function(a,Qr,g){
        Qr/g * (as.numeric(A==a) - g)
      }, SIMPLIFY=FALSE)
      PnDgnStar <- lapply(DngoStar, mean)
    }
    if(verbose){
      cat("TMLE Iteration", ct, "=", round(unlist(eps),5), "\n")
      cat("Mean of IC       =", round(c(unlist(PnDnoStar), unlist(PnDQnStar), unlist(PnDgnStar)), 10),"\n")
    }
  }
  
  # standard tmle fluctuations
  QnStar1Out <- fluctuateQ1(Y=Y,A=A,W=W,Qn=Qn,gn=gn,a_0=a_0)
  QnStar1 <- plyr::llply(QnStar1Out, function(x){unlist(x[[1]])})
  
  # tmle estimates
  psi.t <- lapply(QnStar, mean)
  psi.t1 <- lapply(QnStar1, mean)
  
  # covariance for tmle
  Dno1Star <- mapply(a=split(a_0,1:length(a_0)),
                     Q=QnStar1,g=gn,p=psi.n,
                     FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  
  DnoStar <- mapply(a=split(a_0,1:length(a_0)),
                    Q=QnStar,g=gnStar,p=psi.n,
                    FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  PnDnoStar <- lapply(DnoStar, mean)
  
  DnQoStar <- list(rep(0, length(Y)))
  DngoStar <- list(rep(0, length(Y)))
  
  if("g" %in% guard){
    if(reduction=="univariate"){
      DnQoStar <- mapply(a=split(a_0,1:length(a_0)),
                         Q=QnStar,gr=grnStar,
                         FUN=function(a,Q,gr){
        as.numeric(A==a)/gr$grn2 * gr$grn1 * (Y-Q)
      }, SIMPLIFY=FALSE)
    }else if(reduction=="bivariate"){
      DnQoStar <- mapply(a=split(a_0,1:length(a_0)),Q=QnStar,
                         g=gn,gr=grnStar,FUN=function(a,Q,gr,g){
        as.numeric(A==a)/gr$grn2 * (gr$grn2-g)/g * (Y-Q)
      }, SIMPLIFY=FALSE)
    }
    PnDQnStar <- lapply(DnQoStar, mean)
  }
  
  if("Q" %in% guard){
    DngoStar <- mapply(a=split(a_0,1:length(a_0)),
                       Qr=QrnStar,g=gnStar,
                       FUN=function(a,Qr,g){
      Qr/g * (as.numeric(A==a) - g)
    }, SIMPLIFY=FALSE)
    PnDgnStar <- lapply(DngoStar, mean)
  }
  
  Dno1StarMat <- matrix(unlist(Dno1Star), ncol=length(Y), 
                        nrow=length(a_0), byrow=TRUE)
  DnoStarMat <- matrix(unlist(DnoStar) - unlist(DnQoStar) - unlist(DngoStar), 
                       ncol=length(Y), nrow=length(a_0),byrow=TRUE)
  
  cov.t1 <- (Dno1StarMat - mean(Dno1StarMat))%*%
              t((Dno1StarMat - mean(Dno1StarMat)))/(length(Y)^2)
  cov.t <- (DnoStarMat - mean(DnoStarMat))%*%
              t((DnoStarMat - mean(DnoStarMat)))/(length(Y)^2)
  

  out <- list(drtmle = list(est = unlist(psi.t), cov = cov.t),
              nuisance_drtmle = list(QnStar = QnStar, gnStar = gnStar,
                                     QrnStar = QrnStar, grnStar = grn),
              ic_drtmle = list(eif = PnDnoStar, missQ = PnDgnStar, missg = PnDQnStar),
              aiptw_c = list(est = unlist(psi.o),cov=cov.o),
              nuisance_aiptw = list(Qn = Qn, gn = gn, Qrn = Qrn, grn = grn),
              tmle = list(est=unlist(psi.t1),cov=cov.t1),
              aiptw = list(est=unlist(psi.o1), cov=cov.o1),
              gcomp=list(est=unlist(psi.n), cov=cov.o1),
              QnMod = NULL, gnMod = NULL, QrnMod = NULL, grnMod = NULL,
              a_0 = a_0)

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



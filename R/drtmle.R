#' drtmle 
#' 
#' Main function to estimate TMLE with doubly-robust inference
#' 
#' @param Y A vector of continuous or binary outcomes. 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param a0 A vector of treatment levels at which to compute the adjusted mean outcome. 
#' @param family A \code{character} equal to either \code{"binomial"} or \code{"gaussian"}. This is 
#' passed to the \code{SuperLearner} function internally.
#' @param stratify A \code{boolean} indicating whether to estimate the outcome regression separately
#' for observations with \code{A} equal to 0/1 (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param libraryQ A vector of characters or a list describing the Super Learner library to be used 
#' for the outcome regression. See \code{link{SuperLearner::SuperLearner}} for details.
#' @param libraryg A vector of characters or a list describing the Super Learner library to be used 
#' for the propensity score. See \code{link{SuperLearner::SuperLearner}} for details.
#' @param libraryQr A vector of characters or a list describing the Super Learner library to be used 
#' for the first reduced-dimension regression.
#' @param librarygr A vector of characters or a list describing the Super Learner library to be used 
#' for the second reduced-dimension regression.
#' @param glmQ A character describing a formula to be used in the call to \code{glm} for the outcome regression. Ignored
#' if \code{libraryQ!=NULL}.
#' @param glmg A character describing a formula to be used in the call to \code{glm} for the propensity score. Ignored
#' if \code{libraryg!=NULL}.
#' @param glmQr A character describing a formula to be used in the call to \code{glm} for the first reduced-dimension regression. Ignored
#' if \code{libraryQr!=NULL}.
#' @param glmgr A character describing a formula to be used in the call to \code{glm} for the second reduced-dimension regression. Ignored
#' if \code{librarygr!=NULL}.
#' @param guard A character vector indicating what pattern of misspecifications to guard against. If \code{guard} contains \code{"Q"}, 
#' then the TMLE guards against misspecification of the outcome regression. If \code{guard} contains \code{"g"} then the TMLE 
#' (additionally) guards against misspecification of the propensity score. If \code{NULL}, the usual TMLE is computed.
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
#' @param ... Other options (not currently used)
#' 
#' @importFrom plyr llply laply
#' 
#' 
#' @export 


drtmle <- function(Y,A,W,
                      a0=unique(A),
                      family="binomial",
                      stratify=TRUE,
                      libraryQ=NULL,
                      libraryg=NULL,
                      libraryQr=NULL,
                      librarygr=NULL,
                      glmQ=NULL,
                      glmg=NULL,
                      glmQr=NULL,
                      glmgr=NULL,
                      guard=c("Q","g"),
                      reduction="univariate",
                      returnModels=FALSE,
                      maxIter=100,
                      tolIC=1/length(Y), 
                      tolg=1e-8,
                      verbose=TRUE,
                      Qsteps=2,
                      ...){

  if(tolIC=="default") tolIC <- sqrt(1/length(Y))/50

  # estimate g
  gnOut <- estimateG(A=A, W=W, tolg=tolg, verbose=verbose, returnModels=returnModels,libraryg=libraryg, 
    glmg=glmg, a0=a0)
  if(returnModels){
    gnEst <- gnOut[[1]]
    gnMod <- gnOut[[2]]
    gn <- gnEst
  }else{
    gn <- gnOut
  }
  
  # estimate Q
  QnOut <- estimateQ(Y=Y, A=A, W=W, verbose=verbose, returnModels=returnModels, libraryQ=libraryQ, a0=a0,
                  stratify=stratify, glmQ=glmQ,family=family)
  if(returnModels){
    QnEst <- QnOut[[1]]
    QnMod <- QnOut[[2]]
    Qn <- QnEst
  }else{
    Qn <- QnOut
  }
  
  # naive estimate
  psi.n <- lapply(Qn, mean)
  
  # estimate influence function
  Dno <- mapply(a=split(a0,1:length(a0)),Q=Qn,g=gn,p=psi.n,FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  
  # estimate bias correction
  PnDn <- lapply(Dno, mean)
  
  # additional bias terms
  PnDQn <- PnDgn <- 0
  Dngo <- rep(0, length(Y))
  DnQo <- rep(0, length(Y))
  
  if("Q" %in% guard){
      Qrn <- estimateQrn(Y=Y, A=A, W=W, Qn=Qn, gn=gn, glmQr=glmQr, libraryQr=libraryQr, a0=a0)
      Dngo <- mapply(a=split(a0,1:length(a0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
        Qr/g * (as.numeric(A==a) - g)
      },SIMPLIFY=FALSE)
      PnDgn <- lapply(Dngo, mean)
  }
  if("g" %in% guard){
      grn <- estimategrn(Y=Y,A=A, W=W, tolg=tolg, Qn=Qn, gn=gn, glmgr=glmgr, librarygr=librarygr, a0=a0, reduction=reduction)
      if(reduction=="univariate"){
        DnQo <- mapply(a=split(a0,1:length(a0)),Q=Qn,gr=grn,FUN=function(a,Q,gr){
          as.numeric(A==a)/gr[[2]] * gr[[1]] * (Y-Q)
        },SIMPLIFY=FALSE)
      }else if(reduction=="bivariate"){
        DnQo <- mapply(a=split(a0,1:length(a0)),Q=Qn,g=gn, gr=grn,FUN=function(a,Q,gr,g){
          as.numeric(A==a)/gr[[1]] * (gr[[1]] - g)/g * (Y-Q)
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
  Dno1Mat <- matrix(unlist(Dno),ncol=length(Y), nrow=length(a0), byrow=TRUE)
  DnoMat <- matrix(unlist(Dno) - unlist(DnQo) - unlist(Dngo), ncol=length(Y), nrow=length(a0), byrow=TRUE)
  
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
      gnStarOut <- fluctuateG(Y=Y, A=A, W=W, a0=a0, tolg=tolg, Qn=QnStar, gn=gnStar, Qrn=QrnStar)
      gnStar <- llply(gnStarOut, function(x){unlist(x[[1]])})
      epsg <- laply(gnStarOut, function(x){x$eps})
    }else{
      epsg <- NA
    }
    
    # fluctuate QnStar
    if("g" %in% guard){
      grnStar <- estimategrn(Y=Y, A=A, W=W, reduction=reduction, tolg=tolg, a0=a0, Qn=QnStar, gn=gnStar, glmgr=glmgr, librarygr=librarygr)
      
      if(Qsteps==1){
        QnStarOut <- fluctuateQ(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, grn=grnStar, reduction=reduction)
        QnStar <- llply(QnStarOut, function(x){unlist(x[[1]])})
        epsQ <- laply(QnStarOut, function(x){x$eps})
      }else if(Qsteps==2){
        # do the extra targeting
        QnStarOut2 <- fluctuateQ2(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, grn=grnStar, reduction=reduction)
        QnStar <- llply(QnStarOut2, function(x){unlist(x[[1]])})
        
        # update grn based on new Qn
        # grnStar <- estimategrn(Y=Y, A=A, W=W, reduction=reduction, tolg=tolg, a0=a0, Qn=QnStar, gn=gnStar, glmgr=glmgr, librarygr=librarygr)
        
        # do the usual targeting
        QnStarOut1 <- fluctuateQ1(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar)
        QnStar <- llply(QnStarOut1, function(x){unlist(x[[1]])})
        
        # for later retrieval of fluct coefficients
        epsQ <- mapply(q1=QnStarOut1, q2=QnStarOut2, function(q1,q2){c(q1$eps,q2$eps)})
      }
    }else{
      QnStarOut <- fluctuateQ1(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar)
      QnStar <- llply(QnStarOut, function(x){unlist(x[[1]])})
      epsQ <- laply(QnStarOut, function(x){x$eps})
    }
    
    if("Q" %in% guard){
      QrnStar <- estimateQrn(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, glmQr=glmQr, libraryQr=libraryQr)
    }
    
    # get fluctuation parameters
    eps <- c(epsQ, epsg)
    
    # tmle estimates
    psi.t <- lapply(QnStar, mean)
    
    # calculate influence functions
    DnoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,g=gnStar,p=psi.t,FUN=function(a,Q,g,p){
      as.numeric(A==a)/g * (Y - Q) + Q - p
    },SIMPLIFY=FALSE)
    PnDnoStar <- lapply(DnoStar, mean)
    
    if("g" %in% guard){
      if(reduction=="univariate"){
        DnQoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,gr=grnStar,FUN=function(a,Q,gr){
          as.numeric(A==a)/gr[[2]] * gr[[1]] * (Y-Q)
        }, SIMPLIFY=FALSE)
      }else if(reduction=="bivariate"){
        DnQoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,g=gnStar,gr=grnStar,FUN=function(a,Q,gr,g){
          as.numeric(A==a)/gr[[1]] * (gr[[1]] - g)/g * (Y-Q)
        },SIMPLIFY=FALSE)
      }
      PnDQnStar <- lapply(DnQoStar, mean)
    }
    if("Q" %in% guard){
      DngoStar <- mapply(a=split(a0,1:length(a0)),Qr=QrnStar,g=gnStar,FUN=function(a,Qr,g){
        Qr/g * (as.numeric(A==a) - g)
      }, SIMPLIFY=FALSE)
      PnDgnStar <- lapply(DngoStar, mean)
    }

    cat("TMLE Iteration", ct, "=", round(unlist(eps),5), "\n")
    cat("Mean of IC       =", round(c(unlist(PnDnoStar), unlist(PnDQnStar), unlist(PnDgnStar)), 10),"\n")
  }
  
  # standard tmle fluctuations
  QnStar1Out <- fluctuateQ1(Y=Y,A=A,W=W,Qn=Qn,gn=gn,a0=a0)
  QnStar1 <- llply(QnStar1Out, function(x){unlist(x[[1]])})
  
  # tmle estimates
  psi.t <- lapply(QnStar, mean)
  psi.t1 <- lapply(QnStar1, mean)
  
  # covariance for tmle
  Dno1Star <- mapply(a=split(a0,1:length(a0)),Q=QnStar1,g=gn,p=psi.n,FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  
  DnoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,g=gnStar,p=psi.n,FUN=function(a,Q,g,p){
    as.numeric(A==a)/g * (Y - Q) + Q - p
  },SIMPLIFY=FALSE)
  
  DnQoStar <- list(rep(0, length(Y)))
  DngoStar <- list(rep(0, length(Y)))
  
  if("g" %in% guard){
    if(reduction=="univariate"){
      DnQoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,gr=grnStar,FUN=function(a,Q,gr){
        as.numeric(A==a)/gr[[2]] * gr[[1]] * (Y-Q)
      }, SIMPLIFY=FALSE)
    }else if(reduction=="bivariate"){
      DnQoStar <- mapply(a=split(a0,1:length(a0)),Q=QnStar,g=gn,gr=grnStar,FUN=function(a,Q,gr,g){
        as.numeric(A==a)/gr[[1]] * (gr[[1]]-g)/g * (Y-Q)
      }, SIMPLIFY=FALSE)
    }
   
  }
  
  if("Q" %in% guard){
    DngoStar <- mapply(a=split(a0,1:length(a0)),Qr=QrnStar,g=gnStar,FUN=function(a,Qr,g){
      Qr/g * (as.numeric(A==a) - g)
    }, SIMPLIFY=FALSE)
  }
  
  Dno1StarMat <- matrix(unlist(Dno1Star), ncol=length(Y), nrow=length(a0), byrow=TRUE)
  DnoStarMat <- matrix(unlist(DnoStar) - unlist(DnQoStar) - unlist(DngoStar), ncol=length(Y), nrow=length(a0),byrow=TRUE)
  
  cov.t1 <- (Dno1StarMat - mean(Dno1StarMat))%*%t((Dno1StarMat - mean(Dno1StarMat)))/(length(Y)^2)
  cov.t <- (DnoStarMat - mean(DnoStarMat))%*%t((DnoStarMat - mean(DnoStarMat)))/(length(Y)^2)
  
  # final results
  if(returnModels){
  out <- list(naive=list(est=unlist(psi.n)),
              tmle=list(est=unlist(psi.t1),cov=cov.t1),
              tmle.dral=list(est=unlist(psi.t),cov=cov.t,Qrn=QrnStar,grn=grnStar),
              os=list(est=unlist(psi.o1), cov=cov.o1),
              os.dral=list(est=unlist(psi.o),cov=cov.o, Qrn=Qrn, grn=grn),
              QnMod=QnMod, gnMod=gnMod)
  }else{
    out <- list(naive=list(est=unlist(psi.n)),
                tmle=list(est=unlist(psi.t1),cov=cov.t1),
                tmle.dral=list(est=unlist(psi.t),cov=cov.t),
                os=list(est=unlist(psi.o1), cov=cov.o1),
                os.dral=list(est=unlist(psi.o),cov=cov.o))
  }
  out
}










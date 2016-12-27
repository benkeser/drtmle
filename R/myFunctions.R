#' drtmle 
#' 
#' Main function to estimate TMLE with doubly-robust inference
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
                      tolEps=1e-4, 
                      tolIC="default", 
                      tolg=1e-8,
                      verbose=TRUE,
                      Qsteps=2,
                      ...){

  if(tolIC=="default") tolIC <- sqrt(1/length(Y))/50

  # estimate g
  gnOut <- estimateG(A=A, W=W, tolg=tolg, verbose=verbose, returnModels=returnModels,libraryg=libraryg, 
    glmg=glmg, family="binomial", a0=a0)
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
  
  cov.o1 <- Dno1Mat%*%t(Dno1Mat)/(length(Y)^2)
  cov.o <- DnoMat%*%t(DnoMat)/(length(Y)^2)
  
  # initialize fluctuations
  QnStar <- Qn
  if("g" %in% guard) grnStar <- grn
  if("Q" %in% guard) QrnStar <- Qrn 
  gnStar <- gn; 
  PnDQnStar <- PnDgnStar <- PnDnoStar <- Inf
  eps <- list(Inf)
  ct <- 0

  # fluctuate
  while(max(abs(c(unlist(PnDQnStar),unlist(PnDgnStar),unlist(PnDnoStar)))) > tolIC & max(abs(unlist(eps)), na.rm=TRUE)>tolEps & ct < maxIter){
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
        QnStarOut <- fluctuateQ(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, grn=grnStar, reduction=reduction, family=family)
        QnStar <- llply(QnStarOut, function(x){unlist(x[[1]])})
        epsQ <- laply(QnStarOut, function(x){x$eps})
      }else if(Qsteps==2){
        # do the extra targeting
        QnStarOut2 <- fluctuateQ2(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, grn=grnStar, reduction=reduction, family=family)
        QnStar <- llply(QnStarOut2, function(x){unlist(x[[1]])})
        
        # update grn based on new Qn
        # grnStar <- estimategrn(Y=Y, A=A, W=W, reduction=reduction, tolg=tolg, a0=a0, Qn=QnStar, gn=gnStar, glmgr=glmgr, librarygr=librarygr)
        
        # do the usual targeting
        QnStarOut1 <- fluctuateQ1(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, family=family)
        QnStar <- llply(QnStarOut1, function(x){unlist(x[[1]])})
        
        # for later retrieval of fluct coefficients
        epsQ <- mapply(q1=QnStarOut1, q2=QnStarOut2, function(q1,q2){c(q1$eps,q2$eps)})
      }
    }else{
      QnStarOut <- fluctuateQ1(Y=Y, A=A, W=W, a0=a0, Qn=QnStar, gn=gnStar, family=family)
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
  QnStar1Out <- fluctuateQ1(Y=Y,A=A,W=W,Qn=Qn,gn=gn,a0=a0,family=family)
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
  
  cov.t1 <- Dno1StarMat%*%t(Dno1StarMat)/(length(Y)^2)
  cov.t <- DnoStarMat%*%t(DnoStarMat)/(length(Y)^2)
  
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

#' fluctuateQ1
#' 
#' Function called internally by drtmle to perform the first fluctuation 
#' of the initial estimator of Q (i.e., solves the original EIF estimating eqn)

fluctuateQ1 <- function(Y,A,W, Qn, gn, a0, family){
  QnStar <- mapply(a=a0,Q=Qn,g=gn,FUN=function(x, a, Q, g, gr){
      l <- min(Y); u <- max(Y)
      Yscale <- (Y-l)/(u-l)
      off <- SuperLearner:::trimLogit((Q-l)/(u-l))
      H1 <- as.numeric(A==a)/g
      suppressWarnings(
      fm <- glm(Yscale ~ -1 + offset(off) + H1, start=0,
                data=data.frame(Y=Y, off=off, H1=H1), family="binomial")
      )
      Qnstar <- predict(fm,type="response",
                       newdata=data.frame(off=off, H1=1/g))*(u-l) + l
      list(est=Qnstar,eps=fm$coef)
    }, SIMPLIFY=FALSE)
  QnStar
}

#' fluctuateQ2 
#' 
#' Function called internally by drtmle to perform the second fluctuation 
#' of the initial estimator of Q (i.e., solves the new estimating eqn that results
#' from misspecification of g)

fluctuateQ2 <- function(Y,A,W,Qn,gn,grn,a0,family,reduction,coefTol=1e7){
  QnStar <- mapply(a=a0,Q=Qn,g=gn,gr=grn,FUN=function(a, Q, g, gr){
    l <- min(Y); u <- max(Y)
    Yscale <- (Y-l)/(u-l)
    off <- SuperLearner:::trimLogit((Q-l)/(u-l))
    if(reduction=="univariate") H2 <- as.numeric(A==a)/gr[[2]] * gr[[1]]
    if(reduction=="bivariate") H2 <- as.numeric(A==a)/gr[[1]] * (gr[[1]]-g)/g
    suppressWarnings(
      fm <- glm(Yscale ~ -1 + offset(off) + H2, start=c(0),
                data=data.frame(Y=Y, off=off, H2=H2), family="binomial")
    )
    if(!fm$converged | abs(max(fm$coef)) > coefTol){
      # if it doesn't converge, try with no starting values
      suppressWarnings(
        fm <- glm(Yscale ~ -1 + offset(off) + H2, 
                  data=data.frame(Y=Y, off=off, H2=H2), family="binomial")
      )
      if(!fm$converged | abs(max(fm$coef)) > coefTol){
        warning("No sane fluctuation found. Proceeding using current estimates.")
        if(reduction=="univariate"){
          return(list(est=Q,eps=rep(0,2)))
        }else if(reduction=="bivariate"){
          return(list(est=Q,eps=rep(0,2)))
        }
      }
    }
    
    if(reduction=="univariate"){
      Qnstar <- predict(fm,type="response",newdata=data.frame(off=off, H2=1/gr[[2]] * gr[[1]]))*(u-l) + l
      list(est=Qnstar, 
           eps=fm$coef)
    }else if(reduction=="bivariate"){
      Qnstar <- predict(fm,type="response",newdata=data.frame(off=off, H2=1/gr[[1]] * (gr[[1]]-g)/g))*(u-l) + l
      list(est=Qnstar,
           eps=fm$coef)
    }
  }, SIMPLIFY=FALSE)
  QnStar
}


#' fluctuateG
#' 
#' Function called internally by drtmle to perform the fluctuation 
#' of the initial estimator of g (i.e., solves the new estimating eqn that results
#' from misspecification of Q)

fluctuateG <- function(Y, A, W, a0, Qn, gn, Qrn, tolg, coefTol=1e7){
  gnStar <- mapply(a=a0, Q=Qn, g=gn, Qr=Qrn, FUN=function(x, a, Q, g, Qr){
    H1 <- Qr/g
    off <- SuperLearner:::trimLogit(g, tolg)
    thisA <- as.numeric(A==a)
    suppressWarnings(
      fm <- glm(thisA ~ -1 + offset(off) + H1, start=0,
                data=data.frame(thisA=thisA, off=off, H1=H1), family="binomial")
    )
    if(!fm$converged | abs(fm$coef) > coefTol){
      suppressWarnings(
      fm <- glm(thisA ~ -1 + offset(off) + H1,
                data=data.frame(thisA=thisA, off=off, H1=H1), family="binomial")
      )
      if(!fm$converged | abs(fm$coef) > coefTol){
        warning("No sane fluctuation found for G this iteration. Check mean of IC.") 
      }
    }
    pred <- predict(fm, type="response")
    pred[pred < tolg] <- tolg
    list(est=pred, eps=fm$coef)
  }, SIMPLIFY = FALSE)
  gnStar
}


#' fluctuateQ 
#' 
#' Function called internally by drtmle to perform simultaneous fluctuation 
#' of the initial estimator of Q (i.e., solves both EIF estimating eqn and 
#' the new estimating eqn that results from misspecification of g)

fluctuateQ <- function(Y,A,W,Qn,gn,grn,a0,family,reduction,coefTol=1e7){
  QnStar <- mapply(a=a0,Q=Qn,g=gn,gr=grn,FUN=function(a, Q, g, gr){
    l <- min(Y); u <- max(Y)
    Yscale <- (Y-l)/(u-l)
    off <- SuperLearner:::trimLogit((Q-l)/(u-l))
    H1 <- as.numeric(A==a)/g
    if(reduction=="univariate") H2 <- as.numeric(A==a)/gr[[2]] * gr[[1]]
    if(reduction=="bivariate") H2 <- as.numeric(A==a)/gr[[1]] * (gr[[1]]-g)/g
    suppressWarnings(
    fm <- glm(Yscale ~ -1 + offset(off) + H1 + H2, start=c(0,0),
               data=data.frame(Y=Y, off=off, H1=H1, H2=H2), family="binomial")
    )
    if(!fm$converged | abs(max(fm$coef)) > coefTol){
      # if it doesn't converge, try with no starting values
      suppressWarnings(
      fm <- glm(Yscale ~ -1 + offset(off) + H1 + H2, 
                data=data.frame(Y=Y, off=off, H1=H1, H2=H2), family="binomial")
      )
      if(!fm$converged | abs(max(fm$coef)) > coefTol){
        warning("No sane fluctuation found. Proceeding using current estimates.")
        if(reduction=="univariate"){
          return(list(est=Q,eps=rep(0,2)))
        }else if(reduction=="bivariate"){
          return(list(est=Q,eps=rep(0,2)))
        }
      }
    }
    
    if(reduction=="univariate"){
      Qnstar <- predict(fm,type="response",newdata=data.frame(off=off, H1=1/g, H2=1/gr[[2]] * gr[[1]]))*(u-l) +l 
        list(est=Qnstar,
             eps=fm$coef)
      }else if(reduction=="bivariate"){
        Qnstar <- predict(fm,type="response",newdata=data.frame(off=off, H1=1/g, H2=1/gr[[1]] * (gr[[1]]-g)/g))*(u-l) + l
        list(est=Qnstar,
             eps=fm$coef)
      }
  }, SIMPLIFY=FALSE)
  QnStar
}

#' estimategrn 
#' 
#' Estimates the reduced dimension regressions necessary for the additional 
#' fluctuations. 

estimategrn <- function(Y, A, W, Qn, gn, librarygr, tolg, glmgr, a0, reduction){
  if(is.null(librarygr) & is.null(glmgr)) stop("Specify Super Learner library or GLM formula for gr")
  if(!is.null(librarygr) & !is.null(glmgr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glmgr <- NULL
  }
  # Super Learner
  if(!is.null(librarygr)){
    grn <- mapply(a=a0,Q=Qn,g=gn,SIMPLIFY=FALSE, FUN=function(a,Q,g){
      if(length(unique(Q))==1){
        warning("Only one unique value of Qn0. Proceeding with empicial mean for grn")
        if(reduction=="univariate"){
          grn1 <- rep(mean((as.numeric(A==a)-g)/g), length(Y))
          grn2 <- rep(mean(as.numeric(A==a)), length(Y))
          grn2[grn2 < tolg] <- tolg
        }else if(reduction=="bivariate"){
          grn <- rep(mean(as.numeric(A==a)), length(Y))
          grn[grn < tolg] <- tolg
        }
      }else{
        if(length(librarygr)>1){
          if(reduction=="univariate"){
            fm1 <- SuperLearner(Y=(as.numeric(A==a)-g)/g, X=data.frame(Qn=Q), 
                                family=gaussian(),SL.library=librarygr,  method="method.NNLS2")
            fm2 <- SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q), 
                                family=binomial(),SL.library=librarygr)
            if(!all(fm1$coef==0)){
              grn1 <- predict(fm1, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]              
            }else{
              grn1 <- predict(fm1, newdata=data.frame(Qn=Q), onlySL=FALSE)[[2]][,which(fm1$cvRisk == min(fm1$cvRisk,na.rm=TRUE))]            
            }

            if(!all(fm2$coef==0)){
              grn2 <- predict(fm2, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]  
            }else{
              grn2 <- predict(fm2, newdata=data.frame(Qn=Q), onlySL = FALSE)[[2]][,which(fm2$cvRisk == min(fm2$cvRisk,na.rm=TRUE))]
            }
            
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm1 <- SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q, gn=g), 
                                family=binomial(),SL.library=librarygr)
            if(!all(fm1$coef==0)){
              grn <- predict(fm1, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]
            }else{
              grn <- predict(fm1, newdata=data.frame(Qn=Q), onlySL=FALSE)[[2]][,which(fm1$cvRisk == min(fm1$cvRisk,na.rm=TRUE))]
            }
            grn[grn < tolg] <- tolg
          }
        }else if(length(librarygr)==1){
          if(reduction=="univariate"){
            obj1 <- do.call(librarygr, 
                            args=list(Y=(as.numeric(A==a)-g)/g,X=data.frame(Qn=Q),
                                      obsWeights=rep(1, length(A)),
                                      newX=data.frame(Qn=Q), family=gaussian()))
            grn1 <- predict(object=obj1$fit, newdata=data.frame(Qn=Q))
            obj2 <- do.call(librarygr, args=list(Y=as.numeric(A==a), X=data.frame(Qn=Q),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q), family=binomial()))
            grn2 <- predict(object=obj2$fit, newdata=data.frame(Qn=Q))
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            obj <- do.call(librarygr, args=list(Y=as.numeric(A==a),X=data.frame(Qn=Q,gn=g),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q, gn=g), family=binomial()))
            grn <- predict(object=obj$fit, newdata=data.frame(Qn=Q,gn=g))
            grn[grn < tolg] <- tolg
          }
          
        }
      }
      if(reduction=="univariate") return(list(grn1,grn2))
      if(reduction=="bivariate") return(list(grn))
    })
  }
     
  # GLM
  if(!is.null(glmgr)){
    grn <- mapply(a=a0,Q=Qn,g=gn,SIMPLIFY=FALSE, FUN=function(a,Q,g){
      if(length(unique(Qn[[1]]))==1){
        warning("Only one unique value of Qn0. Proceeding with empicial mean for grn")
        grn1 <- rep(mean((as.numeric(A==a)-g)/g), length(Y))
        grn2 <- rep(mean(as.numeric(A==a)), length(Y))
      }else{
        fm1 <- glm(as.formula(paste0("grn1~",glmgr)), family="gaussian",
                   data=data.frame(grn1=(as.numeric(A==a)-g)/g, Qn=Q))
        fm2 <- glm(as.formula(paste0("A~",glmgr)), family="binomial",
                   data=data.frame(A=A, Qn=Q))
        grn1 <- predict(fm1, type="response")
        grn2 <- predict(fm2, type="response")
        grn2[grn2 < tolg] <- tolg
      }
      list(grn1,grn2)
    })
  }
  
  grn
}

#' estimateQrn 
#' 
#' Estimates the reduced dimension regressions necessary for the  
#' fluctuations of g

estimateQrn  <- function(Y, A, W, Qn, gn, glmQr, libraryQr, a0){
  if(is.null(libraryQr) & is.null(glmQr)) stop("Specify Super Learner library or GLM formula for Qr")
  if(!is.null(libraryQr) & !is.null(glmQr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glmQr <- NULL
  }
  # Super Learner
  if(!is.null(libraryQr)){
    Qrn <- mapply(a=a0, g=gn, Q=Qn, SIMPLIFY=FALSE, FUN=function(a,g,Q){
      if(length(unique(g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        rep(mean((Y-Q)[A==a]), length(Y))
      }else{
        if(length(libraryQr)>1){
          suppressWarnings(
          fm <- SuperLearner(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                             family=gaussian(),SL.library=libraryQr, method="method.NNLS2")
          )
          # if all weights = 0, use discrete SL
          if(!all(fm$coef==0)){
            predict(fm, newdata=data.frame(gn=g), onlySL=TRUE)[[1]]
          }else{
            predict(fm, newdata=data.frame(gn=g), onlySL=FALSE)[[2]][,which(fm$cvRisk == min(fm$cvRisk, na.rm = TRUE))]
          }
        }else if(length(libraryQr)==1){
          obj <- do.call(libraryQr, args=list(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),family=gaussian(),
                                              newX=data.frame(gn=g[A==a]),
                                              obsWeights=rep(1, length(Y[A==a]))))
          pred <- predict(object=obj$fit, newdata=data.frame(gn=g))
          pred
        }
      }
    })
  }
  
  # GLM
  if(!is.null(glmQr)){
    Qrn <- mapply(a=a0, g=gn, Q=Qn, SIMPLIFY=FALSE, FUN=function(a,g,Q){
      if(length(unique(g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        rep(mean((Y-Q)[A==a]), length(Y))
      }else{
        fm <- glm(as.formula(paste0("Qrn ~",glmQr)), 
                  data=data.frame(Qrn=(Y-Q)[A==a], gn=g[A==a]),
                  family="gaussian")
        predict(fm, newdata=data.frame(gn=g),type="response")
      }
    })
  }
  Qrn 
}

#' estimateG
#' 
#' Function to estimate propensity score

estimateG <- function(A, W,libraryg,glmg,family,a0,tolg,verbose=FALSE, returnModels=FALSE){
  if(is.null(libraryg) & is.null(glmg)) stop("Specify Super Learner library or GLM formula for g")
  if(!is.null(libraryg) & !is.null(glmg)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glmg <- NULL
  }
  # Super Learner
  if(!is.null(libraryg)){
    if(length(libraryg)>1 | is.list(libraryg)){
      if(length(a0)==length(unique(A)) & length(unique(A))==2){
        fm <- SuperLearner(Y=as.numeric(A==a0[1]), X=W, family="binomial",SL.library=libraryg,verbose=verbose)
        pred <- predict(fm, onlySL=TRUE)[[1]]
        pred[pred < tolg] <- tolg
        gn <- vector(mode="list",length=2)
        gn[[1]] <- pred; gn[[2]] <- 1-pred
      }else{
        gn <- alply(a0, 1, function(x,A,W,libraryg){
          fm <- SuperLearner(Y=as.numeric(A==x), X=W, family="binomial",SL.library=libraryg,verbose=verbose)
          pred <- predict(fm, onlySL=TRUE)[[1]]
          pred[pred < tolg] <- tolg
          pred
        }, A=A, W=W, libraryg=libraryg)
      }
    }else if(!is.list(libraryg) & length(libraryg)==1){
      gn <- alply(a0, 1, function(x){
        fm <- do.call(libraryg, args=list(Y=as.numeric(A==x), X=W, newX=W, obsWeights=rep(1,length(A)),family=binomial()))
        pred <- predict(object=fm$fit,newdata=W)
        pred[pred < tolg] <- tolg
        list(pred=pred,gnMod=fm)
      })
      fm <- lapply(gn, function(x){x[[2]]})
      gn <- lapply(gn, function(x){x[[1]]})
    }
  }
  
  # GLM
  if(!is.null(glmg)){
    gn  <- alply(a0,1,function(x){
      thisDat <- data.frame(thisA=as.numeric(A==x), W=W)
      colnames(thisDat) <- c("A",colnames(W))
      thisA <- as.numeric(A==x)
      fm <- glm(as.formula(paste0("thisA~",glmg)), data=thisDat, family="binomial")
      pred <- predict(fm, type="response")
      pred[pred < tolg] <- tolg
    })
  }
  if(returnModels){
    return(list(gn, fm))
  }else{
    return(gn)
  }
}

#' estimateQ 
#' 
#' Function to estimate initial outcome regression


estimateQ <- function(Y,A,W,libraryQ,glmQ,a0,stratify,family,verbose=FALSE,returnModels=FALSE,...){
  if(is.null(libraryQ) & is.null(glmQ)) stop("Specify Super Learner library or GLM formula for Q")
  if(!is.null(libraryQ) & !is.null(glmQ)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glmQ <- NULL
  }
  # Super Learner
  if(!is.null(libraryQ)){
    if(!stratify){
      if(length(libraryQ)>1 | is.list(libraryQ)){
        fm <- SuperLearner(Y=Y, X=data.frame(A,W), verbose=verbose,family=family,SL.library=libraryQ,...)
      
        Qn <- alply(a0,1,function(x){
          predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]]
        })
      }else if(length(libraryQ)==1){
        fm <- do.call(libraryQ, args=list(Y=Y, X=data.frame(A,W), verbose=verbose, newX=data.frame(A,W),
                                          obsWeights=rep(1,length(A)),
                                          family=family))
        Qn <- alply(a0,1,function(x){
          predict(object=fm$fit, newdata=data.frame(A=x, W))
        })
      }
   }else{
     if(length(libraryQ)>1 | is.list(libraryQ)){
      Qn <- alply(a0,1,function(x){
        fm <- SuperLearner(Y=Y[A==x], X=W[A==x,], verbose=verbose,family=family,SL.library=libraryQ)
        predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]]
      })
     }else if(length(libraryQ)==1){
       Qn <- alply(a0,1,function(x){
         fm <- do.call(libraryQ, args=list(Y=Y[A==x], X=W[A==x,], newX=W[A==x,],verbose=verbose,obsWeights=rep(1,sum(A==x)),
                                           family=family))
         predict(object=fm$fit, newdata=data.frame(W))
       })
     }
    }
  }
  
  # GLM
  if(!is.null(glmQ)){
    thisDat <- data.frame(Y=Y, A=A, W=W); colnames(thisDat) <- c("Y","A",colnames(W))
    if(!stratify){
      fm <- glm(as.formula(paste0("Y~",glmQ)), data=thisDat, family=family)
      Qn <- alply(matrix(a0),1,function(x,fm){
        predict(fm, newdata=data.frame(A=x,W), type="response")
      },fm=fm)
    }else{
      Qn <- alply(matrix(a0),1,function(x){
        fm <- glm(as.formula(paste0("Y~",glmQ)), data=thisDat[A==x,], family=family)
        predict(fm, newdata=data.frame(A=x,W), type="response")
      })     
    }
  }
  if(returnModels) list(Qn,fm)
  else Qn
}

#' SL.caret1 
#' 
#' A modification of the SL.caret function from the SuperLearner package that
#' suppresses some output when method = "gbm".
#' @export

SL.caret1 <- function (Y, X, newX, family, obsWeights, method = "rf", tuneLength = 3, 
                       trControl = caret::trainControl(method = "cv", number = ifelse(length(Y)<=100, 20, 10), verboseIter = FALSE), 
                       metric,...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,verbose=FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl, verbose = FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        # pass trace==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}

#' SL.rpart.caret1 
#' 
#' Uses SL.caret1 to train rpart
#' @export

SL.rpart.caret1 <- function(...,method="rpart",tuneLength = 8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.rf.caret1
#' 
#' Uses SL.caret1 to train random forest
#' @export

SL.rf.caret1 <- function(...,method="rf",tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.gbm.caret1 
#' 
#' Uses SL.caret1 to train a gbm
#' @export

SL.gbm.caret1 <- function(...,method="gbm",tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.gbm.caret2
#' Uses SL.caret 1 to train a gbm with 10 choices of tuning parameters and 5-fold CV
#' @export

SL.gbm.caret2 <- function (..., method = "gbm", tuneLength = 20, trControl = caret::trainControl(method = "cv", number = 5)) 
{
    SL.caret1(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.rf.caret2
#' Uses SL.caret 1 to train a random forest with 10 choices of tuning parameters and 5-fold CV
#' @export

SL.rf.caret2 <- function (..., method = "rf", tuneLength = 20, trControl = caret::trainControl(method = "cv", number = 5)) 
{
    SL.caret1(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.svmLinear.caret1 
#' 
#' Uses SL.caret1 to train a linear svm 
#' @export

SL.svmLinear.caret1 <- function(...,method="svmLinear",tuneLength = 8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.gamSpline.caret1
#' 
#' Uses SL.caret1 to train a gam
#' @export 
SL.gamSpline.caret1 <- function(...,method="gamSpline",tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}


#' SL.nnet.caret1
#' 
#' Uses SL.caret1 to train a neural net
#' @export 
SL.nnet.caret1 <- function(...,method="nnet", tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.glmnet.caret1
#' 
#' Uses SL.caret1 to train a elastic net regression
#' @export 
SL.glmnet.caret1 <- function(...,method="glmnet", tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.randomGLM.caret1
#' 
#' Uses SL.caret1 to train a random GLM
#' @export 
SL.randomGLM.caret1 <- function(...,method="randomGLM",tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

#' SL.npreg
#' 
#' SuperLearner wrapper to fit kernel regression
#' @export 
SL.npreg <- function (Y, X, newX, family, obsWeights, 
                      rangeThresh=1e-7, ...) 
{
  if(abs(diff(range(Y))) <= rangeThresh){
    thisMod <- glm(Y ~ 1, data=X)
  }else{
  bw <- np::npregbw(as.formula(paste("Y ~", paste(names(X),collapse="+"))), data=X,
                ftol=0.01, tol=0.01, remin=FALSE)
  
  # fit the kernel regression
  thisMod <- np::npreg(bw)
  }
  
  pred <- predict(thisMod, newdata=newX)
  fit <- list(object = thisMod)
  class(fit) <- "SL.npreg"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' predict.SL.npreg
#' 
#' Method for predicting SL.npreg
#' @export
predict.SL.npreg <- function(object, newdata, ...){
  pred <- predict(object=object$object, newdata=newdata)
  pred
}
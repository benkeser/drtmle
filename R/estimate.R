#' estimateG
#' 
#' Function to estimate propensity score
#' 
#' @importFrom plyr alply
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param libraryg A vector of characters or a list describing the Super Learner library to be used 
#' for the propensity score
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' @param glmg A character describing a formula to be used in the call to \code{glm} for the propensity score
#' @param a0 A list of fixed treatment values 
#' 
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 

estimateG <- function(A, W,libraryg,glmg,a0,tolg,verbose=FALSE, returnModels=FALSE){
  if(is.null(libraryg) & is.null(glmg)) stop("Specify Super Learner library or GLM formula for g")
  if(!is.null(libraryg) & !is.null(glmg)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glmg <- NULL
  }
  # Super Learner
  if(!is.null(libraryg)){
    if(length(libraryg)>1 | is.list(libraryg)){
      if(length(a0)==length(unique(A)) & length(unique(A))==2){
        fm <- SuperLearner::SuperLearner(Y=as.numeric(A==a0[1]), X=W, 
                                         family=binomial(),SL.library=libraryg,verbose=verbose)
        pred <- stats::predict(fm, onlySL=TRUE)[[1]]
        pred[pred < tolg] <- tolg
        gn <- vector(mode="list",length=2)
        gn[[1]] <- pred; gn[[2]] <- 1-pred
      }else{
        gn <- plyr::alply(a0, 1, function(x,A,W,libraryg){
          fm <- SuperLearner::SuperLearner(Y=as.numeric(A==x), X=W, 
                                           family=binomial(),SL.library=libraryg,verbose=verbose)
          pred <- stats::predict(fm, onlySL=TRUE)[[1]]
          pred[pred < tolg] <- tolg
          pred
        }, A=A, W=W, libraryg=libraryg)
      }
    }else if(!is.list(libraryg) & length(libraryg)==1){
      gn <- plyr::alply(a0, 1, function(x){
        fm <- do.call(libraryg, args=list(Y=as.numeric(A==x), X=W, 
                                          newX=W, obsWeights=rep(1,length(A)),
                                          family=binomial()))
        pred <- stats::predict(object=fm$fit,newdata=W)
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
      fm <- stats::glm(stats::as.formula(paste0("thisA~",glmg)), data=thisDat, family="binomial")
      pred <- stats::predict(fm, type="response")
      pred[pred < tolg] <- tolg
      pred
    })
  }
  out <- list(est = gn, fm = NULL)
  if(returnModels){
    out$fm <- fm
  }
  return(out)
}



#' estimateQ 
#' 
#' Function to estimate initial outcome regression
#' 
#' @param Y A vector of continuous or binary outcomes. 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param libraryQ A vector of characters or a list describing the Super Learner library to be used 
#' for the outcome regression
#' @param verbose A boolean indicating whether to print status updates
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions
#' @param glmQ A character describing a formula to be used in the call to \code{glm} for the outcome regression
#' @param a0 A list of fixed treatment values 
#' @param family A character passed to \code{SuperLearner}
#' @param stratify A \code{boolean} indicating whether to estimate the outcome regression separately
#' for observations with \code{A} equal to 0/1 (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param ... Additional arguments (not currently used) 
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 


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
        fm <- SuperLearner::SuperLearner(Y=Y, X=data.frame(A,W), verbose=verbose,family=family,SL.library=libraryQ,...)
      
        Qn <- alply(a0,1,function(x){
          stats::predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]]
        })
      }else if(length(libraryQ)==1){
        fm <- do.call(libraryQ, args=list(Y=Y, X=data.frame(A,W), verbose=verbose, newX=data.frame(A,W),
                                          obsWeights=rep(1,length(A)),
                                          family=family))
        Qn <- alply(a0,1,function(x){
          stats::predict(object=fm$fit, newdata=data.frame(A=x, W))
        })
      }
   }else{
     if(length(libraryQ)>1 | is.list(libraryQ)){
      tmp <- plyr::alply(a0,1,function(x){
        fm <- SuperLearner::SuperLearner(Y=Y[A==x], X=W[A==x,], verbose=verbose,family=family,SL.library=libraryQ)
        list(est = predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]],
             fm = fm)
      })
      Qn <- lapply(tmp,"[[",1)
      fm <- lapply(tmp,"[[",2)
     }else if(length(libraryQ)==1){
       tmp <- plyr::alply(a0,1,function(x){
         fm <- do.call(libraryQ, args=list(Y=Y[A==x], X=W[A==x,], newX=W[A==x,],verbose=verbose,obsWeights=rep(1,sum(A==x)),
                                           family=family))
         list(est = stats::predict(object=fm$fit, newdata=data.frame(W)),
              fm = fm)
       })
       Qn <- lapply(tmp,"[[",1)
       fm <- lapply(tmp,"[[",2)
     }
    }
  }
  
  # GLM
  if(!is.null(glmQ)){
    thisDat <- data.frame(Y=Y, A=A, W); colnames(thisDat) <- c("Y","A",colnames(W))
    if(!stratify){
      fm <- stats::glm(stats::as.formula(paste0("Y~",glmQ)), data=thisDat, family=family)
      Qn <- plyr::alply(matrix(a0),1,function(x,fm){
        stats::predict(fm, newdata=data.frame(A=x,W), type="response")
      },fm=fm)
    }else{
      tmp <- plyr::alply(matrix(a0),1,function(x){
        fm <- stats::glm(stats::as.formula(paste0("Y~",glmQ)), data=thisDat[A==x,], family=family)
        list(est = stats::predict(fm, newdata=data.frame(A=x,W), type="response"),
             fm = fm)
      })     
      Qn <- lapply(tmp,"[[",1)
      fm <- lapply(tmp,"[[",2)
    }
  }
  out <- list(est = Qn, fm = NULL)
  if(returnModels) out$fm <- fm
  return(out)
}


#' estimateQrn 
#' 
#' Estimates the reduced dimension regressions necessary for the  
#' fluctuations of g
#' 
#' 
#' @param Y A vector of continuous or binary outcomes. 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param libraryQr A vector of characters or a list describing the Super Learner library to be used 
#' for the first reduced-dimension regression.
#' @param glmQr A character describing a formula to be used in the call to \code{glm} for the first reduced-dimension regression. Ignored
#' if \code{librarygr!=NULL}.
#' @param a0 A list of fixed treatment values 
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimateQrn  <- function(Y, A, W, Qn, gn, glmQr, libraryQr, a0, returnModels){
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
        m1 <- mean((Y-Q)[A==a])
        est <- rep(m1, length(Y))
        fm <- list(object = m1)
        class(fm) <- "SL.mean"
      }else{
        if(length(libraryQr)>1){
          suppressWarnings(
          fm <- SuperLearner::SuperLearner(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                             family=gaussian(),SL.library=libraryQr, method="method.CC_LS")
          )
          # if all weights = 0, use discrete SL
          if(!all(fm$coef==0)){
            est <- stats::predict(fm, newdata=data.frame(gn=g), onlySL=TRUE)[[1]]
          }else{
            est <- stats::predict(fm, newdata=data.frame(gn=g), onlySL=FALSE)[[2]][,which(fm$cvRisk == min(fm$cvRisk, na.rm = TRUE))]
          }
        }else if(length(libraryQr)==1){
          fm <- do.call(libraryQr, args=list(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                                             family=gaussian(),
                                              newX=data.frame(gn=g[A==a]),
                                              obsWeights=rep(1, length(Y[A==a]))))
          est <- stats::predict(object=fm$fit, newdata=data.frame(gn=g))
        }
      }
      out <- list(est = est, fm = NULL)
      if(returnModels) out$fm <- fm
      return(out)
    })
  }
  
  # GLM
  if(!is.null(glmQr)){
    Qrn <- mapply(a=a0, g=gn, Q=Qn, SIMPLIFY=FALSE, FUN=function(a,g,Q){
      if(length(unique(g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        glmQr <- "1"
      }
      fm <- stats::glm(stats::as.formula(paste0("Qrn ~",glmQr)), 
                data=data.frame(Qrn=(Y-Q)[A==a], gn=g[A==a]),
                family="gaussian")
      est <- stats::predict(fm, newdata=data.frame(gn=g),type="response")
      out <- list(est = est, fm = NULL)
      if(returnModels) out$fm <- fm
      return(out)
    })
  }
  tmp1 <- lapply(Qrn,function(x){ x$est })
  tmp2 <- lapply(Qrn,function(x){ fm = x$fm })
  return(list(est = tmp1, fm = tmp2))
  Qrn 
}

#' estimategrn 
#' 
#' Estimates the reduced dimension regressions necessary for the additional 
#' fluctuations. 
#' 
#' @param Y A vector of continuous or binary outcomes. 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param librarygr A vector of characters or a list describing the Super Learner library to be used 
#' for the second reduced-dimension regression.
#' @param glmgr A character describing a formula to be used in the call to \code{glm} for the second reduced-dimension regression. Ignored
#' if \code{librarygr!=NULL}.
#' @param reduction A character equal to \code{"univariate"} for a univariate misspecification correction or \code{"bivariate"}
#' for the bivariate version. 
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param a0 A list of fixed treatment values 
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimategrn <- function(Y, A, W, Qn, gn, librarygr, tolg, glmgr, a0, reduction,returnModels){
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
          m1 <- mean((as.numeric(A==a)-g)/g)
          grn1 <- rep(m1, length(Y))
          m2 <- mean(as.numeric(A==a))
          grn2 <- rep(m2, length(Y))
          grn2[grn2 < tolg] <- tolg
          fm1 <- list(object = m1)
          class(fm1) <- "SL.mean"
          fm2 <- list(object = m2)
          class(fm1) <- "SL.mean"
        }else if(reduction=="bivariate"){
          m2 <- mean(as.numeric(A==a))
          grn2 <- rep(m2, length(Y))
          grn2[grn2 < tolg] <- tolg
          fm2 <- list(object = m2)
          class(fm2) <- "SL.mean"
          fm1 <- NULL; grn1 <- NULL
        }
      }else{
        if(length(librarygr)>1){
          if(reduction=="univariate"){
            fm1 <- SuperLearner::SuperLearner(Y=(as.numeric(A==a)-g)/g, X=data.frame(Qn=Q), 
                                family=gaussian(),SL.library=librarygr,  method="method.NNLS2")
            fm2 <- SuperLearner::SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q), 
                                family=binomial(),SL.library=librarygr)
            if(!all(fm1$coef==0)){
              grn1 <- stats::predict(fm1, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]              
            }else{
              grn1 <- stats::predict(fm1, newdata=data.frame(Qn=Q), onlySL=FALSE)[[2]][,which(fm1$cvRisk == min(fm1$cvRisk,na.rm=TRUE))]            
            }

            if(!all(fm2$coef==0)){
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]  
            }else{
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q), onlySL = FALSE)[[2]][,which(fm2$cvRisk == min(fm2$cvRisk,na.rm=TRUE))]
            }
            
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm2 <- SuperLearner::SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q, gn=g), 
                                family=binomial(),SL.library=librarygr)
            if(!all(fm2$coef==0)){
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q, gn = g), onlySL=TRUE)[[1]]
            }else{
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q, gn = g), onlySL=FALSE)[[2]][,which(fm2$cvRisk == min(fm2$cvRisk,na.rm=TRUE))]
            }
            grn2[grn2 < tolg] <- tolg
            fm1 <- NULL; grn1 <- NULL
          }
        }else if(length(librarygr)==1){
          if(reduction=="univariate"){
            fm1 <- do.call(librarygr, 
                            args=list(Y=(as.numeric(A==a)-g)/g,X=data.frame(Qn=Q),
                                      obsWeights=rep(1, length(A)),
                                      newX=data.frame(Qn=Q), family=gaussian()))
            grn1 <- predict(object=fm1$fit, newdata=data.frame(Qn=Q))
            fm2 <- do.call(librarygr, args=list(Y=as.numeric(A==a), X=data.frame(Qn=Q),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q), family=binomial()))
            grn2 <- predict(object=fm2$fit, newdata=data.frame(Qn=Q))
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm2 <- do.call(librarygr, args=list(Y=as.numeric(A==a),X=data.frame(Qn=Q,gn=g),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q, gn=g), family=binomial()))
            grn2 <- predict(object=fm2$fit, newdata=data.frame(Qn=Q,gn=g))
            grn2[grn2 < tolg] <- tolg
            fm1 <- NULL; grn1 <- NULL
          }
          
        }
      }
      out <- list(grn1 = grn1, grn2 = grn2, fm1 = NULL, fm2 = NULL)
      if(returnModels){
        out$fm1 <- fm1; out$fm2 <- fm2
      }
      return(out)
    })
  }
     
  # GLM
  if(!is.null(glmgr)){
    grn <- mapply(a=a0,Q=Qn,g=gn,SIMPLIFY=FALSE, FUN=function(a,Q,g){
      if(length(unique(Qn[[1]]))==1){
        glmgr <- "1"
      }
      fm2 <- stats::glm(stats::as.formula(paste0("A~",glmgr)), family="binomial",
           data=data.frame(A=A, Qn=Q))
      grn2 <- stats::predict(fm2, type="response")
      grn2[grn2 < tolg] <- tolg
      if(reduction == "univariate"){
        fm1 <- stats::glm(stats::as.formula(paste0("grn1~",glmgr)), family="gaussian",
           data=data.frame(grn1=(as.numeric(A==a)-g)/g, Qn=Q))
        grn1 <- stats::predict(fm1, type="response")
      }else{
        fm1 <- grn1 <- NULL
      }
      out <- list(grn1 = grn1, grn2 = grn2, fm1 = NULL, fm2 = NULL)
      if(returnModels){
        out$fm1 <- fm1; out$fm2 <- fm2
      }
      return(out)
    })
  }
  tmp1 <- lapply(grn,function(x){ list(grn1 = x$grn1, grn2 = x$grn2) })
  tmp2 <- lapply(grn,function(x){ list(fm1 = x$fm1, fm2 = x$fm2) })
  return(list(est = tmp1, fm = tmp2))
}

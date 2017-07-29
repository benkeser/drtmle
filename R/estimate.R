#' estimateG
#' 
#' Function to estimate propensity score
#' 
#' @importFrom plyr alply
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param W A \code{data.frame} of named covariates
#' @param SL_g A vector of characters or a list describing the Super Learner library to be used 
#' for the propensity score
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' @param glm_g A character describing a formula to be used in the call to \code{glm} for the propensity score
#' @param a_0 A list of fixed treatment values 
#' 
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 

estimateG <- function(A, W,SL_g,glm_g,a_0,tolg,verbose=FALSE, returnModels=FALSE){
  if(is.null(SL_g) & is.null(glm_g)) stop("Specify Super Learner library or GLM formula for g")
  if(!is.null(SL_g) & !is.null(glm_g)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_g <- NULL
  }
  # Super Learner
  if(!is.null(SL_g)){
    if(length(SL_g)>1 | is.list(SL_g)){
      if(length(a_0)==length(unique(A)) & length(unique(A))==2){
        fm <- SuperLearner::SuperLearner(Y=as.numeric(A==a_0[1]), X=W, 
                                         family=stats::binomial(),
                                         SL.library=SL_g,verbose=verbose,
                                         method = "method.CC_nloglik")
        pred <- stats::predict(fm, onlySL=TRUE)[[1]]
        pred[pred < tolg] <- tolg
        gn <- vector(mode="list",length=2)
        gn[[1]] <- pred; gn[[2]] <- 1-pred
      }else{
        tmp <- plyr::alply(a_0, 1, function(x,A,W,SL_g){
          fm <- SuperLearner::SuperLearner(Y=as.numeric(A==x), X=W, 
                                           family=stats::binomial(),SL.library=SL_g,verbose=verbose)
          pred <- stats::predict(fm, onlySL=TRUE)[[1]]
          pred[pred < tolg] <- tolg
          list(est = pred, fm = fm)
        }, A=A, W=W, SL_g=SL_g)
        gn <- lapply(tmp,"[[",1)
        fm <- lapply(tmp,"[[",2)
      }
    }else if(!is.list(SL_g) & length(SL_g)==1){
      tmp <- plyr::alply(a_0, 1, function(x){
        fm <- do.call(SL_g, args=list(Y=as.numeric(A==x), X=W, 
                                          newX=W, obsWeights=rep(1,length(A)),
                                          family=stats::binomial()))
        pred <- stats::predict(object=fm$fit,newdata=W)
        pred[pred < tolg] <- tolg
        list(est=pred,fm=fm)
      })
      fm <- lapply(tmp, function(x){x[[2]]})
      gn <- lapply(tmp, function(x){x[[1]]})
    }
  }
  
  # GLM
  if(!is.null(glm_g)){
    if(length(a_0)==length(unique(A)) & length(unique(A))==2){
      thisDat <- data.frame(A = as.numeric(A==a_0[1]), W)
      fm <- stats::glm(stats::as.formula(paste0("A~",glm_g)), data=thisDat, 
                         family=stats::binomial())
      pred <- stats::predict(fm, type = "response")
      pred[pred < tolg] <- tolg; pred[pred > 1 - tolg] <- 1 - tolg
      gn <- vector(mode="list",length=2)
      gn[[1]] <- pred; gn[[2]] <- 1-pred
    }else{
      tmp <- alply(a_0,1,function(x){
        thisDat <- data.frame(thisA=as.numeric(A==x), W=W)
        colnames(thisDat) <- c("A",colnames(W))
        fm <- stats::glm(stats::as.formula(paste0("A~",glm_g)), data=thisDat, 
                         family=stats::binomial())
        pred <- stats::predict(fm, type="response")

        pred[pred < tolg] <- tolg
        list(est = pred, fm = fm)
      })
      gn <- lapply(tmp,"[[",1)
      fm <- lapply(tmp,"[[",2)
    }

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
#' @param SL_Q A vector of characters or a list describing the Super Learner library to be used 
#' for the outcome regression
#' @param verbose A boolean indicating whether to print status updates
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions
#' @param glm_Q A character describing a formula to be used in the call to \code{glm} for the outcome regression
#' @param a_0 A list of fixed treatment values 
#' @param family A character passed to \code{SuperLearner}
#' @param stratify A \code{boolean} indicating whether to estimate the outcome regression separately
#' for observations with \code{A} equal to 0/1 (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param ... Additional arguments (not currently used) 
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 


estimateQ <- function(Y,A,W,SL_Q,glm_Q,a_0,stratify,family,verbose=FALSE,returnModels=FALSE,...){
  if(is.null(SL_Q) & is.null(glm_Q)) stop("Specify Super Learner library or GLM formula for Q")
  if(!is.null(SL_Q) & !is.null(glm_Q)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_Q <- NULL
  }
  # Super Learner
  if(!is.null(SL_Q)){
    if(!stratify){
      if(length(SL_Q)>1 | is.list(SL_Q)){
        fm <- SuperLearner::SuperLearner(Y=Y, X=data.frame(A,W), 
                                         verbose=verbose,family=family,
                                         SL.library=SL_Q,
                                         method =ifelse(family$family=="binomial",
                                                        "method.CC_LS","method.CC_LS"),
                                         ...)
      
        Qn <- alply(a_0,1,function(x){
          stats::predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]]
        })
      }else if(length(SL_Q)==1){
        fm <- do.call(SL_Q, args=list(Y=Y, X=data.frame(A,W), verbose=verbose, newX=data.frame(A,W),
                                          obsWeights=rep(1,length(A)),
                                          family=family))
        Qn <- alply(a_0,1,function(x){
          stats::predict(object=fm$fit, newdata=data.frame(A=x, W))
        })
      }
   }else{
     if(length(SL_Q)>1 | is.list(SL_Q)){
      tmp <- plyr::alply(a_0,1,function(x){
        fm <- SuperLearner::SuperLearner(Y=Y[A==x], X=W[A==x,], 
                                         verbose=verbose,family=family,
                                         SL.library=SL_Q,
                                         method = "method.CC_LS")
        list(est = predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]],
             fm = fm)
      })
      Qn <- lapply(tmp,"[[",1)
      fm <- lapply(tmp,"[[",2)
     }else if(length(SL_Q)==1){
       tmp <- plyr::alply(a_0,1,function(x){
         fm <- do.call(SL_Q, args=list(Y=Y[A==x], X=W[A==x,], newX=W[A==x,],
                                           verbose=verbose,obsWeights=rep(1,sum(A==x)),
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
  if(!is.null(glm_Q)){
    thisDat <- data.frame(Y=Y, A=A, W); colnames(thisDat) <- c("Y","A",colnames(W))
    if(!stratify){
      fm <- stats::glm(stats::as.formula(paste0("Y~",glm_Q)), data=thisDat, family=family)
      Qn <- plyr::alply(matrix(a_0),1,function(x,fm){
        stats::predict(fm, newdata=data.frame(A=x,W), type="response")
      },fm=fm)
    }else{
      tmp <- plyr::alply(matrix(a_0),1,function(x){
        fm <- stats::glm(stats::as.formula(paste0("Y~",glm_Q)), data=thisDat[A==x,], family=family)
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
#' @param Qn A list of outcome regression estimates evaluated on observed data. If 
#' NULL then 0 is used for all Qn (as is needed to estimate reduced dimension regression
#' for sl_iptw)
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param SL_Qr A vector of characters or a list describing the Super Learner library to be used 
#' for the first reduced-dimension regression.
#' @param glm_Qr A character describing a formula to be used in the call to \code{glm} for the first reduced-dimension regression. Ignored
#' if \code{SL_gr!=NULL}.
#' @param a_0 A list of fixed treatment values 
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimateQrn  <- function(Y, A, W, Qn, gn, glm_Qr, SL_Qr, a_0, returnModels){
  if(is.null(Qn)){
    Qn <- vector(mode = "list", length = length(a_0))
    for(i in 1:length(a_0)){ Qn[[i]] <- rep(0, length(Y)) }
  }
  if(is.null(SL_Qr) & is.null(glm_Qr)) stop("Specify Super Learner library or GLM formula for Qr")
  if(!is.null(SL_Qr) & !is.null(glm_Qr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_Qr <- NULL
  }
  # Super Learner
  if(!is.null(SL_Qr)){
    Qrn <- mapply(a=a_0, g=gn, Q=Qn, SIMPLIFY=FALSE, FUN=function(a,g,Q){
      if(length(unique(g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        m1 <- mean((Y-Q)[A==a])
        est <- rep(m1, length(Y))
        fm <- list(object = m1)
        class(fm) <- "SL.mean"
      }else{
        if(length(SL_Qr)>1){
          suppressWarnings(
          fm <- SuperLearner::SuperLearner(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                             family=stats::gaussian(),
                             SL.library=SL_Qr, 
                             method="method.CC_LS")
          )
          # if all weights = 0, use discrete SL
          if(!all(fm$coef==0)){
            est <- stats::predict(fm, newdata=data.frame(gn=g), onlySL=TRUE)[[1]]
          }else{
            est <- stats::predict(fm, newdata=data.frame(gn=g), onlySL=FALSE)[[2]][,which(fm$cvRisk == min(fm$cvRisk, na.rm = TRUE))]
          }
        }else if(length(SL_Qr)==1){
          fm <- do.call(SL_Qr, args=list(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                                             family=stats::gaussian(),
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
  if(!is.null(glm_Qr)){
    Qrn <- mapply(a=a_0, g=gn, Q=Qn, SIMPLIFY=FALSE, FUN=function(a,g,Q){
      if(length(unique(g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        glm_Qr <- "1"
      }
      fm <- stats::glm(stats::as.formula(paste0("Qrn ~",glm_Qr)), 
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
#' @param SL_gr A vector of characters or a list describing the Super Learner library to be used 
#' for the second reduced-dimension regression.
#' @param glm_gr A character describing a formula to be used in the call to \code{glm} for the second reduced-dimension regression. Ignored
#' if \code{SL_gr!=NULL}.
#' @param reduction A character equal to \code{"univariate"} for a univariate misspecification correction or \code{"bivariate"}
#' for the bivariate version. 
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param a_0 A list of fixed treatment values 
#' @param returnModels A boolean indicating whether to return model fits for the outcome regression, propensity score,
#' and reduced-dimension regressions.
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimategrn <- function(Y, A, W, Qn, gn, SL_gr, tolg, glm_gr, a_0, reduction,returnModels){
  if(is.null(SL_gr) & is.null(glm_gr)) stop("Specify Super Learner library or GLM formula for gr")
  if(!is.null(SL_gr) & !is.null(glm_gr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_gr <- NULL
  }
  # Super Learner
  if(!is.null(SL_gr)){
    grn <- mapply(a=a_0,Q=Qn,g=gn,SIMPLIFY=FALSE, FUN=function(a,Q,g){
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
        if(length(SL_gr)>1){
          if(reduction=="univariate"){
            fm1 <- SuperLearner::SuperLearner(Y=(as.numeric(A==a)-g)/g, X=data.frame(Qn=Q), 
                                family=stats::gaussian(),SL.library=SL_gr,  
                                method="method.CC_LS")
            fm2 <- SuperLearner::SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q), 
                                family=stats::binomial(),SL.library=SL_gr,
                                method="method.CC_nloglik")
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
                                family=stats::binomial(),SL.library=SL_gr,
                                method = "method.CC_nloglik")
            if(!all(fm2$coef==0)){
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q, gn = g), onlySL=TRUE)[[1]]
            }else{
              grn2 <- stats::predict(fm2, newdata=data.frame(Qn=Q, gn = g), onlySL=FALSE)[[2]][,which(fm2$cvRisk == min(fm2$cvRisk,na.rm=TRUE))]
            }
            grn2[grn2 < tolg] <- tolg
            fm1 <- NULL; grn1 <- NULL
          }
        }else if(length(SL_gr)==1){
          if(reduction=="univariate"){
            fm1 <- do.call(SL_gr, 
                            args=list(Y=(as.numeric(A==a)-g)/g,X=data.frame(Qn=Q),
                                      obsWeights=rep(1, length(A)),
                                      newX=data.frame(Qn=Q), family=stats::gaussian()))
            grn1 <- predict(object=fm1$fit, newdata=data.frame(Qn=Q))
            fm2 <- do.call(SL_gr, args=list(Y=as.numeric(A==a), X=data.frame(Qn=Q),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q), family=stats::binomial()))
            grn2 <- predict(object=fm2$fit, newdata=data.frame(Qn=Q))
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm2 <- do.call(SL_gr, args=list(Y=as.numeric(A==a),X=data.frame(Qn=Q,gn=g),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q, gn=g), family=stats::binomial()))
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
  if(!is.null(glm_gr)){
    grn <- mapply(a=a_0,Q=Qn,g=gn,SIMPLIFY=FALSE, FUN=function(a,Q,g){
      if(length(unique(Qn[[1]]))==1){
        glm_gr <- "1"
      }
      fm2 <- stats::glm(stats::as.formula(paste0("A~",glm_gr)), family="binomial",
           data=data.frame(A=A, Qn=Q))
      grn2 <- stats::predict(fm2, type="response")
      grn2[grn2 < tolg] <- tolg
      if(reduction == "univariate"){
        fm1 <- stats::glm(stats::as.formula(paste0("grn1~",glm_gr)), family="gaussian",
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

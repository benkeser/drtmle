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
#' @param a_0 A vector of fixed treatment values
#' @param validRows A \code{list} of length \code{cvFolds} containing the row indexes
#' of observations to include in validation fold. 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 

estimateG <- function(A, W, SL_g, glm_g, a_0, tolg, validRows = NULL,
                      verbose=FALSE, returnModels=FALSE){
  if(is.null(SL_g) & is.null(glm_g)) stop("Specify Super Learner library or GLM formula for g")
  if(!is.null(SL_g) & !is.null(glm_g)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_g <- NULL
  }

  # subset data into training and validation sets
  if(!is.null(validRows)){
    trainA <- A[-validRows]
    trainW <- W[-validRows,,drop=FALSE]
    validW <- W[validRows,,drop=FALSE]
    validA <- A[validRows]
  }else{
    trainA <- A
    trainW <- W
    validW <- W
    validA <- A
  }
  # if a super learner library is specified, 
  # fit the super learner
  if(!is.null(SL_g)){
    # if the library is of length > 1, then call SuperLearner
    if(length(SL_g)>1 | is.list(SL_g)){
      # if there are only two unique values of A, then only need one fit
      if(length(a_0)==length(unique(A)) & length(unique(A))==2){
        fm <- SuperLearner::SuperLearner(Y=as.numeric(trainA==a_0[1]), 
                                         X=trainW, 
                                         family=stats::binomial(),
                                         SL.library=SL_g,verbose=verbose,
                                         method = "method.CC_nloglik")
        pred <- stats::predict(fm, newdata = validW, onlySL=TRUE)[[1]]
        gn <- vector(mode="list",length=2)
        gn[[1]] <- pred; gn[[2]] <- 1-pred
      # if there are more than two unique values of A, then we need more
      # than one call to super learner
      }else{
        a_ct <- 0
        gn <- vector(mode = "list", length = length(a_0))
        fm <- vector(mode = "list", length = length(a_0) - 1)
        for(a in a_0[1:(length(a_0)-1)]){
          # determine who to include in the regression for this outcome
          if(a_ct == 0){  
            include <- rep(TRUE, length(trainA))
          }else{
            include <- !(trainA %in% a_0[1:a_ct])          
          }
          # fit super learner
          tmp_fm <- SuperLearner::SuperLearner(Y=as.numeric(trainA[include]==a), 
                                               X=trainW[include,,drop=FALSE], 
                                               newX = validW,
                                               family=stats::binomial(),
                                               SL.library=SL_g,
                                               verbose=verbose)
          # get predictions
          tmp_pred <- tmp_fm$SL.pred
          if(a_ct != 0){
            gn[[a_ct + 1]] <- tmp_pred * 
                    Reduce("*",lapply(gn[1:a_ct], function(x){ 1 - x }))
          }else{
            gn[[a_ct + 1]] <- tmp_pred
          }
          fm[[a_ct + 1]] <- tmp_fm
          a_ct <- a_ct + 1
        }
        # add in final predictions
        gn[[a_ct+1]] <- 1 - Reduce("+",gn[1:a_ct])
      }
    }else if(!is.list(SL_g) & length(SL_g)==1){
      if(length(a_0)==length(unique(A)) & length(unique(A))==2){
        gn <- vector(mode = "list", length = 2)
        fm <- do.call(SL_g, args=list(
                        Y=as.numeric(trainA==a_0[1]), X=trainW, 
                        newX=validW, obsWeights=rep(1,length(trainA)),
                        family=stats::binomial()
                        ))
        pred <- fm$pred
        gn[[1]] <- pred; gn[[2]] <- 1 - pred
      }else{
        a_ct <- 0
        gn <- vector(mode = "list", length = length(a_0))
        fm <- vector(mode = "list", length = length(a_0) - 1)
        for(a in a_0[1:(length(a_0)-1)]){
          # determine who to include in the regression for this outcome
          if(a_ct == 0){  
            include <- rep(TRUE, length(trainA))
          }else{
            include <- !(trainA %in% a_0[1:a_ct])          
          }
          # fit super learner
          tmp_fm <- do.call(SL_g, 
                            args=list(Y=as.numeric(trainA[include]==a),
                                      X=trainW[include,,drop=FALSE],
                                      newX=validW, obsWeights=rep(1,length(A)),
                                      family=stats::binomial()))
          # get predictions
          tmp_pred <- tmp_fm$pred
          if(a_ct != 0){
            gn[[a_ct + 1]] <- tmp_pred * 
                    Reduce("*",lapply(gn[1:a_ct], function(x){ 1 - x }))
          }else{
            gn[[a_ct + 1]] <- tmp_pred
          }
          fm[[a_ct + 1]] <- tmp_fm
          a_ct <- a_ct + 1
        }
        # add in final predictions
        gn[[a_ct+1]] <- 1 - Reduce("+",gn[1:a_ct])
      }
    }
  }
  #----------------------------------------------------------------------
  # GLM
  #----------------------------------------------------------------------
  if(!is.null(glm_g)){
    if(length(a_0)==length(unique(A)) & length(unique(A))==2){
      thisDat <- data.frame(A = as.numeric(trainA==a_0[1]), trainW)
      fm <- stats::glm(stats::as.formula(paste0("A~",glm_g)), data=thisDat, 
                         family=stats::binomial())
      pred <- stats::predict(fm, type = "response")
      gn <- vector(mode="list",length=2)
      gn[[1]] <- pred; gn[[2]] <- 1-pred
    }else{
      a_ct <- 0
      gn <- vector(mode = "list", length = length(a_0))
      fm <- vector(mode = "list", length = length(a_0) - 1)
      for(a in a_0[1:(length(a_0)-1)]){
        # determine who to include in the regression for this outcome
        if(a_ct == 0){  
          include <- rep(TRUE, length(A))
        }else{
          include <- !(A %in% a_0[1:a_ct])          
        }
        # fit super learner
        thisDat <- data.frame(as.numeric(trainA[include]==a), 
                              trainW[include,,drop=FALSE])
        colnames(thisDat) <- c("A",colnames(W))
        tmp_fm <- stats::glm(stats::as.formula(paste0("A~",glm_g)), data=thisDat, 
                         family=stats::binomial())
        tmp_pred <- stats::predict(tmp_fm, newdata = data.frame(A=validA, validW),
                                   type="response")
        # get predictions
        if(a_ct != 0){
          gn[[a_ct + 1]] <- tmp_pred * 
                  Reduce("*",lapply(gn[1:a_ct], function(x){ 1 - x }))
        }else{
          gn[[a_ct + 1]] <- tmp_pred
        }
        fm[[a_ct + 1]] <- tmp_fm
        a_ct <- a_ct + 1
      } # end for loop over treatment levels
      # add in final predictions
      gn[[a_ct+1]] <- 1 - Reduce("+",gn[1:a_ct])
    } # end multi-level treatment if
  } # end glm_g if

  # truncate too small predictions
  gn <- lapply(gn, function(g){ g[g < tolg] <- tolg; g })

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
#' @param validRows A \code{list} of length \code{cvFolds} containing the row indexes
#' of observations to include in validation fold.
#' @param ... Additional arguments (not currently used) 
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 


estimateQ <- function(Y,A,W,SL_Q,glm_Q,a_0,stratify,family,
                      verbose=FALSE,returnModels=FALSE,
                      validRows=NULL,...){
  if(is.null(SL_Q) & is.null(glm_Q)) stop("Specify Super Learner library or GLM formula for Q")
  if(!is.null(SL_Q) & !is.null(glm_Q)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_Q <- NULL
  }
  # subset data into training and validation sets
  if(!is.null(validRows)){
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows,,drop=FALSE]
    validW <- W[validRows,,drop=FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
  }else{
    trainA <- A
    trainW <- W
    trainY <- Y
    validW <- W
    validA <- A
    validY <- Y
  }
  # Super Learner
  if(!is.null(SL_Q)){
    if(!stratify){
      if(length(SL_Q)>1 | is.list(SL_Q)){
        fm <- SuperLearner::SuperLearner(
                Y=trainY, X=data.frame(A=trainA,trainW),
                verbose=verbose,family=family,
                SL.library=SL_Q,
                method = ifelse(family$family=="binomial",
                                "method.CC_nloglik","method.CC_LS"),
                ...)
      
        Qn <- alply(a_0,1,function(x){
          stats::predict(fm, newdata=data.frame(A=x,validW), onlySL=TRUE)[[1]]
        })
      }else if(length(SL_Q)==1){
        fm <- do.call(SL_Q, args=list(Y=trainY, X=data.frame(A = trainA,trainW), 
                                      verbose=verbose, 
                                      newX=data.frame(A=validA,validW),
                                      obsWeights=rep(1,length(trainA)),
                                      family=family))
        Qn <- alply(a_0,1,function(x){
          stats::predict(object=fm$fit, newdata=data.frame(A=x, validW))
        })
      }
   }else{
     if(length(SL_Q)>1 | is.list(SL_Q)){
      tmp <- plyr::alply(a_0,1,function(x){
        fm <- SuperLearner::SuperLearner(
                Y=trainY[trainA==x], 
                X=trainW[trainA==x,,drop=FALSE],
                newX = validW, 
                verbose=verbose,family=family,
                SL.library=SL_Q,
                ifelse(family$family=="binomial",
                       "method.CC_nloglik","method.CC_LS"))
        list(est = fm$SL.predict,fm = fm)
      })
      Qn <- lapply(tmp,"[[",1)
      fm <- lapply(tmp,"[[",2)
     }else if(length(SL_Q)==1){
       tmp <- plyr::alply(a_0,1,function(x){
         fm <- do.call(SL_Q, args=list(Y=trainY[trainA==x], 
                                       X=trainW[trainA==x,,drop=FALSE], 
                                       newX=validW,
                                       verbose=verbose,
                                       obsWeights=rep(1,sum(trainA==x)),
                                       family=family))
         list(est = fm$pred, fm = fm)
       })
       Qn <- lapply(tmp,"[[",1)
       fm <- lapply(tmp,"[[",2)
     }
    }
  }
  
  # GLM
  if(!is.null(glm_Q)){
    thisDat <- data.frame(Y=trainY, A=trainA, trainW)
    if(!stratify){
      fm <- stats::glm(stats::as.formula(paste0("Y~",glm_Q)), 
                       data=thisDat, family=family)
      Qn <- plyr::alply(matrix(a_0),1,function(x,fm){
        stats::predict(fm, newdata=data.frame(A=x,validW), type="response")
      },fm=fm)
    }else{
      tmp <- plyr::alply(matrix(a_0),1,function(x){
        fm <- stats::glm(stats::as.formula(paste0("Y~",glm_Q)), 
                         data=thisDat[trainA==x,,drop=FALSE], family=family)
        list(est = stats::predict(fm, newdata=data.frame(A=x,validW), 
                                  type="response"),fm = fm)
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
#' @param validRows A \code{list} of length \code{cvFolds} containing the row indexes
#' of observations to include in validation fold.
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimateQrn  <- function(Y, A, W, Qn, gn, glm_Qr, SL_Qr, a_0, returnModels,
                         validRows = NULL){
  # if estimateQrn is called in islptw, then Qn will enter as 
  # NULL. Here we fill its value to 0 so that we estimate
  # the correct nuisance parameter for islptw
  if(is.null(Qn)){
    Qn <- vector(mode = "list", length = length(a_0))
    for(i in 1:length(a_0)){ 
      Qn[[i]] <- rep(0, length(Y)) 
    }
  }

  # subset data into training and validation sets
  if(!is.null(validRows)){
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows,,drop=FALSE]
    train_gn <- lapply(gn, "[", -validRows)
    train_Qn <- lapply(Qn, "[", -validRows)
    validW <- W[validRows,,drop=FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
    valid_gn <- lapply(gn, "[", validRows)
    valid_Qn <- lapply(Qn, "[", validRows)
  }else{
    trainA <- A
    trainW <- W
    trainY <- Y
    validW <- W
    validA <- A
    validY <- Y
    train_gn <- gn; train_Qn <- Qn
    valid_gn <- gn; valid_Qn <- Qn
  }

  if(is.null(SL_Qr) & is.null(glm_Qr)) stop("Specify Super Learner library or GLM formula for Qr")
  if(!is.null(SL_Qr) & !is.null(glm_Qr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_Qr <- NULL
  }
  # Super Learner
  if(!is.null(SL_Qr)){
    Qrn <- mapply(a=a_0, train_g=train_gn, train_Q=train_Qn, 
                  valid_g=valid_gn, valid_Q=valid_Qn,
                  SIMPLIFY=FALSE, 
                  FUN=function(a,train_g,train_Q,valid_g,valid_Q){
      if(length(unique(train_g))==1){
        warning(paste0("Only one unique value of gn",a,". Using empirical average as Qr estimate."))
        m1 <- mean((trainY-train_Q)[trainA==a])
        est <- rep(m1, length(validY))
        fm <- list(object = m1)
        class(fm) <- "SL.mean"
      }else{
        if(length(SL_Qr)>1){
          suppressWarnings(
          fm <- SuperLearner::SuperLearner(
                  Y=(trainY-train_Q)[trainA==a], 
                  X=data.frame(gn=train_g[trainA==a]),
                  newX=data.frame(gn=valid_g),
                  family=stats::gaussian(),
                  SL.library=SL_Qr, 
                  method="method.CC_LS"))
          # if all weights = 0, use discrete SL
          if(!all(fm$coef==0)){
            est <- fm$SL.predict
          }else{
            est <- fm$library.predict[,which.min(fm$cvRisk)]
          }
        }else if(length(SL_Qr)==1){
          fm <- do.call(SL_Qr, args=list(
                    Y=(trainY-train_Q)[trainA==a], 
                    X=data.frame(gn=train_g[trainA==a]),
                    family=stats::gaussian(),
                    newX=data.frame(gn=valid_g),
                    obsWeights=rep(1, length(trainY[trainA==a]))))
          est <- fm$pred
        }
      }
      out <- list(est = est, fm = NULL)
      if(returnModels) out$fm <- fm
      return(out)
    })
  }
  
  # GLM
  if(!is.null(glm_Qr)){
    Qrn <- mapply(a=a_0, train_g=train_gn, train_Q=train_Qn, 
                  valid_g = valid_gn, valid_Q = valid_Qn,  
                  SIMPLIFY=FALSE, FUN=function(a,train_g,train_Q,valid_g,valid_Q){
      if(length(unique(train_g))==1){
        warning(paste0("Only one unique value of gn",a,
                       ". Using empirical average as Qr estimate."))
        glm_Qr <- "1"
      }
      fm <- stats::glm(stats::as.formula(paste0("Qrn ~",glm_Qr)), 
                data=data.frame(Qrn=(trainY-train_Q)[trainA==a], 
                                gn=train_g[trainA==a]),
                family="gaussian")
      est <- stats::predict(fm, newdata=data.frame(gn=valid_g),type="response")
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
#' @param validRows A \code{list} of length \code{cvFolds} containing the row indexes
#' of observations to include in validation fold.
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimategrn <- function(Y, A, W, Qn, gn, SL_gr, tolg, 
                        glm_gr, a_0, reduction,returnModels,
                        validRows){
  if(!is.null(validRows)){
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows,,drop=FALSE]
    train_gn <- lapply(gn, "[", -validRows)
    train_Qn <- lapply(Qn, "[", -validRows)
    validW <- W[validRows,,drop=FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
    valid_gn <- lapply(gn, "[", validRows)
    valid_Qn <- lapply(Qn, "[", validRows)
  }else{
    trainA <- A
    trainW <- W
    trainY <- Y
    validW <- W
    validA <- A
    validY <- Y
    train_gn <- gn; train_Qn <- Qn
    valid_gn <- gn; valid_Qn <- Qn
  }

  if(is.null(SL_gr) & is.null(glm_gr)){
    stop("Specify Super Learner library or GLM formula for gr")
  }
  if(!is.null(SL_gr) & !is.null(glm_gr)){
    warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
    glm_gr <- NULL
  }
  # Super Learner
  if(!is.null(SL_gr)){
    grn <- mapply(a=a_0,train_Q=train_Qn,train_g=train_gn,
                  valid_Q = valid_Qn, valid_g = valid_gn,
                  SIMPLIFY=FALSE, 
                  FUN=function(a,train_Q,train_g,valid_Q,valid_g){
      if(length(unique(train_Q))==1){
        warning("Only one unique value of Qn. Proceeding with empirical mean for grn")
        if(reduction=="univariate"){
          m1 <- mean((as.numeric(trainA==a)-train_g)/train_g)
          grn1 <- rep(m1, length(validY))
          m2 <- mean(as.numeric(trainA==a))
          grn2 <- rep(m2, length(validY))
          grn2[grn2 < tolg] <- tolg
          fm1 <- list(object = m1)
          class(fm1) <- "SL.mean"
          fm2 <- list(object = m2)
          class(fm1) <- "SL.mean"
        }else if(reduction=="bivariate"){
          m2 <- mean(as.numeric(trainA==a))
          grn2 <- rep(m2, length(validY))
          grn2[grn2 < tolg] <- tolg
          fm2 <- list(object = m2)
          class(fm2) <- "SL.mean"
          fm1 <- NULL; grn1 <- rep(NA,length(validY))
        }
      }else{
        if(length(SL_gr)>1){
          if(reduction=="univariate"){
            fm1 <- SuperLearner::SuperLearner(
                    Y=(as.numeric(trainA==a)-train_g)/train_g, 
                    X=data.frame(Qn=train_Q), 
                    newX = data.frame(Qn = valid_Q),
                    family=stats::gaussian(),
                    SL.library=SL_gr,  
                    method="method.CC_LS")
            fm2 <- SuperLearner::SuperLearner(
                    Y=as.numeric(trainA==a), 
                    X=data.frame(Qn=train_Q), 
                    newX = data.frame(Qn = valid_Q),
                    family=stats::binomial(),SL.library=SL_gr,
                    method="method.CC_nloglik")
            if(!all(fm1$coef==0)){
              grn1 <- fm1$SL.predict
            }else{
              grn1 <- fm1$library.predict[,which.min(fm1$cvRisk)]
            }

            if(!all(fm2$coef==0)){
              grn2 <- fm2$SL.predict
            }else{
              grn2 <- fm2$library.predict[,which.min(fm2$cvRisk)]
            }
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm2 <- SuperLearner::SuperLearner(
                      Y=as.numeric(trainA==a), 
                      X=data.frame(Qn=train_Q, gn=train_g), 
                      newX = data.frame(Qn = valid_Q, gn = valid_g),
                      family=stats::binomial(),
                      SL.library=SL_gr,
                      method = "method.CC_nloglik")
            if(!all(fm2$coef==0)){
              grn2 <- fm2$SL.predict
            }else{
              grn2 <- fm2$library.predict[,which.min(fm2$cvRisk)]
            }
            grn2[grn2 < tolg] <- tolg
            fm1 <- NULL; grn1 <- rep(NA,length(validY))
          }
        }else if(length(SL_gr)==1){
          if(reduction=="univariate"){
            fm1 <- do.call(SL_gr, 
                            args=list(Y=(as.numeric(trainA==a)-train_g)/train_g,
                                      X=data.frame(Qn=train_Q),
                                      obsWeights=rep(1, length(trainA)),
                                      newX=data.frame(Qn=valid_Q), 
                                      family=stats::gaussian()))
            grn1 <- fm1$pred
            fm2 <- do.call(SL_gr, args=list(Y=as.numeric(trainA==a), 
                                            X=data.frame(Qn=train_Q),
                                            obsWeights=rep(1, length(trainA)),
                                            newX=data.frame(Qn=valid_Q), 
                                            family=stats::binomial()))
            grn2 <- fm2$pred
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            fm2 <- do.call(SL_gr, args=list(Y=as.numeric(trainA==a),
                                            X=data.frame(Qn=train_Q,gn=train_g),
                                            obsWeights=rep(1, length(trainA)),
                                            newX=data.frame(Qn=valid_Q, gn=valid_g), 
                                            family=stats::binomial()))
            grn2 <- fm2$pred
            grn2[grn2 < tolg] <- tolg
            fm1 <- NULL; grn1 <- rep(NA,length(validY))
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
    grn <- mapply(a=a_0,train_Q = train_Qn, train_g=train_gn,
                  valid_Q = valid_Qn, valid_g = valid_gn,
                  SIMPLIFY=FALSE, FUN=function(a,train_Q,train_g,valid_Q,valid_g){
      if(length(unique(train_Q))==1){
        glm_gr <- "1"
      }
      if(reduction == "univariate"){
        fm1 <- stats::glm(stats::as.formula(paste0("grn1~",glm_gr)), family="gaussian",
           data=data.frame(grn1=(as.numeric(trainA==a)-train_g)/train_g, Qn=train_Q))
        grn1 <- stats::predict(fm1, newdata = data.frame(
                  grn1 = (as.numeric(validA==a)-valid_g)/valid_g, 
                  Qn=valid_Q),type="response")
        fm2 <- stats::glm(stats::as.formula(paste0("A~",glm_gr)), family="binomial",
                          data=data.frame(A=as.numeric(trainA==a), Qn=train_Q))
        grn2 <- stats::predict(fm2, newdata = data.frame(A=validA, Qn = valid_Q), 
                               type="response")
      }else if(reduction == "bivariate"){
        fm1 <- NULL; grn1 <- rep(NA,length(validY))
        fm2 <- stats::glm(stats::as.formula(paste0("A~",glm_gr)), family="binomial",
                          data=data.frame(A=as.numeric(trainA==a), Qn=train_Q, gn = train_g))
        grn2 <- stats::predict(fm2, newdata = data.frame(A=validA, Qn = valid_Q, gn = valid_g), 
                               type="response")
      }
      grn2[grn2 < tolg] <- tolg
      out <- list(grn1 = grn1, grn2 = grn2, fm1 = NULL, fm2 = NULL)
      if(returnModels){
        out$fm1 <- fm1; out$fm2 <- fm2
      }
      return(out)
    })
  }
  tmp1 <- lapply(grn,function(x){ data.frame(grn1 = x$grn1, grn2 = x$grn2) })
  tmp2 <- lapply(grn,function(x){ list(fm1 = x$fm1, fm2 = x$fm2) })
  return(list(est = tmp1, fm = tmp2))
}

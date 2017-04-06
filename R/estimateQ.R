
#' estimateQ 
#' @importsFrom SuperLearner SuperLearner trimLogit
#' @importsFrom stats predict glm as.formula
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
#' 
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
      Qn <- plyr::alply(a0,1,function(x){
        fm <- SuperLearner::SuperLearner(Y=Y[A==x], X=W[A==x,], verbose=verbose,family=family,SL.library=libraryQ)
        predict(fm, newdata=data.frame(A=x,W), onlySL=TRUE)[[1]]
      })
     }else if(length(libraryQ)==1){
       Qn <- plyr::alply(a0,1,function(x){
         fm <- do.call(libraryQ, args=list(Y=Y[A==x], X=W[A==x,], newX=W[A==x,],verbose=verbose,obsWeights=rep(1,sum(A==x)),
                                           family=family))
         stats::predict(object=fm$fit, newdata=data.frame(W))
       })
     }
    }
  }
  
  # GLM
  if(!is.null(glmQ)){
    thisDat <- data.frame(Y=Y, A=A, W=W); colnames(thisDat) <- c("Y","A",colnames(W))
    if(!stratify){
      fm <- stats::glm(stats::as.formula(paste0("Y~",glmQ)), data=thisDat, family=family)
      Qn <- plyr::alply(matrix(a0),1,function(x,fm){
        stats::predict(fm, newdata=data.frame(A=x,W), type="response")
      },fm=fm)
    }else{
      Qn <- plyr::alply(matrix(a0),1,function(x){
        fm <- stats::glm(stats::as.formula(paste0("Y~",glmQ)), data=thisDat[A==x,], family=family)
        stats::predict(fm, newdata=data.frame(A=x,W), type="response")
      })     
    }
  }
  if(returnModels) list(Qn,fm)
  else Qn
}


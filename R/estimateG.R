
#' estimateG
#' 
#' Function to estimate propensity score
#' @importsFrom SuperLearner SuperLearner trimLogit
#' @importsFrom stats predict glm as.formula
#' @importsFrom plyr alply
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
#' Estimates the reduced dimension regressions necessary for the additional 
#' fluctuations. 

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
        fm <- SuperLearner::SuperLearner(Y=as.numeric(A==a0[1]), X=W, family="binomial",SL.library=libraryg,verbose=verbose)
        pred <- stats::predict(fm, onlySL=TRUE)[[1]]
        pred[pred < tolg] <- tolg
        gn <- vector(mode="list",length=2)
        gn[[1]] <- pred; gn[[2]] <- 1-pred
      }else{
        gn <- plyr::alply(a0, 1, function(x,A,W,libraryg){
          fm <- SuperLearner::SuperLearner(Y=as.numeric(A==x), X=W, family="binomial",SL.library=libraryg,verbose=verbose)
          pred <- stats::predict(fm, onlySL=TRUE)[[1]]
          pred[pred < tolg] <- tolg
          pred
        }, A=A, W=W, libraryg=libraryg)
      }
    }else if(!is.list(libraryg) & length(libraryg)==1){
      gn <- plyr::alply(a0, 1, function(x){
        fm <- do.call(libraryg, args=list(Y=as.numeric(A==x), X=W, newX=W, obsWeights=rep(1,length(A)),family=data.frame(family="binomial")))
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
    })
  }
  if(returnModels){
    return(list(gn, fm))
  }else{
    return(gn)
  }
}

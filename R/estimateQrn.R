
#' estimateQrn 
#' 
#' Estimates the reduced dimension regressions necessary for the  
#' fluctuations of g
#' @importsFrom SuperLearner SuperLearner trimLogit
#' @importsFrom stats predict glm as.formula
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
#' Estimates the reduced dimension regressions necessary for the additional 
#' fluctuations. 
#' 
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
          fm <- SuperLearner::SuperLearner(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),
                             family=gaussian(),SL.library=libraryQr, method="method.NNLS2")
          )
          # if all weights = 0, use discrete SL
          if(!all(fm$coef==0)){
            stats::predict(fm, newdata=data.frame(gn=g), onlySL=TRUE)[[1]]
          }else{
            stats::predict(fm, newdata=data.frame(gn=g), onlySL=FALSE)[[2]][,which(fm$cvRisk == min(fm$cvRisk, na.rm = TRUE))]
          }
        }else if(length(libraryQr)==1){
          obj <- do.call(libraryQr, args=list(Y=(Y-Q)[A==a], X=data.frame(gn=g[A==a]),family=gaussian(),
                                              newX=data.frame(gn=g[A==a]),
                                              obsWeights=rep(1, length(Y[A==a]))))
          pred <- stats::predict(object=obj$fit, newdata=data.frame(gn=g))
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
        fm <- stats::glm(stats::as.formula(paste0("Qrn ~",glmQr)), 
                  data=data.frame(Qrn=(Y-Q)[A==a], gn=g[A==a]),
                  family="gaussian")
        stats::predict(fm, newdata=data.frame(gn=g),type="response")
      }
    })
  }
  Qrn 
}

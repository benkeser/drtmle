
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
#' 
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#' 
#' 
#' 
#' 

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
            fm1 <- SuperLearner::SuperLearner(Y=(as.numeric(A==a)-g)/g, X=data.frame(Qn=Q), 
                                family="gaussian",SL.library=librarygr,  method="method.NNLS2")
            fm2 <- SuperLearner::SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q), 
                                family=data.frame(family="binomial"),SL.library=librarygr)
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
            fm1 <- SuperLearner::SuperLearner(Y=as.numeric(A==a), X=data.frame(Qn=Q, gn=g), 
                                family=data.frame(family="binomial"),SL.library=librarygr)
            if(!all(fm1$coef==0)){
              grn <- stats::predict(fm1, newdata=data.frame(Qn=Q), onlySL=TRUE)[[1]]
            }else{
              grn <- stats::predict(fm1, newdata=data.frame(Qn=Q), onlySL=FALSE)[[2]][,which(fm1$cvRisk == min(fm1$cvRisk,na.rm=TRUE))]
            }
            grn[grn < tolg] <- tolg
          }
        }else if(length(librarygr)==1){
          if(reduction=="univariate"){
            obj1 <- do.call(librarygr, 
                            args=list(Y=(as.numeric(A==a)-g)/g,X=data.frame(Qn=Q),
                                      obsWeights=rep(1, length(A)),
                                      newX=data.frame(Qn=Q), family=data.frame(family="gaussian")))
            grn1 <- predict(object=obj1$fit, newdata=data.frame(Qn=Q))
            obj2 <- do.call(librarygr, args=list(Y=as.numeric(A==a), X=data.frame(Qn=Q),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q), family=data.frame(family="binomial")))
            grn2 <- predict(object=obj2$fit, newdata=data.frame(Qn=Q))
            grn2[grn2 < tolg] <- tolg
          }else if(reduction=="bivariate"){
            obj <- do.call(librarygr, args=list(Y=as.numeric(A==a),X=data.frame(Qn=Q,gn=g),obsWeights=rep(1, length(A)),
                                                 newX=data.frame(Qn=Q, gn=g), family=data.frame(family="binomial")))
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
        fm1 <- stats::glm(stats::as.formula(paste0("grn1~",glmgr)), family="gaussian",
                   data=data.frame(grn1=(as.numeric(A==a)-g)/g, Qn=Q))
        fm2 <- stats::glm(stats::as.formula(paste0("A~",glmgr)), family="binomial",
                   data=data.frame(A=A, Qn=Q))
        grn1 <- stats::predict(fm1, type="response")
        grn2 <- stats::predict(fm2, type="response")
        grn2[grn2 < tolg] <- tolg
      }
      list(grn1,grn2)
    })
  }
  
  grn
}

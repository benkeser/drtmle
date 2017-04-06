#' fluctuateQ2 
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param grn A list of reduced-dimension regression estimates evaluated on observed data
#' @param coefTol A tolerance level on the magnitude of the coefficient that flags the
#' result as potentially the result of numeric instability.
#' @param reduction A character indicating what reduced dimension regression was used. 
#' @param a0 A list of fixed treatment values 
#' 
#' @importsFrom SuperLearner trimLogit
#' @importsFrom stats predict glm
#' 
#' 
#' Function called internally by drtmle to perform the second fluctuation 
#' of the initial estimator of Q (i.e., solves the new estimating eqn that results
#' from misspecification of g)

fluctuateQ2 <- function(Y,A,W,Qn,gn,grn,a0,family,reduction,coefTol=1e5){
  QnStar <- mapply(a=a0,Q=Qn,g=gn,gr=grn,FUN=function(a, Q, g, gr){
    l <- min(Y); u <- max(Y)
    Yscale <- (Y-l)/(u-l)
    off <- SuperLearner::trimLogit((Q-l)/(u-l))
    if(reduction=="univariate") H2 <- as.numeric(A==a)/gr[[2]] * gr[[1]]
    if(reduction=="bivariate") H2 <- as.numeric(A==a)/gr[[1]] * (gr[[1]]-g)/g
    suppressWarnings(
      fm <- stats::glm(Yscale ~ -1 + offset(off) + H2, start=c(0),
                data=data.frame(Y=Y, off=off, H2=H2), family="binomial")
    )
    if(!fm$converged | abs(max(fm$coef)) > coefTol){
      # if it doesn't converge, try with no starting values
      suppressWarnings(
        fm <- stats::glm(Yscale ~ -1 + offset(off) + H2, 
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
      Qnstar <- stats::predict(fm,type="response",newdata=data.frame(off=off, H2=1/gr[[2]] * gr[[1]]))*(u-l) + l
      list(est=Qnstar, 
           eps=fm$coef)
    }else if(reduction=="bivariate"){
      Qnstar <- stats::predict(fm,type="response",newdata=data.frame(off=off, H2=1/gr[[1]] * (gr[[1]]-g)/g))*(u-l) + l
      list(est=Qnstar,
           eps=fm$coef)
    }
  }, SIMPLIFY=FALSE)
  QnStar
}

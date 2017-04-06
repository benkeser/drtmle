#' fluctuateQ1
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param a0 A list of fixed treatment values 
#' 
#' @importsFrom SuperLearner trimLogit
#' @importsFrom stats predict glm
#' 
#' Function called internally by drtmle to perform the first fluctuation 
#' of the initial estimator of Q (i.e., solves the original EIF estimating eqn)

fluctuateQ1 <- function(Y,A,W, Qn, gn, a0){
  QnStar <- mapply(a=a0,Q=Qn,g=gn,FUN=function(x, a, Q, g){
      l <- min(Y); u <- max(Y)
      Yscale <- (Y-l)/(u-l)
      off <- SuperLearner::trimLogit((Q-l)/(u-l))
      H1 <- as.numeric(A==a)/g
      suppressWarnings(
      fm <- stats::glm(Yscale ~ -1 + offset(off) + H1, start=0,
                data=data.frame(Y=Y, off=off, H1=H1), family="binomial")
      )
      Qnstar <- stats::predict(fm,type="response",
                       newdata=data.frame(off=off, H1=1/g))*(u-l) + l
      list(est=Qnstar,eps=fm$coef)
    }, SIMPLIFY=FALSE)
  QnStar
}


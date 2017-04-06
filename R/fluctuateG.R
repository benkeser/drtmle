#' fluctuateG
#' 
#' Function called internally by drtmle to perform the fluctuation 
#' of the initial estimator of g (i.e., solves the new estimating eqn that results
#' from misspecification of Q)
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param Qn A list of outcome regression estimates evaluated on observed data
#' @param gn A list of propensity regression estimates evaluated on observed data
#' @param Qrn A list of reduced-dimension regression estimates evaluated on observed data
#' @param coefTol A tolerance level on the magnitude of the coefficient that flags the
#' result as potentially the result of numeric instability.
#' @param tolg The lower bound on propensity score estimates
#' @param a0 A list of fixed treatment values 
#' 
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#' 

fluctuateG <- function(Y, A, W, a0, Qn, gn, Qrn, tolg, coefTol=1e5){
  gnStar <- mapply(a=a0, Q=Qn, g=gn, Qr=Qrn, FUN=function(x, a, Q, g, Qr){
    H1 <- Qr/g
    off <- SuperLearner::trimLogit(g, tolg)
    thisA <- as.numeric(A==a)
    suppressWarnings(
      fm <- stats::glm(thisA ~ -1 + offset(off) + H1, start=0,
                data=data.frame(thisA=thisA, off=off, H1=H1), family="binomial")
    )
    if(!fm$converged | abs(fm$coef) > coefTol){
      suppressWarnings(
      fm <- stats::glm(thisA ~ -1 + offset(off) + H1,
                data=data.frame(thisA=thisA, off=off, H1=H1), family="binomial")
      )
      if(!fm$converged | abs(fm$coef) > coefTol){
        warning("No sane fluctuation found for G this iteration. Check mean of IC.") 
      }
    }
    pred <- stats::predict(fm, type="response")
    pred[pred < tolg] <- tolg
    list(est=pred, eps=fm$coef)
  }, SIMPLIFY = FALSE)
  gnStar
}


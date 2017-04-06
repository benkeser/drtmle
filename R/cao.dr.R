#' cao.dr 
#' Compute the Cao 2009 estimator
#' 
#' @param R Missing indicator
#' @param Y Outcome
#' @param cov Covariates
#' @param family Character
#' @param nBoot Number of bootstrap resamples 
#' 
#' @export 
#' 
#' @importFrom stats quantile
cao.dr <- function(R,Y,cov, family = "gaussian",nBoot=500){
    # get estimate
    est <- getCaoEst(R=R,Y=Y,cov=cov,family=family)
    # get bootstrap CI
    doOneBoot <- function(){
        samp <- sample(1:length(R), replace = TRUE)
        tmpEst <- getCaoEst(R=R[samp],Y=Y[samp],cov=cov[samp,],
            family=family)
        return(tmpEst)
    }
    estVec <- replicate(nBoot, doOneBoot())
    ci <- stats::quantile(estVec, p = c(0.025, 0.975))
    return(list(est = est, ci = ci))
}

#' SL.npreg
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param family Not used by the function, but needed by \code{SuperLearner}
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param rangeThresh If the the range of the outcomes is smaller than this number, the method
#' returns the empirical average of the outcomes. 
#' @param ... Other arguments (not currently used)
#' 
#' @importsFrom np npregbw npreg
#' @importsFrom stats predict as.formula
#' 
#' SuperLearner wrapper to fit kernel regression
#' @export 

SL.npreg <- function (Y, X, newX, family, obsWeights, 
                      rangeThresh=1e-7, ...) 
{
  if(abs(diff(range(Y))) <= rangeThresh){
    thisMod <- glm(Y ~ 1, data=X)
  }else{
  bw <- np::npregbw(stats::as.formula(paste("Y ~", paste(names(X),collapse="+"))), data=X,
                ftol=0.01, tol=0.01, remin=FALSE)
  
  # fit the kernel regression
  thisMod <- np::npreg(bw)
  }
  
  pred <- stats::predict(thisMod, newdata=newX)
  fit <- list(object = thisMod)
  class(fit) <- "SL.npreg"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' SL.svmLinear.caretMod
#' Uses SL.caretMod to train a linear svm with 15 choices of 
#' tuning parameter
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.svmLinear.caretMod <- function(...,method="svmLinear",tuneLength = 15, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.rpart.caretMod 
#' 
#' Uses SL.caretMod with five-fold CV to train regression tree 
#' with ten choices of tuning parameter.
#' 
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param ... Other arguments passed to SL.caretMod
#' 
#' @export

SL.rpart.caretMod <- function(...,method="rpart",tuneLength = 10, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength, trControl = trControl)
}

#' SL.rf.caretMod
#' 
#' Uses SL.caretMod with five-fold CV to train random forest
#' with ten choices of tuning parameter.
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.rf.caretMod <- function (..., method = "rf", tuneLength = 10, 
                            trControl = caret::trainControl(method = "cv", number = 2)) 
{
    SL.caretMod(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.randomGLM.caretMod
#' Uses SL.caretMod with five-fold CV to train a random GLM
#' with ten choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.randomGLM.caretMod <- function(...,method="randomGLM",tuneLength=10, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}


#' SL.npreg
#' 
#' Kernel regression based on the \href{https://cran.r-project.org/web/packages/np/}{np}
#' package. Uses leave-one-out cross-validation to fit a kernel regression. 
#' See \code{?npreg} for more details.
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param family Not used by the function, but needed for compatibility with \code{SuperLearner}
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param rangeThresh For stability: if the the range of the outcomes is smaller than this number, the method
#' returns the empirical average of the outcomes. 
#' @param ... Other arguments (not currently used)
#' 
#' @importFrom np npregbw npreg
#' @importFrom stats predict as.formula
#' 
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

#' predict.SL.npreg
#' 
#' Method for predicting SL.npreg objects. 
#' 
#' @param object An object of class \code{"SL.npreg"}
#' @param newdata The new data used to obtain predictions 
#' @param ... Other arguments passed to predict
#' 
#' @importFrom stats predict
#' @method predict SL.npreg
#' @export
predict.SL.npreg <- function(object, newdata, ...){
  pred <- stats::predict(object=object$object, newdata=newdata)
  pred
}

#' SL.nnet.caretMod
#' 
#' Uses SL.caretMod to train with five-fold CV to train a neural network
#' with twenty choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.nnet.caretMod <- function(...,method="nnet", tuneLength=20, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.glmnet.caretMod
#' 
#' Uses SL.caretMod with five-fold CV to train a glmnet regression
#' with ten choices of tuning parameter.
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export 
SL.glmnet.caretMod <- function(...,method="glmnet", tuneLength=10, trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}


#' SL.gbm.caretMod
#' 
#' Uses SL.caretMod to train a gbm with 10 choices of tuning parameters and 5-fold CV
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#' 
#' @export

SL.gbm.caretMod <- function (..., method = "gbm", tuneLength = 20, 
                           trControl = caret::trainControl(method = "cv", number = 2)) 
{
    SL.caretMod(..., method = method, tuneLength = tuneLength, trControl = trControl)
}

#' SL.gamSpline.caretMod
#' 
#' Uses SL.caretMod with 10 choices of tuning parameters and 5-fold CV
#' to train a generalized additive model.
#' 
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caretMod}
#'  
#' @export 
SL.gamSpline.caretMod <- function(...,method="gamSpline",tuneLength=10, 
                                  trControl = caret::trainControl(method = "cv", number = 2)){
  SL.caretMod(...,method=method,tuneLength=tuneLength)
}

#' SL.caretMod 
#' 
#' A modification of the SL.caret function from the SuperLearner package that
#' suppresses some output when method = "gbm" or "nnet". See \code{?SL.caret}
#' for more information.
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param family Character indicating family argument (ignored)
#' @param ... Other arguments (not currently used)
#' 
#' @importFrom caret train
#' @importFrom stats predict
#' 
#' @export

SL.caretMod <- function (Y, X, newX, family, obsWeights, method, tuneLength = 10, 
                       trControl = caret::trainControl(method = "cv", number = 2, verboseIter = FALSE), 
                       ...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,verbose=FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl, verbose = FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        # pass trace==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}


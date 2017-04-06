#' SL.rpart.caret1 
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param ... Other arguments (not currently used)
#' 
#' Uses SL.caret1 to train regression tree. 
#' @export

SL.rpart.caret1 <- function(...,method="rpart",tuneLength = 8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength, trControl = trControl)
}
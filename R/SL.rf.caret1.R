#' SL.rf.caret1
#' 
#' Uses SL.caret 1 to train a gbm with 10 choices of tuning parameters and 5-fold CV
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caret1}
#' 
#' Uses SL.caret1 to train random forest
#' @export

SL.rf.caret1 <- function(...,method="rf",tuneLength=8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength, trControl = trControl)
}

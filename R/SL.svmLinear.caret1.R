#' SL.svmLinear.caret1 
#' Uses SL.caret1 to train a linear svm 
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caret1}
#' 
#' @export

SL.svmLinear.caret1 <- function(...,method="svmLinear",tuneLength = 8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

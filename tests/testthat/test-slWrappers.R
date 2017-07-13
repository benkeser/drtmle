library(drtmle)
library(caret)
library(SuperLearner)
library(randomForest)
library(rpart)
library(gam)
library(gbm)
library(kernlab)
library(glmnet)
library(np)
library(randomGLM)

context("Testing super learner wrappers")

test_that("SL.npreg works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.npreg(Y=Y,X=X,newX=X,
	                obsWeights = rep(1,n))
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.npreg")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.svmLinear.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.svmLinear.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.rpart.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.rpart.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.rpart.caretMod works as expected with binary outcome",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rbinom(n,1,plogis(X$X1))
	fit <- SL.rpart.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.rf.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.rf.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.randomGLM.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.randomGLM.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.nnet.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.nnet.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.glmnet.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.glmnet.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10),X2 = rep(0,10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.gbm.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.gbm.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})
test_that("SL.gamSpline.caretMod works as expected",{
	n <- 100
	X <- data.frame(X1 = rnorm(n))
	Y <- rnorm(n,X$X1,1)
	fit <- SL.gamSpline.caretMod(
     Y=Y,X=X,newX=X,tuneLength = 2, 
     trControl = caret::trainControl(method = "cv", number = 2),
     obsWeights = rep(1,n)
    )
	expect_true(is.numeric(fit$pred))
	expect_true(all(!is.na(fit$pred)))
	expect_true(class(fit$fit)=="SL.caret")
	# test predict function
	newX <- data.frame(X1=seq(-1,1,length=10))
	pred <- predict(fit$fit, newdata=newX)
	expect_true(is.numeric(pred))
	expect_true(all(!is.na(pred)))
})


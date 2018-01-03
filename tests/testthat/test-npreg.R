library(drtmle)
library(SuperLearner)
library(np)

context("Testing super learner wrappers")

test_that("SL.npreg works as expected", {
  set.seed(1234)
  n <- 100
  X <- data.frame(X1 = rnorm(n))
  Y <- rnorm(n, X$X1, 1)
  fit <- SL.npreg(
    Y = Y, X = X, newX = X,
    obsWeights = rep(1, n)
  )
  expect_true(is.numeric(fit$pred))
  expect_true(all(!is.na(fit$pred)))
  expect_true(class(fit$fit) == "SL.npreg")
  # test predict function
  newX <- data.frame(X1 = seq(-1, 1, length = 10))
  pred <- predict(fit$fit, newdata = newX)
  expect_true(is.numeric(pred))
  expect_true(all(!is.na(pred)))
})

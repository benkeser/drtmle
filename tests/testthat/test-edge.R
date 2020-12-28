# install.packages("~/Dropbox/R/drtmle",repos=NULL,type="source")

library(drtmle)
library(SuperLearner)
context("Testing edge cases.")

test_that("drtmle executes as expected when only one value of gn", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    cvFolds = 1, maxIter = 2,
    family = gaussian(),
    stratify = TRUE,
    glm_Q = "W1 + W2",
    glm_g = "1",
    glm_Qr = "gn",
    glm_gr = "Qn",
    guard = c("Q", "g"),
    reduction = "univariate"
  )
  expect_true(is.numeric(fit1$gcomp$est))
  expect_true(is.numeric(fit1$tmle$est))
  expect_true(is.numeric(fit1$tmle$est))
  expect_true(is.numeric(fit1$tmle$cov))
  expect_true(is.numeric(fit1$drtmle$est))
  expect_true(is.numeric(fit1$drtmle$cov))
  expect_true(is.numeric(fit1$aiptw$est))
  expect_true(is.numeric(fit1$aiptw$cov))
  expect_true(is.numeric(fit1$aiptw_c$est))
  expect_true(is.numeric(fit1$aiptw_c$cov))
})

test_that("drtmle executes when glm and SL are specified", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  expect_warning(fit1 <- drtmle(
    W = W, A = A, Y = Y,
    cvFolds = 1, maxIter = 2,
    family = gaussian(),
    returnModels = TRUE,
    stratify = FALSE,
    glm_Q = "W1 + W2",
    SL_Q = c("SL.glm", "SL.mean"),
    glm_g = "1",
    SL_g = c("SL.glm", "SL.mean"),
    glm_Qr = "gn",
    SL_Qr = c("SL.glm", "SL.mean"),
    glm_gr = "Qn",
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate"
  ))
  expect_true(is.numeric(fit1$gcomp$est))
  expect_true(is.numeric(fit1$tmle$est))
  expect_true(is.numeric(fit1$tmle$est))
  expect_true(is.numeric(fit1$tmle$cov))
  expect_true(is.numeric(fit1$drtmle$est))
  expect_true(is.numeric(fit1$drtmle$cov))
  expect_true(is.numeric(fit1$aiptw$est))
  expect_true(is.numeric(fit1$aiptw$cov))
  expect_true(is.numeric(fit1$aiptw_c$est))
  expect_true(is.numeric(fit1$aiptw_c$cov))
  expect_true(class(fit1$QnMod[[1]]) == "SuperLearner")
  expect_true(class(fit1$gnMod[[1]]$A[[1]]) == "SuperLearner")
  expect_true(class(fit1$QrnMod[[1]][[1]]) == "SuperLearner")
  expect_true(class(fit1$grnMod[[1]][[1]]$fm1) == "SuperLearner")
  expect_true(class(fit1$grnMod[[1]][[1]]$fm2) == "SuperLearner")
})

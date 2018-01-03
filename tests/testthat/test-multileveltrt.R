# install.packages("~/Dropbox/R/drtmle",repos=NULL,type="source")

library(drtmle)
library(SuperLearner)
context("Testing drtmle function with > 2 treatment levels")

test_that("drtmle executes as expected with multiple treatment levels", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2)) + rbinom(n, 1, plogis(W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- drtmle(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    family = gaussian(),
    stratify = TRUE,
    glm_Q = "W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn",
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE
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

  # check ci works
  ci <- ci(fit1)
  expect_true(length(fit1$drtmle$est) == 3)
  # check contrasts work
  ci2 <- ci(fit1, contrast = c(1, -1, 0))
  expect_true(row.names(ci2$drtmle) == "E[Y(0)]-E[Y(1)]")
  ci3 <- ci(fit1, contrast = c(1, 0, -1))
  expect_true(row.names(ci3$drtmle) == "E[Y(0)]-E[Y(2)]")


  # same thing but with super learner for g
  fit1 <- drtmle(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    family = gaussian(),
    stratify = TRUE,
    glm_Q = "W1 + W2",
    SL_g = c("SL.step", "SL.glm"),
    glm_Qr = "gn",
    glm_gr = "Qn",
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE
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

  # check ci works
  ci <- ci(fit1)
  expect_true(length(fit1$drtmle$est) == 3)
  # check contrasts work
  ci2 <- ci(fit1, contrast = c(1, -1, 0))
  expect_true(row.names(ci2$drtmle) == "E[Y(0)]-E[Y(1)]")
  ci3 <- ci(fit1, contrast = c(1, 0, -1))
  expect_true(row.names(ci3$drtmle) == "E[Y(0)]-E[Y(2)]")
})


test_that("adaptive_iptw executes as expected with multiple treatment levels", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2)) + rbinom(n, 1, plogis(W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- adaptive_iptw(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    glm_g = "W1 + W2",
    returnModels = TRUE,
    glm_Qr = "gn"
  )

  expect_true(is.numeric(fit1$iptw$est))
  expect_true(is.numeric(fit1$iptw_tmle$est))
  expect_true(is.numeric(fit1$iptw_os$est))
  expect_true(is.numeric(fit1$iptw_os$cov))
  expect_true(is.numeric(fit1$iptw_tmle$cov))
  # check ci works
  ci <- ci(fit1)
  expect_true(length(fit1$iptw_tmle$est) == 3)
  # check contrasts work
  ci2 <- ci(fit1, contrast = c(1, -1, 0))
  expect_true(row.names(ci2$iptw_tmle) == "E[Y(0)]-E[Y(1)]")
  ci3 <- ci(fit1, contrast = c(1, 0, -1))
  expect_true(row.names(ci3$iptw_tmle) == "E[Y(0)]-E[Y(2)]")


  # same thing but with super learner for g
  # univariate reduction with
  # all GLMs + stratify
  fit1 <- adaptive_iptw(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    SL_g = c("SL.step", "SL.step.interaction"),
    returnModels = TRUE,
    glm_Qr = "gn"
  )

  expect_true(is.numeric(fit1$iptw$est))
  expect_true(is.numeric(fit1$iptw_tmle$est))
  expect_true(is.numeric(fit1$iptw_os$est))
  expect_true(is.numeric(fit1$iptw_os$cov))
  expect_true(is.numeric(fit1$iptw_tmle$cov))
  # check ci works
  ci <- ci(fit1)
  expect_true(length(fit1$iptw_tmle$est) == 3)
  # check contrasts work
  ci2 <- ci(fit1, contrast = c(1, -1, 0))
  expect_true(row.names(ci2$iptw_tmle) == "E[Y(0)]-E[Y(1)]")
  ci3 <- ci(fit1, contrast = c(1, 0, -1))
  expect_true(row.names(ci3$iptw_tmle) == "E[Y(0)]-E[Y(2)]")
})

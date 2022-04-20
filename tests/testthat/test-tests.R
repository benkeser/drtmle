library(drtmle)
library(SuperLearner)

context("Testing wald_test.drtmle method ")
test_that("wald_test.drtmle works as expected", {
  # simulate data
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
  # fit a drtmle
  fit1 <- drtmle(
    W = W, A = A, Y = Y, a_0 = c(1, 0),
    family = binomial(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm"
  )

  # get test for each mean for only drtmle
  tmp <- wald_test(fit1)
  # correct class
  expect_true(inherits(tmp,  "wald_test.drtmle"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # get test for each mean for drtmle and tmle
  tmp <- wald_test(fit1, est = c("drtmle", "tmle", "aiptw", "aiptw_c", "gcomp"))

  # correct class
  expect_true(inherits(tmp, "wald_test.drtmle"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
  expect_true(length(tmp) == 5)
  expect_true(all(names(tmp) == c("drtmle", "tmle", "aiptw", "aiptw_c", "gcomp")))
  expect_true(nrow(tmp$drtmle) == 2)
  expect_true(nrow(tmp$tmle) == 2)

  # get test for ATE
  tmp <- wald_test(fit1, contrast = c(1, -1))
  # correct class
  expect_true(inherits(tmp, "wald_test.drtmle"))
  # correct row name
  expect_true(row.names(tmp$drtmle) == "H0:E[Y(1)]-E[Y(0)]=0")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # throws error if crazy contrast is put in
  expect_error(wald_test(fit1, contrast = c(10214, NA)))

  # get test for risk ratio
  # by inputting own contrast function
  # this computes CI on log scale and back transforms
  myContrast <- list(
    f = function(eff) {
      log(eff)
    },
    f_inv = function(eff) {
      exp(eff)
    },
    h = function(est) {
      est[1] / est[2]
    },
    fh_grad = function(est) {
      c(1 / est[1], -1 / est[2])
    }
  )
  tmp <- wald_test(fit1, contrast = myContrast, null = 1)
  expect_true(inherits(tmp, "wald_test.drtmle"))
  expect_true(row.names(tmp$drtmle) == "H0: user contrast = 1")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
})

context("Testing ci.adaptive_iptw method ")
test_that("wald_test.adaptive_iptw works as expected", {
  # simulate data
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
  # fit a adaptive_iptw
  fit1 <- adaptive_iptw(
    W = W, A = A, Y = Y, a_0 = c(1, 0),
    SL_g = c("SL.glm", "SL.mean", "SL.step"),
    SL_Qr = "SL.glm"
  )

  # get test for each mean for only drtmle
  tmp <- wald_test(fit1)
  # correct class
  expect_true(inherits(tmp, "wald_test.adaptive_iptw"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # get test for each mean for drtmle and tmle
  tmp <- wald_test(fit1, est = c("iptw_tmle", "iptw_os"))

  # correct class
  expect_true(inherits(tmp, "wald_test.adaptive_iptw"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
  expect_true(length(tmp) == 2)
  expect_true(all(names(tmp) == c("iptw_tmle", "iptw_os")))
  expect_true(nrow(tmp$iptw_tmle) == 2)
  expect_true(nrow(tmp$iptw_os) == 2)

  # get test for ATE
  tmp <- wald_test(fit1, contrast = c(1, -1))
  # correct class
  expect_true(inherits(tmp, "wald_test.adaptive_iptw"))
  # correct row name
  expect_true(row.names(tmp$iptw_tmle) == "H0:E[Y(1)]-E[Y(0)]=0")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # throws error if crazy contrast is put in
  expect_error(wald_test(fit1, contrast = c(10214, NA)))

  # get test for risk ratio
  # by inputting own contrast function
  # this computes CI on log scale and back transforms
  myContrast <- list(
    f = function(eff) {
      log(eff)
    },
    f_inv = function(eff) {
      exp(eff)
    },
    h = function(est) {
      est[1] / est[2]
    },
    fh_grad = function(est) {
      c(1 / est[1], -1 / est[2])
    }
  )
  tmp <- wald_test(fit1, contrast = myContrast, null = 1)
  expect_true(inherits(tmp, "wald_test.adaptive_iptw"))
  expect_true(row.names(tmp$iptw_tmle) == "H0: user contrast = 1")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
})

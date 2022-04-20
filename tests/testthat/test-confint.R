library(drtmle)
library(SuperLearner)

context("Testing ci.drtmle method ")
test_that("ci.drtmle works as expected", {
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

  # get confidence intervals for each mean for only drtmle
  tmp <- ci(fit1)
  # correct class
  expect_true(inherits(tmp, "ci.drtmle"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # get confidence intervals for each mean for drtmle and tmle
  tmp <- ci(fit1, est = c("drtmle", "tmle", "aiptw", "aiptw_c", "gcomp"))

  # correct class
  expect_true(inherits(tmp, "ci.drtmle"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
  expect_true(length(tmp) == 5)
  expect_true(all(names(tmp) == c("drtmle", "tmle", "aiptw", "aiptw_c", "gcomp")))
  expect_true(nrow(tmp$drtmle) == 2)
  expect_true(nrow(tmp$tmle) == 2)

  # get confidence intervals for ATE
  tmp <- ci(fit1, contrast = c(1, -1))
  # correct class
  expect_true(inherits(tmp, "ci.drtmle"))
  # correct row name
  expect_true(row.names(tmp$drtmle) == "E[Y(1)]-E[Y(0)]")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # throws error if crazy contrast is put in
  expect_error(ci(fit1, contrast = c(10214, NA)))

  # get confidence intervals for risk ratio
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
  tmp <- ci(fit1, contrast = myContrast)
  expect_true(inherits(tmp, "ci.drtmle"))
  expect_true(row.names(tmp$drtmle) == "user contrast")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
})

##### ADD TEST FOR ci.adaptive_iptw

context("Testing ci.adaptive_iptw method ")
test_that("ci.adaptive_iptw works as expected", {
  # simulate data
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))
  # fit a drtmle
  fit1 <- adaptive_iptw(
    W = W, A = A, Y = Y, a_0 = c(1, 0),
    SL_g = c("SL.glm", "SL.mean", "SL.step"),
    SL_Qr = "SL.glm"
  )

  # get confidence intervals for each
  tmp <- ci(fit1)
  # correct class
  expect_true(inherits(tmp, "ci.adaptive_iptw"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # get confidence intervals for each mean
  tmp <- ci(fit1, est = c("iptw_tmle", "iptw_os"))

  # correct class
  expect_true(inherits(tmp, "ci.adaptive_iptw"))
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
  expect_true(length(tmp) == 2)
  expect_true(all(names(tmp) == c("iptw_tmle", "iptw_os")))
  expect_true(nrow(tmp$iptw_tmle) == 2)
  expect_true(nrow(tmp$iptw_os) == 2)

  # get confidence intervals for ATE
  tmp <- ci(fit1, contrast = c(1, -1))
  # correct class
  expect_true(inherits(tmp,  "ci.adaptive_iptw"))
  # correct row name
  expect_true(row.names(tmp$iptw_tmle) == "E[Y(1)]-E[Y(0)]")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)

  # throws error if crazy contrast is put in
  expect_error(ci(fit1, contrast = c(10214, NA)))

  # get confidence intervals for risk ratio
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
  tmp <- ci(fit1, contrast = myContrast)
  expect_true(inherits(tmp, "ci.adaptive_iptw"))
  expect_true(row.names(tmp$iptw_tmle) == "user contrast")
  # no NAs
  expect_true(sum(is.na(unlist(tmp))) == 0)
})

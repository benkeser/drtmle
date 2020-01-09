# devtools::load_all("~/Dropbox/R/drtmle")
library(drtmle)
library(SuperLearner)
context("Testing drtmle function with repeated super learners")

test_that("drtmle executes as expected with stratify = TRUE", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    glm_Q = "W1 + W2",
    glm_g = "W1 + W2",
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

  # bivariate reduction with
  # all GLMs + stratify
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    glm_Q = "W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn + gn",
    guard = c("Q", "g"),
    reduction = "bivariate"
  )

  expect_true(is.numeric(fit2$gcomp$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$cov))
  expect_true(is.numeric(fit2$drtmle$est))
  expect_true(is.numeric(fit2$drtmle$cov))
  expect_true(is.numeric(fit2$aiptw$est))
  expect_true(is.numeric(fit2$aiptw$cov))
  expect_true(is.numeric(fit2$aiptw_c$est))
  expect_true(is.numeric(fit2$aiptw_c$cov))

  # univariate reduction with
  # all SL + stratify
  fit3 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    n_SL = 2,
    reduction = "univariate"
  )
  expect_true(is.numeric(fit3$gcomp$est))
  expect_true(is.numeric(fit3$tmle$est))
  expect_true(is.numeric(fit3$tmle$est))
  expect_true(is.numeric(fit3$tmle$cov))
  expect_true(is.numeric(fit3$drtmle$est))
  expect_true(is.numeric(fit3$drtmle$cov))
  expect_true(is.numeric(fit3$aiptw$est))
  expect_true(is.numeric(fit3$aiptw$cov))
  expect_true(is.numeric(fit3$aiptw_c$est))
  expect_true(is.numeric(fit3$aiptw_c$cov))

  # bivariate reduction with
  # all SL + stratify
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    n_SL = 2,
    reduction = "bivariate"
  )
  expect_true(is.numeric(fit4$gcomp$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$cov))
  expect_true(is.numeric(fit4$drtmle$est))
  expect_true(is.numeric(fit4$drtmle$cov))
  expect_true(is.numeric(fit4$aiptw$est))
  expect_true(is.numeric(fit4$aiptw$cov))
  expect_true(is.numeric(fit4$aiptw_c$est))
  expect_true(is.numeric(fit4$aiptw_c$cov))
  # bivariate reduction with
  # single SL + stratify
  # admittedly this is somewhat of an edge case
  # should be same result as fit1
  fit5 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    guard = c("Q", "g"),
    n_SL = 2,
    reduction = "bivariate",
    returnModels = TRUE, use_future = FALSE
  )
  expect_true(is.numeric(fit5$gcomp$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$cov))
  expect_true(is.numeric(fit5$drtmle$est))
  expect_true(is.numeric(fit5$drtmle$cov))
  expect_true(is.numeric(fit5$aiptw$est))
  expect_true(is.numeric(fit5$aiptw$cov))
  expect_true(is.numeric(fit5$aiptw_c$est))
  expect_true(is.numeric(fit5$aiptw_c$cov))
  expect_true(all(abs(fit5$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit5$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  # univariate reduction with
  # single SL + stratify
  fit6 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 2,
    guard = c("Q", "g"),
    reduction = "univariate"
  )
  expect_true(is.numeric(fit6$gcomp$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$cov))
  expect_true(is.numeric(fit6$drtmle$est))
  expect_true(is.numeric(fit6$drtmle$cov))
  expect_true(is.numeric(fit6$aiptw$est))
  expect_true(is.numeric(fit6$aiptw$cov))
  expect_true(is.numeric(fit6$aiptw_c$est))
  expect_true(is.numeric(fit6$aiptw_c$cov))
  expect_true(all(abs(fit6$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit6$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))
})


# --------------------------------------------------------------------

test_that("drtmle executes as expected with stratify = FALSE", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    glm_Q = "A + W1 + W2",
    glm_g = "W1 + W2",
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

  # bivariate reduction with
  # all GLMs + stratify
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    glm_Q = "A + W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn + gn",
    guard = c("Q", "g"),
    reduction = "bivariate"
  )

  expect_true(is.numeric(fit2$gcomp$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$cov))
  expect_true(is.numeric(fit2$drtmle$est))
  expect_true(is.numeric(fit2$drtmle$cov))
  expect_true(is.numeric(fit2$aiptw$est))
  expect_true(is.numeric(fit2$aiptw$cov))
  expect_true(is.numeric(fit2$aiptw_c$est))
  expect_true(is.numeric(fit2$aiptw_c$cov))

  # bivariate reduction with
  # all SL + stratify
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    n_SL = 2,
    reduction = "bivariate"
  )
  expect_true(is.numeric(fit4$gcomp$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$cov))
  expect_true(is.numeric(fit4$drtmle$est))
  expect_true(is.numeric(fit4$drtmle$cov))
  expect_true(is.numeric(fit4$aiptw$est))
  expect_true(is.numeric(fit4$aiptw$cov))
  expect_true(is.numeric(fit4$aiptw_c$est))
  expect_true(is.numeric(fit4$aiptw_c$cov))
  # bivariate reduction with
  # single SL + stratify
  fit5 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 2,
    guard = c("Q", "g"),
    reduction = "bivariate"
  )
  expect_true(is.numeric(fit5$gcomp$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$cov))
  expect_true(is.numeric(fit5$drtmle$est))
  expect_true(is.numeric(fit5$drtmle$cov))
  expect_true(is.numeric(fit5$aiptw$est))
  expect_true(is.numeric(fit5$aiptw$cov))
  expect_true(is.numeric(fit5$aiptw_c$est))
  expect_true(is.numeric(fit5$aiptw_c$cov))
  expect_true(all(abs(fit5$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit5$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  # univariate reduction with
  # single SL + stratify
  fit6 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 2,
    guard = c("Q", "g"),
    reduction = "univariate"
  )
  expect_true(is.numeric(fit6$gcomp$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$cov))
  expect_true(is.numeric(fit6$drtmle$est))
  expect_true(is.numeric(fit6$drtmle$cov))
  expect_true(is.numeric(fit6$aiptw$est))
  expect_true(is.numeric(fit6$aiptw$cov))
  expect_true(is.numeric(fit6$aiptw_c$est))
  expect_true(is.numeric(fit6$aiptw_c$cov))
  expect_true(all(abs(fit6$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit6$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))
})


# --------------------------------------------------------------------

test_that("drtmle executes as expected with stratify = FALSE and multiple cvFolds", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  set.seed(123456)
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    glm_Q = "A + W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn",
    guard = c("Q", "g"),
    reduction = "univariate",
    cvFolds = 3
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

  # bivariate reduction with
  # all GLMs + stratify
  set.seed(123456)
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    glm_Q = "A + W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn + gn",
    guard = c("Q", "g"),
    reduction = "bivariate",
    cvFolds = 3, use_future = FALSE, returnModels = TRUE
  )

  expect_true(is.numeric(fit2$gcomp$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$est))
  expect_true(is.numeric(fit2$tmle$cov))
  expect_true(is.numeric(fit2$drtmle$est))
  expect_true(is.numeric(fit2$drtmle$cov))
  expect_true(is.numeric(fit2$aiptw$est))
  expect_true(is.numeric(fit2$aiptw$cov))
  expect_true(is.numeric(fit2$aiptw_c$est))
  expect_true(is.numeric(fit2$aiptw_c$cov))

  # bivariate reduction with
  # all SL + stratify
  set.seed(123456)
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    n_SL = 2,
    reduction = "bivariate",
    cvFolds = 3
  )
  expect_true(is.numeric(fit4$gcomp$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$est))
  expect_true(is.numeric(fit4$tmle$cov))
  expect_true(is.numeric(fit4$drtmle$est))
  expect_true(is.numeric(fit4$drtmle$cov))
  expect_true(is.numeric(fit4$aiptw$est))
  expect_true(is.numeric(fit4$aiptw$cov))
  expect_true(is.numeric(fit4$aiptw_c$est))
  expect_true(is.numeric(fit4$aiptw_c$cov))
  # bivariate reduction with
  # single SL + stratify
  set.seed(123456)
  fit5 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 2,
    guard = c("Q", "g"),
    reduction = "bivariate",
    cvFolds = 3, returnModels = TRUE, use_future = FALSE
  )
  expect_true(is.numeric(fit5$gcomp$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$est))
  expect_true(is.numeric(fit5$tmle$cov))
  expect_true(is.numeric(fit5$drtmle$est))
  expect_true(is.numeric(fit5$drtmle$cov))
  expect_true(is.numeric(fit5$aiptw$est))
  expect_true(is.numeric(fit5$aiptw$cov))
  expect_true(is.numeric(fit5$aiptw_c$est))
  expect_true(is.numeric(fit5$aiptw_c$cov))
  expect_true(all(abs(fit5$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit5$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit5$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  # univariate reduction with
  # single SL + stratify
  set.seed(123456)
  fit6 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 2,
    guard = c("Q", "g"),
    reduction = "univariate",
    cvFolds = 3
  )
  expect_true(is.numeric(fit6$gcomp$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$est))
  expect_true(is.numeric(fit6$tmle$cov))
  expect_true(is.numeric(fit6$drtmle$est))
  expect_true(is.numeric(fit6$drtmle$cov))
  expect_true(is.numeric(fit6$aiptw$est))
  expect_true(is.numeric(fit6$aiptw$cov))
  expect_true(is.numeric(fit6$aiptw_c$est))
  expect_true(is.numeric(fit6$aiptw_c$cov))
  expect_true(all(abs(fit6$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit6$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit6$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))
})

# --------------------------------------------------------------------

test_that("drtmle executes as expected with multiple sl and multiple tx lvl", {
  set.seed(123456)
  n <- 400
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 2, plogis(W$W1 - W$W2))
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

  # univariate reduction with
  # all GLMs + stratify
  fit2 <- drtmle(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    family = gaussian(),
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 3,
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE
  )
  expect_true(all(abs(fit2$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit2$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit2$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit2$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))
})

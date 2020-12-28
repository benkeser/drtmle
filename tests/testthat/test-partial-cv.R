# devtools::load_all("~/Dropbox/R/drtmle")
library(drtmle)
library(SuperLearner)
context("Testing drtmle function with partial CV")

test_that("Same point estimates, different variance estimates, single treatment", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )

  set.seed(1234)
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit1$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit1$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit1$drtmle$cov - fit2$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$tmle$cov - fit2$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw$cov - fit2$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw_c$cov - fit2$aiptw_c$cov) < 1e-7)))

  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit3 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )

  set.seed(1234)
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit3$drtmle$est - fit4$drtmle$est) < 1e-5))
  expect_true(all(abs(fit3$tmle$est - fit4$tmle$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw$est - fit4$aiptw$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw_c$est - fit4$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit3$drtmle$cov - fit4$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$tmle$cov - fit4$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw$cov - fit4$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw_c$cov - fit4$aiptw_c$cov) < 1e-7)))
})

test_that("Same point estimates, different variance estimates, multi treatment", {
  set.seed(123456)
  n <- 300
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  A[1:50] <- 2
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )

  set.seed(1234)
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit1$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit1$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit1$drtmle$cov - fit2$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$tmle$cov - fit2$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw$cov - fit2$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw_c$cov - fit2$aiptw_c$cov) < 1e-7)))

  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit3 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )

  set.seed(1234)
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit3$drtmle$est - fit4$drtmle$est) < 1e-5))
  expect_true(all(abs(fit3$tmle$est - fit4$tmle$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw$est - fit4$aiptw$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw_c$est - fit4$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit3$drtmle$cov - fit4$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$tmle$cov - fit4$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw$cov - fit4$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw_c$cov - fit4$aiptw_c$cov) < 1e-7)))
})


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
    returnModels = TRUE,
    use_future = FALSE
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
    n_SL = 2, avg_over = "drtmle", 
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )
  expect_true(all(abs(fit2$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit2$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit2$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit2$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))

  # make sure same as fit 1
  fit3 <- drtmle(
    W = W, A = A, Y = Y, a_0 = c(0, 1, 2),
    family = gaussian(),
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    n_SL = 1, avg_over = "drtmle", 
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )
  expect_true(all(abs(fit3$drtmle$est - fit1$drtmle$est) < 1e-5))
  expect_true(all(abs(fit3$tmle$est - fit1$tmle$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw$est - fit1$aiptw$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw_c$est - fit1$aiptw_c$est) < 1e-5))
})

test_that("Same point estimates, different variance estimates, missing data treatment", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)
  A[1:20] <- NA; Y[21:40] <- NA
  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
  )

  set.seed(1234)
  fit2 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = TRUE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit1$drtmle$est - fit2$drtmle$est) < 1e-5))
  expect_true(all(abs(fit1$tmle$est - fit2$tmle$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw$est - fit2$aiptw$est) < 1e-5))
  expect_true(all(abs(fit1$aiptw_c$est - fit2$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit1$drtmle$cov - fit2$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$tmle$cov - fit2$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw$cov - fit2$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit1$aiptw_c$cov - fit2$aiptw_c$cov) < 1e-7)))

  # univariate reduction with
  # all SL + stratify
  set.seed(1234)
  fit3 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    reduction = "univariate",
    returnModels = FALSE,
    use_future = FALSE
  )

  set.seed(1234)
  fit4 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.mean"),
    SL_g = c("SL.glm", "SL.mean"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
    se_cv = "partial", se_cvFolds = 10,
    targeted_se = FALSE, 
    reduction = "univariate",
    returnModels = FALSE,
    use_future = FALSE
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

  # check for same point estimates
  expect_true(all(abs(fit3$drtmle$est - fit4$drtmle$est) < 1e-5))
  expect_true(all(abs(fit3$tmle$est - fit4$tmle$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw$est - fit4$aiptw$est) < 1e-5))
  expect_true(all(abs(fit3$aiptw_c$est - fit4$aiptw_c$est) < 1e-5))

  expect_true(all(!(abs(fit3$drtmle$cov - fit4$drtmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$tmle$cov - fit4$tmle$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw$cov - fit4$aiptw$cov) < 1e-7)))
  expect_true(all(!(abs(fit3$aiptw_c$cov - fit4$aiptw_c$cov) < 1e-7)))
})

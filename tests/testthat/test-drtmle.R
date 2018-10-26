# install.packages("~/Dropbox/R/drtmle",repos=NULL,type="source")
library(drtmle)
library(SuperLearner)
context("Testing drtmle function")
test_that("drtmle executes as expected with parallel = TRUE", {
  skip_on_os("windows") # Windows doesn't support multicore (esp. Appveyor CI)
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)

  # univariate reduction with
  # all GLMs + stratify
  fit1 <- drtmle(
    W = W, A = A, Y = Y,
    parallel = TRUE,
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
})
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
    glm_gr = "Qn",
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
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
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
    fit5 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = TRUE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    guard = c("Q", "g"),
    reduction = "bivariate",
    use_future = FALSE
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
    stratify = FALSE,
    glm_Q = "W1 + W2",
    glm_g = "W1 + W2",
    glm_Qr = "gn",
    glm_gr = "Qn",
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
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
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
    stratify = FALSE,
    SL_Q = c("SL.glm", "SL.step"),
    SL_g = c("SL.glm", "SL.step"),
    SL_Qr = c("SL.glm", "SL.mean"),
    SL_gr = c("SL.glm", "SL.mean"),
    guard = c("Q", "g"),
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

  # univariate reduction with
  # single SL + stratify + Qsteps = 1
  fit7 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    guard = c("Q", "g"),
    reduction = "univariate",
    Qsteps = 1
  )
  expect_true(is.numeric(fit7$gcomp$est))
  expect_true(is.numeric(fit7$tmle$est))
  expect_true(is.numeric(fit7$tmle$est))
  expect_true(is.numeric(fit7$tmle$cov))
  expect_true(is.numeric(fit7$drtmle$est))
  expect_true(is.numeric(fit7$drtmle$cov))
  expect_true(is.numeric(fit7$aiptw$est))
  expect_true(is.numeric(fit7$aiptw$cov))
  expect_true(is.numeric(fit7$aiptw_c$est))
  expect_true(is.numeric(fit7$aiptw_c$cov))

  # bivariate reduction with
  # single SL + stratify + Qsteps = 1
  fit8 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    guard = c("Q", "g"),
    reduction = "bivariate",
    Qsteps = 1
  )
  expect_true(is.numeric(fit8$gcomp$est))
  expect_true(is.numeric(fit8$tmle$est))
  expect_true(is.numeric(fit8$tmle$est))
  expect_true(is.numeric(fit8$tmle$cov))
  expect_true(is.numeric(fit8$drtmle$est))
  expect_true(is.numeric(fit8$drtmle$cov))
  expect_true(is.numeric(fit8$aiptw$est))
  expect_true(is.numeric(fit8$aiptw$cov))
  expect_true(is.numeric(fit8$aiptw_c$est))
  expect_true(is.numeric(fit8$aiptw_c$cov))

  # give one a go with family = binomial()
  # bivariate reduction with
  # single SL + stratify + Qsteps = 1
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rbinom(n, 1, plogis(W$W1 * W$W2 * A))

  fit9 <- drtmle(
    W = W, A = A, Y = Y,
    family = binomial(),
    stratify = FALSE,
    SL_Q = "SL.glm",
    SL_g = "SL.glm",
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    guard = c("Q", "g"),
    reduction = "bivariate",
    Qsteps = 1
  )
  expect_true(is.numeric(fit9$gcomp$est))
  expect_true(is.numeric(fit9$tmle$est))
  expect_true(is.numeric(fit9$tmle$est))
  expect_true(is.numeric(fit9$tmle$cov))
  expect_true(is.numeric(fit9$drtmle$est))
  expect_true(is.numeric(fit9$drtmle$cov))
  expect_true(is.numeric(fit9$aiptw$est))
  expect_true(is.numeric(fit9$aiptw$cov))
  expect_true(is.numeric(fit9$aiptw_c$est))
  expect_true(is.numeric(fit9$aiptw_c$cov))
})


# --------------------------------------------------------------------

test_that("drtmle executes when user inputs Qn and gn and returnModels = TRUE", {
  set.seed(123456)
  n <- 100
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)
  Qn <- list(rnorm(n), rnorm(n))
  tmp <- runif(n)
  gn <- list(tmp, 1 - tmp)

    fit9 <- drtmle(
    W = W, A = A, Y = Y,
    family = gaussian(),
    stratify = FALSE,
    Qn = Qn, gn = gn,
    SL_Qr = "SL.glm",
    SL_gr = "SL.glm",
    returnModels = TRUE,
    guard = c("Q", "g"),
    reduction = "univariate",
    Qsteps = 2
  )
  expect_true(is.numeric(fit9$gcomp$est))
  expect_true(is.numeric(fit9$tmle$est))
  expect_true(is.numeric(fit9$tmle$est))
  expect_true(is.numeric(fit9$tmle$cov))
  expect_true(is.numeric(fit9$drtmle$est))
  expect_true(is.numeric(fit9$drtmle$cov))
  expect_true(is.numeric(fit9$aiptw$est))
  expect_true(is.numeric(fit9$aiptw$cov))
  expect_true(is.numeric(fit9$aiptw_c$est))
  expect_true(is.numeric(fit9$aiptw_c$cov))
})

# --------------------------------------------------------------------

test_that("GitHub error #16 resolves", {
  set.seed(123456)
  X <- runif(100, 0, 1)

Q <- X

g <- exp(X) / (1 + exp(X))

A <- rbinom(100, 1, g)

Y <- runif(Q, -0.1, 0.1)

X <- as.data.frame(X)

a <- drtmle(W = X, A = A, Y = Y, a_0 = 1, glm_Q = 'X', 
            glm_g = 'X', SL_Qr = 'SL.npreg', 
            guard = 'Q', returnModel = TRUE)
  expect_true(is.numeric(a$gcomp$est))
  expect_true(is.numeric(a$tmle$est))
  expect_true(is.numeric(a$tmle$est))
  expect_true(is.numeric(a$tmle$cov))
  expect_true(is.numeric(a$drtmle$est))
  expect_true(is.numeric(a$drtmle$cov))
  expect_true(is.numeric(a$aiptw$est))
  expect_true(is.numeric(a$aiptw$cov))
  expect_true(is.numeric(a$aiptw_c$est))
  expect_true(is.numeric(a$aiptw_c$cov))
# and the converse 
b <- drtmle(W = X, A = A, Y = Y, a_0 = 1, glm_Q = 'X', 
            glm_g = 'X', SL_gr = 'SL.npreg', 
            guard = 'g', returnModel = TRUE)
  expect_true(is.numeric(b$gcomp$est))
  expect_true(is.numeric(b$tmle$est))
  expect_true(is.numeric(b$tmle$est))
  expect_true(is.numeric(b$tmle$cov))
  expect_true(is.numeric(b$drtmle$est))
  expect_true(is.numeric(b$drtmle$cov))
  expect_true(is.numeric(b$aiptw$est))
  expect_true(is.numeric(b$aiptw$cov))
  expect_true(is.numeric(b$aiptw_c$est))
  expect_true(is.numeric(b$aiptw_c$cov))

})
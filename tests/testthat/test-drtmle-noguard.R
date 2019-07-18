# install.packages("~/Dropbox/R/drtmle",repos=NULL,type="source")
library(drtmle)
library(SuperLearner)
context("Testing drtmle function with no guarding")
test_that("drtmle/tmle and aiptw/aiptw_c give same results", {
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
    guard = NULL,
    reduction = "univariate"
  )

  expect_true(all(fit1$drtmle$est == fit1$tmle$est))
  expect_true(all(fit1$aiptw$est == fit1$aiptw_c$est))
})

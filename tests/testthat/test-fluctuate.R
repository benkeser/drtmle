library(drtmle)
library(SuperLearner)
library(np)

context("Testing fail safes in fluctuation")
test_that("Fail safe kicks in", {
  set.seed(123456)
  n <- 200
  W <- data.frame(W1 = runif(n), W2 = rnorm(n))
  A <- rbinom(n, 1, plogis(W$W1 - W$W2))
  Y <- rnorm(n, W$W1 * W$W2 * A, 2)
  n <- length(Y)
  tolg <- 1e-3
  SL_Q <- SL_g <- SL_gr <- NULL
  glm_Q <- glm_g <- "W1 + W2"
  glm_gr <- "Qn"
  cvFolds <- 1
  parallel <- FALSE
  a_0 <- 1
  DeltaA <- rep(1, n)
  DeltaY <- rep(1, n)
  returnModels <- FALSE
  stratify <- TRUE
  family <- gaussian()
  reduction <- "univariate"
  # guts of drtmle function
  if (cvFolds != 1) {
    validRows <- split(sample(1:n), rep(1:cvFolds, length = n))
  } else {
    validRows <- list(1:n)
  }
  # -------------------------------
  # estimate propensity score
  # -------------------------------
  if (!parallel) {
    gnOut <- lapply(
      X = validRows, FUN = estimateG,
      A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
      tolg = tolg, verbose = verbose, stratify = stratify,
      returnModels = returnModels, SL_g = SL_g,
      glm_g = glm_g, a_0 = a_0
    )
  } else {
    gnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
      estimateG(
        A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, verbose = verbose, stratify = stratify,
        returnModels = returnModels, SL_g = SL_g,
        glm_g = glm_g, a_0 = a_0, validRows = validRows[[v]]
      )
    }
  }

  # # re-order predictions
  gnValid <- unlist(gnOut, recursive = FALSE, use.names = FALSE)
  gnUnOrd <- do.call(Map, c(c, gnValid[seq(1, length(gnValid), 2)]))
  gn <- vector(mode = "list", length = length(a_0))
  for (i in 1:length(a_0)) {
    gn[[i]] <- rep(NA, n)
    gn[[i]][unlist(validRows)] <- gnUnOrd[[i]]
  }

  # -------------------------------
  # estimate outcome regression
  # -------------------------------
  if (!parallel) {
    QnOut <- lapply(
      X = validRows, FUN = estimateQ,
      Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
      verbose = verbose, returnModels = returnModels,
      SL_Q = SL_Q, a_0 = a_0, stratify = stratify,
      glm_Q = glm_Q, family = family
    )
  } else {
    QnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
      estimateQ(
        Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
        verbose = verbose, returnModels = returnModels,
        SL_Q = SL_Q, a_0 = a_0, glm_Q = glm_Q, family = family,
        stratify = stratify, validRows = validRows[[v]]
      )
    }
  }
  # re-order predictions
  QnValid <- unlist(QnOut, recursive = FALSE, use.names = FALSE)
  QnUnOrd <- do.call(Map, c(c, QnValid[seq(1, length(QnValid), 2)]))
  Qn <- vector(mode = "list", length = length(a_0))
  for (i in 1:length(a_0)) {
    Qn[[i]] <- rep(NA, n)
    Qn[[i]][unlist(validRows)] <- QnUnOrd[[i]]
  }

  if (!parallel) {
    grnOut <- lapply(
      X = validRows, FUN = estimategrn,
      Y = Y, A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY,
      tolg = tolg, Qn = Qn, gn = gn,
      glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
      reduction = reduction, returnModels = returnModels
    )
  } else {
    grnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
      estimategrn(
        Y = Y, A = A, W = W,
        DeltaA = DeltaA, DeltaY = DeltaY,
        tolg = tolg, Qn = Qn, gn = gn,
        glm_gr = glm_gr, SL_gr = SL_gr, a_0 = a_0,
        reduction = reduction, returnModels = returnModels,
        validRows = validRows[[v]]
      )
    }
  }
  # re-order predictions
  grnValid <- unlist(grnOut, recursive = FALSE, use.names = FALSE)
  grnUnOrd <- do.call(Map, c(rbind, grnValid[seq(1, length(grnValid), 2)]))
  grn <- vector(mode = "list", length = length(a_0))
  for (i in 1:length(a_0)) {
    grn[[i]] <- data.frame(grn1 = rep(NA, n), grn2 = rep(NA, n))
    grn[[i]][unlist(validRows), ] <- cbind(grnUnOrd[[i]])
  }



  # ----------------------------------------------
  # now call fluctuateQ's setting coefTol = 0,
  # which should trigger fail-safes
  # ----------------------------------------------
  grbg <- fluctuateQ2(
    Y, A, W, DeltaY, DeltaA,
    Qn, gn, grn, a_0, reduction,
    coefTol = 0
  )
  # set tolerance threshold
  tol <- 1e-4
  # should just return Qn
  expect_true(all(grbg[[1]]$est - Qn[[1]] < tol))

  grbg2 <- fluctuateQ(
    Y, A, W, DeltaY, DeltaA,
    Qn, gn, grn, a_0, reduction,
    coefTol = 0
  )
  expect_true(all(grbg2[[1]]$est - Qn[[1]] < tol))
})

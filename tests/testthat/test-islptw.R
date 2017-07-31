library(drtmle)
library(SuperLearner)
library(np)

# TO DO: Add tests for multiple treatment levels
context("Testing islptw works")
test_that("islptw works as expected with parallel = TRUE",{
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	fit1 <- islptw(W = W, A = A, Y = Y, 
	               a_0 = c(0,1),
                  glm_g="W1 + W2",
                  glm_Qr="gn",parallel = TRUE)
	expect_true(all(!is.na(fit1$islptw_tmle$est)))
	expect_true(all(!is.na(fit1$islptw_tmle$cov)))
	expect_true(all(!is.na(fit1$islptw_os$est)))
	expect_true(all(!is.na(fit1$islptw_os$cov)))
	expect_true(all(!is.na(fit1$iptw$est)))
	expect_true(all(!is.na(fit1$iptw$cov)))
	expect_true(class(fit1)=="islptw")
})
test_that("islptw works as expected",{
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	fit1 <- islptw(W = W, A = A, Y = Y, 
	               a_0 = c(0,1),
                  glm_g="W1 + W2",
                  glm_Qr="gn")
	expect_true(all(!is.na(fit1$islptw_tmle$est)))
	expect_true(all(!is.na(fit1$islptw_tmle$cov)))
	expect_true(all(!is.na(fit1$islptw_os$est)))
	expect_true(all(!is.na(fit1$islptw_os$cov)))
	expect_true(all(!is.na(fit1$iptw$est)))
	expect_true(all(!is.na(fit1$iptw$cov)))
	expect_true(class(fit1)=="islptw")
})	

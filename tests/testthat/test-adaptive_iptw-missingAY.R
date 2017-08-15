library(drtmle)
library(SuperLearner)
library(np)
library(testthat)

# TO DO: Add tests for multiple treatment levels
context("Testing adaptive_iptw works")
test_that("adaptive_iptw works as expected with missing data",{
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	DeltaA <- 1-rbinom(n,1,plogis(-4 + W$W1 - W$W2))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	DeltaY <- 1-rbinom(n,1,plogis(-4 + W$W1 - A))
	Y <- rnorm(n, W$W1*W$W2*A, 2)
	Y[DeltaY == 0] <- NA
	A[DeltaA == 0] <- NA

	fit1 <- adaptive_iptw(W = W, A = A, Y = Y, 
	               a_0 = c(0,1),
                  glm_g="W1 + W2",
                  glm_Qr="gn")
	expect_true(all(!is.na(fit1$iptw_tmle$est)))
	expect_true(all(!is.na(fit1$iptw_tmle$cov)))
	expect_true(all(!is.na(fit1$iptw_os$est)))
	expect_true(all(!is.na(fit1$iptw_os$cov)))
	expect_true(all(!is.na(fit1$iptw$est)))
	expect_true(class(fit1)=="adaptive_iptw")
})

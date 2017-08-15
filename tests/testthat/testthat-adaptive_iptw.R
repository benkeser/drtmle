library(drtmle)
library(SuperLearner)
library(np)
library(testthat)
# TO DO: Add tests for multiple treatment levels
context("Testing adaptive_iptw works with cv")

test_that("adaptive_iptw works as expected with cv",{
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	fit1 <- adaptive_iptw(W = W, A = A, Y = Y,
					cvFolds = 2,  
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

library(drtmle)
library(SuperLearner)

context("Testing confint.drtmle method ")
test_that("confint.drtmle works as expected",{
	# simulate data
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
	# fit a drtmle
	fit1 <- drtmle(W = W, A = A, Y = Y, a_0 = c(1,0),
					  family=binomial(),
	               stratify=FALSE,
	               SL_Q="SL.glm",
	               SL_g="SL.glm",
	               SL_Qr="SL.glm",
	               SL_gr="SL.glm")

	# get confidence intervals for each mean for only drtmle
	tmp <- confint(fit1)
	# correct class
	expect_true(class(tmp)=="confint.drtmle")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)

	# get confidence intervals for each mean for drtmle and tmle
	tmp <- confint(fit1, parm = c("drtmle","tmle","aiptw","aiptw_c","gcomp"))

	# correct class
	expect_true(class(tmp)=="confint.drtmle")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)
	expect_true(length(tmp) == 5)
	expect_true(all(names(tmp) == c("drtmle","tmle","aiptw","aiptw_c","gcomp")))
	expect_true(nrow(tmp$drtmle) == 2)
	expect_true(nrow(tmp$tmle) == 2)

	# get confidence intervals for ATE
	tmp <- confint(fit1, contrast = c(1,-1))
	# correct class
	expect_true(class(tmp)=="confint.drtmle")
	# correct row name
	expect_true(row.names(tmp$drtmle) == "E[Y(1)]-E[Y(0)]")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)

	# throws error if crazy contrast is put in 
	expect_error(confint(fit1, contrast = c(10214,NA)))

	# get confidence intervals for risk ratio 
	# by inputting own contrast function
	# this computes CI on log scale and back transforms
	myContrast <- list(f = function(eff){ log(eff) },
	                   f_inv = function(eff){ exp(eff) },
	                   h = function(parm){ parm[1]/parm[2] },
	                   h_grad =  function(parm){ c(1/parm[1],-1/parm[2]) })
	tmp <- confint(fit1, contrast = myContrast)
	expect_true(class(tmp)=="confint.drtmle")
	expect_true(row.names(tmp$drtmle) == "user contrast")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)
})

##### ADD TEST FOR confint.islptw

context("Testing confint.drtmle method ")
test_that("confint.islptw works as expected",{
	# simulate data
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
	# fit a drtmle
	fit1 <- islptw(W = W, A = A, Y = Y, a_0 = c(1,0),
	               SL_g=c("SL.glm","SL.mean","SL.step"),
	               SL_Qr="SL.glm")

	# get confidence intervals for each 
	tmp <- confint(fit1)
	# correct class
	expect_true(class(tmp)=="confint.islptw")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)

	# get confidence intervals for each mean 
	tmp <- confint(fit1, parm = c("islptw_tmle","islptw_os"))

	# correct class
	expect_true(class(tmp)=="confint.islptw")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)
	expect_true(length(tmp) == 2)
	expect_true(all(names(tmp) == c("islptw_tmle","islptw_os")))
	expect_true(nrow(tmp$islptw_tmle) == 2)
	expect_true(nrow(tmp$islptw_os) == 2)

	# get confidence intervals for ATE
	tmp <- confint(fit1, contrast = c(1,-1))
	# correct class
	expect_true(class(tmp)=="confint.islptw")
	# correct row name
	expect_true(row.names(tmp$islptw_tmle) == "E[Y(1)]-E[Y(0)]")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)

	# throws error if crazy contrast is put in 
	expect_error(confint(fit1, contrast = c(10214,NA)))

	# get confidence intervals for risk ratio 
	# by inputting own contrast function
	# this computes CI on log scale and back transforms
	myContrast <- list(f = function(eff){ log(eff) },
	                   f_inv = function(eff){ exp(eff) },
	                   h = function(parm){ parm[1]/parm[2] },
	                   h_grad =  function(parm){ c(1/parm[1],-1/parm[2]) })
	tmp <- confint(fit1, contrast = myContrast)
	expect_true(class(tmp)=="confint.islptw")
	expect_true(row.names(tmp$islptw_tmle) == "user contrast")
	# no NAs
	expect_true(sum(is.na(unlist(tmp)))==0)
})

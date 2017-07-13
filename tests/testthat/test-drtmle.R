# install.packages("~/Dropbox/R/drtmle",repos=NULL,type="source")
library(drtmle)

context("Testing drtmle function")

test_that("drtmle executes as expected with stratify = TRUE", {
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	# univariate reduction with
	# all GLMs + stratify
	fit1 <- drtmle(W = W, A = A, Y = Y, 
                   family=gaussian(),
                      stratify=TRUE,
                      glmQ="W1 + W2",
                      glmg="W1 + W2",
                      glmQr="gn",
                      glmgr="Qn",
                      guard=c("Q","g"),
                      reduction="univariate")

	expect_true(is.numeric(fit1$naive$est))
	expect_true(is.numeric(fit1$tmle$est))
	expect_true(is.numeric(fit1$tmle$est))
	expect_true(is.numeric(fit1$tmle$cov))
	expect_true(is.numeric(fit1$tmle.dral$est))
	expect_true(is.numeric(fit1$tmle.dral$cov))
	expect_true(is.numeric(fit1$os$est))
	expect_true(is.numeric(fit1$os$cov))
	expect_true(is.numeric(fit1$os.dral$est))
	expect_true(is.numeric(fit1$os.dral$cov))

	# bivariate reduction with 
	# all GLMs + stratify 
	fit2 <- drtmle(W = W, A = A, Y = Y, 
                   family=gaussian(),
                      stratify=TRUE,
                      glmQ="W1 + W2",
                      glmg="W1 + W2",
                      glmQr="gn",
                      glmgr="Qn",
                      guard=c("Q","g"),
                      reduction="bivariate")

	expect_true(is.numeric(fit2$naive$est))
	expect_true(is.numeric(fit2$tmle$est))
	expect_true(is.numeric(fit2$tmle$est))
	expect_true(is.numeric(fit2$tmle$cov))
	expect_true(is.numeric(fit2$tmle.dral$est))
	expect_true(is.numeric(fit2$tmle.dral$cov))
	expect_true(is.numeric(fit2$os$est))
	expect_true(is.numeric(fit2$os$cov))
	expect_true(is.numeric(fit2$os.dral$est))
	expect_true(is.numeric(fit2$os.dral$cov))

	# univariate reduction with 
	# all SL + stratify
	fit3 <- drtmle(W = W, A = A, Y = Y, 
               family=gaussian(),
                  stratify=TRUE,
                  libraryQ=c("SL.glm","SL.step"),
                  libraryg=c("SL.glm","SL.step"),
                  libraryQr=c("SL.glm","SL.mean"),
                  librarygr=c("SL.glm","SL.mean"),
                  guard=c("Q","g"),
                  reduction="univariate")
	expect_true(is.numeric(fit3$naive$est))
	expect_true(is.numeric(fit3$tmle$est))
	expect_true(is.numeric(fit3$tmle$est))
	expect_true(is.numeric(fit3$tmle$cov))
	expect_true(is.numeric(fit3$tmle.dral$est))
	expect_true(is.numeric(fit3$tmle.dral$cov))
	expect_true(is.numeric(fit3$os$est))
	expect_true(is.numeric(fit3$os$cov))
	expect_true(is.numeric(fit3$os.dral$est))
	expect_true(is.numeric(fit3$os.dral$cov))

	#bivariate reduction with 
	# all SL + stratify
	fit4 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=TRUE,
              libraryQ=c("SL.glm","SL.step"),
              libraryg=c("SL.glm","SL.step"),
              libraryQr=c("SL.glm","SL.mean"),
              librarygr=c("SL.glm","SL.mean"),
              guard=c("Q","g"),
              reduction="bivariate")
	expect_true(is.numeric(fit4$naive$est))
	expect_true(is.numeric(fit4$tmle$est))
	expect_true(is.numeric(fit4$tmle$est))
	expect_true(is.numeric(fit4$tmle$cov))
	expect_true(is.numeric(fit4$tmle.dral$est))
	expect_true(is.numeric(fit4$tmle.dral$cov))
	expect_true(is.numeric(fit4$os$est))
	expect_true(is.numeric(fit4$os$cov))
	expect_true(is.numeric(fit4$os.dral$est))
	expect_true(is.numeric(fit4$os.dral$cov))
	# bivariate reduction with 
	# single SL + stratify
	fit5 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=TRUE,
              libraryQ="SL.glm",
              libraryg="SL.glm",
              libraryQr="SL.glm",
              librarygr="SL.glm",
              guard=c("Q","g"),
              reduction="bivariate")
	expect_true(is.numeric(fit5$naive$est))
	expect_true(is.numeric(fit5$tmle$est))
	expect_true(is.numeric(fit5$tmle$est))
	expect_true(is.numeric(fit5$tmle$cov))
	expect_true(is.numeric(fit5$tmle.dral$est))
	expect_true(is.numeric(fit5$tmle.dral$cov))
	expect_true(is.numeric(fit5$os$est))
	expect_true(is.numeric(fit5$os$cov))
	expect_true(is.numeric(fit5$os.dral$est))
	expect_true(is.numeric(fit5$os.dral$cov))
	# univariate reduction with 
	# single SL + stratify
	fit6 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=TRUE,
              libraryQ="SL.glm",
              libraryg="SL.glm",
              libraryQr="SL.glm",
              librarygr="SL.glm",
              guard=c("Q","g"),
              reduction="univariate")
	expect_true(is.numeric(fit6$naive$est))
	expect_true(is.numeric(fit6$tmle$est))
	expect_true(is.numeric(fit6$tmle$est))
	expect_true(is.numeric(fit6$tmle$cov))
	expect_true(is.numeric(fit6$tmle.dral$est))
	expect_true(is.numeric(fit6$tmle.dral$cov))
	expect_true(is.numeric(fit6$os$est))
	expect_true(is.numeric(fit6$os$cov))
	expect_true(is.numeric(fit6$os.dral$est))
	expect_true(is.numeric(fit6$os.dral$cov))

})


#--------------------------------------------------------------------

test_that("drtmle executes as expected with stratify = FALSE", {
	set.seed(123456)
	n <- 200
	W <- data.frame(W1 = runif(n), W2 = rnorm(n))
	A <- rbinom(n,1,plogis(W$W1 - W$W2))
	Y <- rnorm(n, W$W1*W$W2*A, 2)

	# univariate reduction with
	# all GLMs + stratify
	fit1 <- drtmle(W = W, A = A, Y = Y, 
                   family=gaussian(),
                      stratify=FALSE,
                      glmQ="W1 + W2",
                      glmg="W1 + W2",
                      glmQr="gn",
                      glmgr="Qn",
                      guard=c("Q","g"),
                      reduction="univariate")

	expect_true(is.numeric(fit1$naive$est))
	expect_true(is.numeric(fit1$tmle$est))
	expect_true(is.numeric(fit1$tmle$est))
	expect_true(is.numeric(fit1$tmle$cov))
	expect_true(is.numeric(fit1$tmle.dral$est))
	expect_true(is.numeric(fit1$tmle.dral$cov))
	expect_true(is.numeric(fit1$os$est))
	expect_true(is.numeric(fit1$os$cov))
	expect_true(is.numeric(fit1$os.dral$est))
	expect_true(is.numeric(fit1$os.dral$cov))

	# bivariate reduction with 
	# all GLMs + stratify 
	fit2 <- drtmle(W = W, A = A, Y = Y, 
                   family=gaussian(),
                      stratify=FALSE,
                      glmQ="W1 + W2",
                      glmg="W1 + W2",
                      glmQr="gn",
                      glmgr="Qn",
                      guard=c("Q","g"),
                      reduction="bivariate")

	expect_true(is.numeric(fit2$naive$est))
	expect_true(is.numeric(fit2$tmle$est))
	expect_true(is.numeric(fit2$tmle$est))
	expect_true(is.numeric(fit2$tmle$cov))
	expect_true(is.numeric(fit2$tmle.dral$est))
	expect_true(is.numeric(fit2$tmle.dral$cov))
	expect_true(is.numeric(fit2$os$est))
	expect_true(is.numeric(fit2$os$cov))
	expect_true(is.numeric(fit2$os.dral$est))
	expect_true(is.numeric(fit2$os.dral$cov))

	# univariate reduction with 
	# all SL + stratify
	fit3 <- drtmle(W = W, A = A, Y = Y, 
               family=gaussian(),
                  stratify=FALSE,
                  libraryQ=c("SL.glm","SL.step"),
                  libraryg=c("SL.glm","SL.step"),
                  libraryQr=c("SL.glm","SL.mean"),
                  librarygr=c("SL.glm","SL.mean"),
                  guard=c("Q","g"),
                  reduction="univariate")
	expect_true(is.numeric(fit3$naive$est))
	expect_true(is.numeric(fit3$tmle$est))
	expect_true(is.numeric(fit3$tmle$est))
	expect_true(is.numeric(fit3$tmle$cov))
	expect_true(is.numeric(fit3$tmle.dral$est))
	expect_true(is.numeric(fit3$tmle.dral$cov))
	expect_true(is.numeric(fit3$os$est))
	expect_true(is.numeric(fit3$os$cov))
	expect_true(is.numeric(fit3$os.dral$est))
	expect_true(is.numeric(fit3$os.dral$cov))

	#bivariate reduction with 
	# all SL + stratify
	fit4 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=FALSE,
              libraryQ=c("SL.glm","SL.step"),
              libraryg=c("SL.glm","SL.step"),
              libraryQr=c("SL.glm","SL.mean"),
              librarygr=c("SL.glm","SL.mean"),
              guard=c("Q","g"),
              reduction="bivariate")
	expect_true(is.numeric(fit4$naive$est))
	expect_true(is.numeric(fit4$tmle$est))
	expect_true(is.numeric(fit4$tmle$est))
	expect_true(is.numeric(fit4$tmle$cov))
	expect_true(is.numeric(fit4$tmle.dral$est))
	expect_true(is.numeric(fit4$tmle.dral$cov))
	expect_true(is.numeric(fit4$os$est))
	expect_true(is.numeric(fit4$os$cov))
	expect_true(is.numeric(fit4$os.dral$est))
	expect_true(is.numeric(fit4$os.dral$cov))
	# bivariate reduction with 
	# single SL + stratify
	fit5 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=FALSE,
              libraryQ="SL.glm",
              libraryg="SL.glm",
              libraryQr="SL.glm",
              librarygr="SL.glm",
              guard=c("Q","g"),
              reduction="bivariate")
	expect_true(is.numeric(fit5$naive$est))
	expect_true(is.numeric(fit5$tmle$est))
	expect_true(is.numeric(fit5$tmle$est))
	expect_true(is.numeric(fit5$tmle$cov))
	expect_true(is.numeric(fit5$tmle.dral$est))
	expect_true(is.numeric(fit5$tmle.dral$cov))
	expect_true(is.numeric(fit5$os$est))
	expect_true(is.numeric(fit5$os$cov))
	expect_true(is.numeric(fit5$os.dral$est))
	expect_true(is.numeric(fit5$os.dral$cov))
	# univariate reduction with 
	# single SL + stratify
	fit6 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=FALSE,
              libraryQ="SL.glm",
              libraryg="SL.glm",
              libraryQr="SL.glm",
              librarygr="SL.glm",
              guard=c("Q","g"),
              reduction="univariate")
	expect_true(is.numeric(fit6$naive$est))
	expect_true(is.numeric(fit6$tmle$est))
	expect_true(is.numeric(fit6$tmle$est))
	expect_true(is.numeric(fit6$tmle$cov))
	expect_true(is.numeric(fit6$tmle.dral$est))
	expect_true(is.numeric(fit6$tmle.dral$cov))
	expect_true(is.numeric(fit6$os$est))
	expect_true(is.numeric(fit6$os$cov))
	expect_true(is.numeric(fit6$os.dral$est))
	expect_true(is.numeric(fit6$os.dral$cov))

	# univariate reduction with 
	# single SL + stratify + Qsteps = 1
	fit7 <- drtmle(W = W, A = A, Y = Y, 
           family=gaussian(),
              stratify=FALSE,
              libraryQ="SL.glm",
              libraryg="SL.glm",
              libraryQr="SL.glm",
              librarygr="SL.glm",
              guard=c("Q","g"),
              reduction="univariate",
              Qsteps = 1)
	expect_true(is.numeric(fit7$naive$est))
	expect_true(is.numeric(fit7$tmle$est))
	expect_true(is.numeric(fit7$tmle$est))
	expect_true(is.numeric(fit7$tmle$cov))
	expect_true(is.numeric(fit7$tmle.dral$est))
	expect_true(is.numeric(fit7$tmle.dral$cov))
	expect_true(is.numeric(fit7$os$est))
	expect_true(is.numeric(fit7$os$cov))
	expect_true(is.numeric(fit7$os.dral$est))
	expect_true(is.numeric(fit7$os.dral$cov))

})
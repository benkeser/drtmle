#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

# load libraries #
library(hal)
library(drtmle)
library(np)
library(methods)
library(caret)
library(gbm)

# parameters for simulation
n <- c(5000)
est <- c("adaptive1","adaptive2")
seed <- 1:500
parm <- expand.grid(seed=seed, n=n, est=est)

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  #   parm.red <- expand.grid(seed=seed, n=n)
  #   for(i in 1:nrow(parm.red)){
  #       set.seed(parm.red$seed[i])
  #       dat <- drtmle:::makeDataBiv(n=parm.red$n[i])
  #       save(dat, file=paste0("~/dral/scratch/dat2_n=",parm.red$n[i],
  #                             "_seed=",parm.red$seed[i],".RData"))
  #   }
  # print(paste0('initial datasets saved to: ~/dral/scratch/inFile ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    
      # load data
      load(paste0("~/dral/scratch/dat2_n=",parm$n[i],
                "_seed=",parm$seed[i],".RData"))
      
      # compute truth
      truth <- drtmle:::getTruth()
      # fit estimator 
      if(parm$est[i]=="parametric"){
          out <- vector(mode = "list")
          # cao
          out1 <- cao.dr(R=dat$A,Y=dat$Y,cov=dat$W,nBoot=100)
          names(out1$ci) <- NULL
          out$cao <- list(
              est = out1$est, ci = out1$ci,
              cov = out1$ci[1] < truth & out1$ci[2] > truth,
              err = out1$est - truth
          )
          
          # original vermeulen
          out2 <- m.biasreducedDR.identity(R=dat$A,Y=dat$Y,cov=data.matrix(dat$W))
          out2$ci <- c(out2$mn.Y - 1.96*out2$se.mn.Y, out2$mn.Y + 1.96*out2$se.mn.Y)
          out$verm1 <- list(
              est = out2$mn.Y, ci = out2$ci,
              cov = out2$ci[1] < truth & out2$ci[2] > truth,
              err = out2$mn.Y - truth
          )
          
          # tan
          library(iWeigReg)
          Qmod <- glm(dat$Y[dat$A==1] ~ ., data = dat$W[dat$A==1,])
          Qn <- predict(Qmod, type="response", newdata=dat$W)
          gmod <- glm(dat$A ~ ., data=dat$W, family="binomial")
          gn <- predict(gmod, type="response")
          
          out3 <- mn.clik(y = dat$Y, tr = dat$A, p = gn, g = cbind(1, Qn), X = model.matrix(gmod))
          out3$ci <- c(out3$mu - 1.96*sqrt(out3$v), out3$mu + 1.96*sqrt(out3$v))
          out$tan = list(
              est = out3$mu, ci = out3$ci, 
              cov = out3$ci[1] < truth & out3$ci[2] > truth,
              err = out3$mu - truth
          )
      }else if(parm$est[i] == "vv2"){
          out <- vector(mode = "list")
          # data-adaptive vermeulen
          out1 <- data.adaptive.biasreduced.DR(
              R=dat$A,Y=dat$Y,cov=data.matrix(dat$W),type.initQ=c("npreg"),
              zeta=0.005,fluc=c("unweighted"),
              alpha=0.05,psi.tilde=0)
          out$vv2 <- list(
              est = out1$est, ci = out1$ci, 
              cov = out1$ci[1] < truth & out1$ci[2] > truth,
              err = out1$est - truth
          )
      }else{
          out <- vector(mode = "list")
          
          # define GBM function
          out1 <- drtmle(Y=dat$Y, A=dat$A, W=dat$W, family=gaussian(),
                        a0=1,
                        maxIter = 3,
                        libraryQ=c("SL.npreg"),
                        libraryg=c("SL.npreg"),
                        librarygr=c("SL.npreg"),
                        libraryQr=c("SL.npreg"),
                        reduction=ifelse(parm$est[i]=="adaptive1","univariate","bivariate"),
                        tolIC=0.001, 
                        tolg=0.025
          )


          # add confidence intervals to output
          out1$os$ci <- c(
              out1$os$est - 1.96*sqrt(out1$os$cov),
              out1$os$est + 1.96*sqrt(out1$os$cov)
          )
          out1$os.dral$ci <- c(
              out1$os.dral$est - 1.96*sqrt(out1$os.dral$cov),
              out1$os.dral$est + 1.96*sqrt(out1$os.dral$cov)
          )
          out1$tmle$ci <- c(
              out1$tmle$est - 1.96*sqrt(out1$tmle$cov),
              out1$tmle$est + 1.96*sqrt(out1$tmle$cov)
          )
          
          out1$tmle.dral$ci <- c(
              out1$tmle.dral$est - 1.96*sqrt(out1$tmle.dral$cov),
              out1$tmle.dral$est + 1.96*sqrt(out1$tmle.dral$cov)
          )
          
          # only include for adaptive 1
          if(parm$est[i] == "adaptive1"){
          out$os <- list(
              est = out1$os$est, ci = out1$os$ci,
              cov = out1$os$ci[1] < truth & out1$os$ci[2] > truth,
              err = out1$os$est - truth
          )
          }
          out$os.dral <- list(
              est = out1$os.dral$est, ci = out1$os.dral$ci,
              cov = out1$os.dral$ci[1] < truth & out1$os.dral$ci[2] > truth,
              err = out1$os.dral$est - truth
          )
          if(parm$est[i] == "adaptive1"){
          out$tmle <- list(
              est = out1$tmle$est, ci = out1$tmle$ci,
              cov = out1$tmle$ci[1] < truth & out1$tmle$ci[2] > truth,
              err = out1$tmle$est - truth
          )
          }
          out$tmle.dral <- list(
              est = out1$tmle.dral$est, ci = out1$tmle.dral$ci,
              cov = out1$tmle.dral$ci[1] < truth & out1$tmle.dral$ci[2] > truth,
              err = out1$tmle.dral$est - truth
          )
          
      }
    save(out,file=paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData.tmp"))
    file.rename(paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData.tmp"),
                paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData"))
    print("file saved!")
  }
}

# merge job ###########################
if (args[1] == 'merge') {
    n <- c(250,1000,5000)
    seed <- 1:500
    parm <- expand.grid(seed=seed, n=n)
    truth <- drtmle:::getTruth()

    allOut <- NULL
    for(i in 1:nrow(parm)){
        allOut <- tryCatch({
        out.parm <- get(load(paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=parametric.RData")))
        out.vv2 <- get(load(paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=vv2.RData")))
        
        if(parm$n[i]==5000){
          out.adapt1 <- get(load(paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=adaptive1.RData")))
          out.adapt2 <- get(load(paste0("~/dral/out/newSim2_n=",parm$n[i],"_seed=",parm$seed[i],"_est=adaptive2.RData")))
        }else{
          #---------------------------------------------------
          # need to add old adaptive results to this output!
          #----------------------------------------------------
          out2 <- get(load(paste0("~/dral/out/outN_Biv_20160302_n=",parm$n[i],"_seed=",parm$seed[i],".RData")))
          out1 <- get(load(paste0("~/dral/out/outN_Uni_20160302_n=",parm$n[i],"_seed=",parm$seed[i],".RData")))
          out.adapt1 <- out.adapt2 <- vector(mode="list",length=0)
          out.adapt1$os <- list(
               est = out1$os$est, ci = out1$os$est +c(-1.96,1.96)*sqrt(out1$os$cov) ,
               cov = out1$os$est - 1.96*sqrt(out1$os$cov) < truth & out1$os$est + 1.96*sqrt(out1$os$cov) > truth,
               err = out1$os$est - truth
          )
          out.adapt1$tmle <- list(
              est = out1$tmle$est, ci =  out1$tmle$est + c(-1.96,1.96)*sqrt(out1$tmle$cov),
              cov = out1$tmle$est - 1.96*sqrt(out1$tmle$cov) < truth & out1$tmle$est + 1.96*sqrt(out1$tmle$cov) > truth,
              err = out1$tmle$est - truth
          )
          out.adapt1$tmle.dral <- list(
              est = out1$tmle.dral$est, ci = out1$tmle.dral$est + c(-1.96,1.96)*sqrt(out1$tmle.dral$cov),
              cov = out1$tmle.dral$est - 1.96*sqrt(out1$tmle.dral$cov) < truth & out1$tmle.dral$est + 1.96*sqrt(out1$tmle.dral$cov) > truth,
              err = out1$tmle.dral$est - truth
          )
          out.adapt1$os.dral <- list(
              est = out1$os.dral$est, ci = out1$os.dral$est +c(-1.96,1.96)*sqrt(out1$os.dral$cov),
              cov = out1$os.dral$est - 1.96*sqrt(out1$os.dral$cov) < truth & out1$os.dral$est + 1.96*sqrt(out1$os.dral$cov) > truth,
              err = out1$os.dral$est - truth
          )
          out.adapt2$tmle.dral <- list(
              est = out2$tmle.dral$est, ci = out2$tmle.dral$est + c(-1.96,1.96)*sqrt(out2$tmle.dral$cov),
              cov = out2$tmle.dral$est - 1.96*sqrt(out2$tmle.dral$cov) < truth & out2$tmle.dral$est + 1.96*sqrt(out2$tmle.dral$cov) > truth,
              err = out2$tmle.dral$est - truth
          )
          out.adapt2$os.dral <- list(
              est = out2$os.dral$est, ci = out2$os.dral$est + c(-1.96,1.96)*sqrt(out2$os.dral$cov) ,
              cov = out2$os.dral$est - 1.96*sqrt(out2$os.dral$cov) < truth & out2$os.dral$est + 1.96*sqrt(out2$os.dral$cov) > truth,
              err = out2$os.dral$est - truth
          )
        }

        rbind(allOut,
          #c(unlist(out.parm),unlist(out.vv2)))
          c(unlist(out.parm),unlist(out.vv2),unlist(out.adapt1),unlist(out.adapt2)))
          #c(unlist(out.adapt1),unlist(out.adapt2)))
        },error=function(e){
            rbind(allOut, rep(NA, ncol(allOut)))
        })
    }
    allOut <- cbind(parm,allOut)
    # get bias and coverage by looking at column means
    by(allOut, factor(allOut$n), colMeans, na.rm=TRUE)
        
    # get CI width
    by(allOut, factor(allOut$n), function(x){
        up <- grep("ci2", colnames(x))
        down <- grep("ci1", colnames(x))
        colMeans(x[,up] - x[,down],na.rm=TRUE)        
    })

    # get summary
    by(allOut, factor(allOut$n), function(x){
      summary(x[,grep("est",colnames(x))])
    })

    # get variance
    by(allOut, factor(allOut$n), function(x){
      apply(x[,grep("est",colnames(x))],2,var,na.rm=TRUE)*10
    })


    save(allOut, file = "~/dral/out/newSim2_allOut.RData")
}
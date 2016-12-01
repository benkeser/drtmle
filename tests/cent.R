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

# parameters for simulation
n <- c(500,2000)
est <- c("parametric","vv2","adaptive")
seed <- 1:500
parm <- expand.grid(seed=seed, n=n, est=est)

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
    for(i in 1:nrow(parm)){
        set.seed(parm$seed[i])
        dat <- makeData(n=parm$n[i])
        save(dat, file=paste0("~/dral/scratch/dat_n=",parm$n[i],
                              "_seed=",parm$seed[i],".RData"))
    }
  print(paste0('initial datasets saved to: ~/dral/scratch/inFile ... .RData'))
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
      load(paste0("~/dral/scratch/dat_n=",parm$n[i],
                "_seed=",parm$seed[i],".RData"))
      
      # fit estimator 
      if(parm$est[i]=="parametric"){
          out <- vector(mode = "list")
          # cao
          out1 <- cao.dr(R=dat$A,Y=dat$Y,cov=dat$W,nBoot=100)
          
          out$cao <- list(
              est = out1$est, ci = out1$ci,
              cov = out1$ci[1] < 210 & out1$ci[2] > 210,
              err = out1$est - 210
          )
          
          # original vermeulen
          out2 <- m.biasreducedDR.identity(R=dat$A,Y=dat$Y,cov=data.matrix(dat$W))
          out2$ci <- c(out2$mn.Y - 1.96*out2$se.mn.Y, out2$mn.Y + 1.96*out2$se.mn.Y)
          out$verm1 <- list(
              est = out2$mn.Y, ci = out2$ci,
              cov = out2$ci[1] < 210 & out2$ci[2] > 210,
              err = out2$mn.Y - 210
          )
          
          # tan
          library(iWeigReg)
          Qmod <- glm(dat$Y[dat$A==1] ~ ., data = dat$W[dat$A==1,])
          Qn <- predict(Qmod, type="response", newdata=dat$W)
          gmod <- glm(dat$A ~ ., data=dat$W, family="binomial")
          gn <- predict(gmod, type="response")
          
          out3 <- mn.clik(y = dat$Y, tr = dat$A, p = gn, g = cbind(1, Qn), X = model.matrix(gmod))
          out3$ci <- c(out3$est - 1.96*sqrt(out3$v), out3$mu + 1.96*sqrt(out3$v))
          out$tan = list(
              est = out5$mu, ci = out3$ci, 
              cov = out3$ci[1] < 210 & out3$ci[2] > 210,
              err = out3$mu - 210
          )
      }else if(parm$est[i] == "vv2"){
          out <- vector(mode = "list")
          # data-adaptive vermeulen
          out1 <- data.adaptive.biasreduced.DR(
              R=dat$A,Y=dat$Y,cov=data.matrix(dat$W),type.initQ=c("SL"),
              zeta=0.005,fluc=c("unweighted"),
              alpha=0.05,psi.tilde=0)
          out$vv2 <- list(
              est = out1$est, ci = out1$ci, 
              cov = out1$ci[1] < 210 & out1$ci[2] > 210,
              err = out$est - 210
          )
      }else{
          out <- vector(mode = "list")
          
          out1 <- drtmle(Y=dat$Y, A=dat$A, W=dat$W, family=gaussian(),
                        a0=1,
                        libraryQ=c("SL.hal"),
                        libraryg=c("SL.hal"),
                        librarygr=c("SL.npreg"),
                        libraryQr=c("SL.npreg"),
                        reduction="univariate",
                        tolIC="default", 
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
          
          out$os <- list(
              est = out1$os$est, ci = out1$os$ci,
              cov = out1$os$ci[1] < 210 & out1$os$ci[2] > 210,
              err = out1$os$est - 210
          )
          
          out$os.dral <- list(
              est = out1$os.dral$est, ci = out1$os.dral$ci,
              cov = out1$os.dral$ci[1] < 210 & out1$os.dral$ci[2] > 210,
              err = out1$os.dral$est - 210
          )
          
          out$tmle <- list(
              est = out1$tmle$est, ci = out1$tmle$ci,
              cov = out1$tmle$ci[1] < 210 & out1$tmle$ci[2] > 210,
              err = out1$tmle$est - 210
          )
          
          out$tmle.dral <- list(
              est = out1$tmle.dral$est, ci = out1$tmle.dral$ci,
              cov = out1$tmle.dral$ci[1] < 210 & out1$tmle.dral$ci[2] > 210,
              err = out1$tmle.dral$est - 210
          )
          
      }
    save(out,file=paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData.tmp"))
    file.rename(paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData.tmp"),
                paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=",parm$est[i], ".RData"))
    print("file saved!")
  }
}

# merge job ###########################
if (args[1] == 'merge') {
    n <- c(500,2000)
    est <- c("parametric","vv2","adaptive")
    seed <- 1:500
    parm <- expand.grid(seed=seed, n=n)
    
    allOut <- NULL
    for(i in 1:nrow(parm)){
        out.parm <- get(load(paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=parametric.RData")))
        out.vv2 <- get(load(paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=vv2.RData")))
        out.adapt <- get(load(paste0("~/dral/out/KS_n=",parm$n[i],"_seed=",parm$seed[i],"_est=adaptive.RData")))
        
        allOut <- rbind(allOut,
                        unlist(out.parm),unlist(out.vv2),unlist(out.adapt))
    }
    save(allOut, file = "~/dral/out/KS_allOut.RData")
}
#--------------------
# make tables
#--------------------
setwd("~/Dropbox/Dissertation/one step dr inference/computerCode/computerOut")
library(xtable)

allOut6 <- get(load("newSim6_allOut.RData"))
allOut2 <- get(load("newSim2_allOut.RData"))

# get bias results
bias2 <- Reduce(cbind,by(allOut2[,grep("err",colnames(allOut2))], factor(allOut$n), colMeans, na.rm=TRUE))*10
bias6 <- Reduce(cbind,by(allOut6[,grep("err",colnames(allOut6))], factor(allOut$n), colMeans, na.rm=TRUE))*10

row.names(bias2) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
row.names(bias6) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
colnames(bias2) <- c("n=250","n=1000","n=5000")
colnames(bias6) <- c("n=250","n=1000","n=5000")
xtable(bias2,digits=3)
xtable(bias6,digits=3)

# root-n times bias
rootnbias2 <- Reduce(cbind,by(allOut2[,grep("err",colnames(allOut2))], factor(allOut$n), colMeans, na.rm=TRUE))*sqrt(c(250,1000,5000))
rootnbias6 <- Reduce(cbind,by(allOut6[,grep("err",colnames(allOut6))], factor(allOut$n), colMeans, na.rm=TRUE))*sqrt(c(250,1000,5000))

row.names(rootnbias2) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
row.names(rootnbias6) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
colnames(rootnbias2) <- c("250","1000","5000")
colnames(rootnbias6) <- c("250","1000","5000")
xtable(rootnbias2,digits=3)
xtable(rootnbias6,digits=3)

# get coverage results
cov2 <- Reduce(cbind,by(allOut2[,grep("cov",colnames(allOut2))], factor(allOut$n), colMeans, na.rm=TRUE))
cov6 <- Reduce(cbind,by(allOut6[,grep("cov",colnames(allOut6))], factor(allOut$n), colMeans, na.rm=TRUE))

row.names(cov2) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
row.names(cov6) <- c("Cao","Verm-1","Tan","Verm-2","OS","TMLE","TMLE-1","OS-1","TMLE-2","OS-2")
colnames(cov2) <- c("n=250","n=1000","n=5000")
colnames(cov6) <- c("n=250","n=1000","n=5000")

xtable(cov2, digits=3)
xtable(cov6, digits=3)

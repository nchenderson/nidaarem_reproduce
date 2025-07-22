
library(softImpute)
library(SQUAREM)
library(nidaarem)

setwd("~/Documents/nidaarem_reproduce")
source("SimulationCode/MatrixCompletion/MatrixCompletionFunctions.R")
load("SimulationCode/Data/movie_data.RData")

rm(VV)
n <- nrow(XX)
p <- ncol(XX)

ind <- which(XX != 0)

Y <- XX[ind]
XNA <- XX
XNA[-ind] <- NA

## Perform matrix completion using the fixed-point iteration (may take about 10 minutes)
set.seed(2356)
B.init <- matrix(0.0, nrow=n, ncol=p)

tols <- 1e-7
maxit <- 15000
log.lambda.seq <- seq(log(641/100), log(641), length.out=10)

lam.seq <- rev(exp(log.lambda.seq))
n.lambda <- length(lam.seq)

NIDAARAMNI <- RNIDAARAMNI <- DAAREMNI <- RDAARAMNI <- SQNI <- RENESTNI <- rep(0, n.lambda)
NESTNI <- PGDNI <- rep(0, n.lambda)
NIDAARAMObj <- RNIDAARAMObj <- DAAREMObj <- RDAARAMObj <- SQObj <- RENESTObj <- rep(0, n.lambda)
NESTObj <- PGDObj <- rep(0, n.lambda)
RNIDAARAMConv <- NIDAARAMConv <- RDAARAMConv <- DAAREMConv <- SQConv <- RENESTConv <- rep(0, n.lambda)
NESTConv <- PGDConv <- rep(0, n.lambda)
NIDAARAMTime <- RNIDAARAMTime <- DAAREMTime <- RDAARAMTime <- SQTime <- RENESTTime <- rep(0, n.lambda)
NESTTime <- PGDTime <- rep(0, n.lambda)
SIObj <- SITime <- rep(0, n.lambda)

for(k in 1:n.lambda) {
    ans.pgd.time <- system.time(ans.pgd <- fpiter(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                  A=Y, PA=XX, ind=ind, lambda=lam.seq[k], control=list(maxiter=maxit, tol=tols)))
    
    ans.nest.time <- system.time(ans.nest <- Nesterov(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, A=Y, PA=XX, 
                                                      ind=ind, lambda=lam.seq[k], control=list(maxiter=maxit, tol=tols)))
    
    ans.renest.time <- system.time(ans.renest <- Nesterov(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, restart=TRUE, 
                                                          A=Y, PA=XX, ind=ind, lambda=lam.seq[k], control=list(maxiter=maxit, tol=tols)))
    
    ans.daarem.time <- system.time(ans.daarem <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                          nesterov.init=FALSE, A=Y, PA=XX, ind=ind, lambda=lam.seq[k], 
                                                          control=list(maxiter=maxit, mon.tol=1, order=5, tol=tols)))
    
    ans.nidaarem.time <- system.time(ans.nidaarem <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                              A=Y, PA=XX, ind=ind, lambda=lam.seq[k], control=list(maxiter=maxit, order=5,
                                                                                                                   mon.tol=c(1,0), tol=tols)))
  
    ans.nidaarem.resid.time <- system.time(ans.nidaarem.resid <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                                         A=Y, PA=XX, ind=ind, lambda=lam.seq[k], control=list(tol=tols, 
                                                                        order=5, maxiter=maxit, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.sq.time <- system.time(ans.sq <- squarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=NegMatrixCompleteObj, A=Y, PA=XX, 
                                                ind=ind, lambda=lam.seq[k], control=list(maxiter=maxit, tol=tols)))
    
    if(k > 1) {
         ans.si.time <- system.time(ans.si <- softImpute(XNA, rank.max=942, lambda=lam.seq[k], thresh=1e-12, maxit=maxit, trace.it=TRUE,
                                                           warm.start=ans.si))
    } else {
         ans.si <- list()
         ans.si$u <- matrix(0.0, nrow=943, ncol=164)
         ans.si$v <- matrix(0.0, nrow=1682, ncol=164)
         ans.si$d <- rep(0.0, 164)
         ans.si.time <- system.time(ans.si <- softImpute(XNA, rank.max=942, lambda=lam.seq[k], thresh=1e-12, 
                                                         maxit=maxit, trace.it=TRUE, warm.start=ans.si))
    }
    
    ld <- length(ans.si$d)
    Xsi <- matrix(0, nrow=nrow(XNA), ncol=ncol(XNA))
    if(ld > 1) {
        for(h in 1:ld) {
            Xsi <- Xsi + ans.si$d[h]*outer(ans.si$u[,h], ans.si$v[,h])
        }
    } else if (ld == 1) {
        Xsi <- ans.si$d*outer(ans.si$u, ans.si$v)
    }
    si.obj <- MatrixCompleteObj(c(Xsi), lambda=lam.seq[k], A=Y, PA=XX, ind=ind)

    DAAREMNI[k] <- ans.daarem$fpevals
    NIDAARAMNI[k] <- ans.nidaarem$fpevals
    RNIDAARAMNI[k] <- ans.nidaarem.resid$fpevals
    SQNI[k] <- ans.sq$fpevals
    RENESTNI[k] <- ans.renest$fpevals
    NESTNI[k] <- ans.nest$fpevals
    PGDNI[k] <- ans.pgd$fpevals

    DAAREMObj[k] <- ans.daarem$value.objfn
    NIDAARAMObj[k] <- ans.nidaarem$value.objfn
    RNIDAARAMObj[k] <- ans.nidaarem.resid$value.objfn
    SQObj[k] <- (-1)*ans.sq$value.objfn
    RENESTObj[k] <- ans.renest$value.objfn
    NESTObj[k] <- ans.nest$value.objfn
    PGDObj[k] <- ans.pgd$value.objfn
    SIObj[k] <- si.obj
    
    DAAREMConv[k] <- ans.daarem$convergence
    NIDAARAMConv[k] <- ans.nidaarem$convergence
    RNIDAARAMConv[k] <- ans.nidaarem.resid$convergence
    SQConv[k] <- ans.sq$convergence
    RENESTConv[k] <- ans.renest$convergence
    NESTConv[k] <- ans.nest$convergence
    PGDConv[k] <- ans.pgd$convergence
    
    DAAREMTime[k] <- ans.daarem.time[3]
    NIDAARAMTime[k] <- ans.nidaarem.time[3]
    RNIDAARAMTime[k] <- ans.nidaarem.resid.time[3]
    SQTime[k] <- ans.sq.time[3]
    RENESTTime[k] <- ans.renest.time[3]
    NESTTime[k] <- ans.nest.time[3]
    PGDTime[k] <- ans.pgd.time[3]
    SITime[k] <- ans.si.time[3]
    
    B.init <- ans.nidaarem$par
    print(k)
}

save(DAAREMNI, NIDAARAMNI, RNIDAARAMNI, SQNI, RENESTNI, NESTNI, PGDNI,
     DAAREMObj, NIDAARAMObj, RNIDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, SIObj,
     DAAREMTime, NIDAARAMTime, RNIDAARAMTime, SQTime, RENESTTime, NESTTime, PGDTime, SITime,
     DAAREMConv, NIDAARAMConv, RNIDAARAMConv, SQConv, RENESTConv, NESTConv, PGDConv,
     file="SimulationResults/MatrixCompletion/MatrixCompleteWarmStarts.RData")







## Code to run simulations with penalized logistic regression on
## the "madelon" data. This computes the solution across
## all 10 lambda values using "warm starts". An initial 
## value of all zeros is used for the first lambda.

library(SQUAREM)
library(nidaarem)

setwd("~/Documents/nidaarem_reproduce")
source("SimulationCode/LogisticLasso/LogisticRegressionFunctions.R")


ni <- 200000  ## maximum of 200,000 iterations.
tols <- 1e-7
n.reps <- 1  ## Run 50 times each with a different random starting value

## Load madelon data and get responses y and scaled design matrix XX
data(madelon)
yy <- madelon$y
XX <- madelon$X
XX <- scale(XX)
rm(madelon)

n <- nrow(XX)
p <- ncol(XX)

# Subtract off mean if you want
for(h in 1:p) {
  XX[,h] <- XX[,h] - mean(XX[,h])
}

eps <- .01
K <- 10
ww <- (1000/n)*(yy - mean(yy)*(1 - mean(yy)))
Xtw <- crossprod(XX, ww)
log.lambda.max <- log(max(abs(Xtw)))
log.lambda.min <- log(eps) + log.lambda.max
log.lambdas <- seq(log.lambda.max, log.lambda.min, length.out=K)
lambda.seq <- exp(log.lambdas)
n.lambda <- length(lambda.seq)

XtX <- crossprod(XX, XX)
Xty <- as.vector(crossprod(XX, yy))
VV <- eigen(XtX)
LL <- max(VV$values)/4

NIDAAREMNI <- NIDAARAMNI <- DAAREMNI <- DAARAMNI <- SQNI <- RENESTNI <- matrix(0, n.lambda, n.reps)
NESTNI <- PGDNI <- MPENI <- matrix(0, n.lambda, n.reps)
NIDAAREMObj <- NIDAARAMObj <- DAAREMObj <- DAARAMObj <- SQObj <- RENESTObj <- matrix(0, n.lambda, n.reps)
NESTObj <- PGDObj <- MPEObj <- matrix(0, n.lambda, n.reps)
NIDAAREMConv <- NIDAARAMConv <- DAAREMConv <- DAARAMConv <- SQConv <- RENESTConv <- matrix(0, n.lambda, n.reps)
NESTConv <- PGDConv <- MPEConv <- matrix(0, n.lambda, n.reps)
NIDAAREMTime <- NIDAARAMTime <- DAAREMTime <- DAARAMTime <- SQTime <- RENESTTime <- matrix(0, n.lambda, n.reps)
NESTTime <- PGDTime <- MPETime <- matrix(0, n.lambda, n.reps)

for(k in 1:n.lambda) {
  lam <- lambda.seq[k]
  if(k==1) {
      beta.init <- rep(0, p)
  } else if(k > 1) {
      beta.init <- ans.daarem$par
  }
  for(j in 1:n.reps) {
    
    ans.pgd.time <- system.time(ans.pgd <- fpiter(par=beta.init, fixptfn=GDLogisticStep, objfn=LogisticObjFn, X=XX, Xty=Xty, 
                                                  lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols)))
    
    ans.nest.time<- system.time(ans.nest <- Nesterov(beta.init, fixptfn=GDLogisticStep, objfn=LogisticObjFn, restart=FALSE, 
                                                     X=XX, Xty=Xty, lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols)))
    
    ans.renest.time <- system.time(ans.renest <- Nesterov(beta.init, fixptfn=GDLogisticStep, objfn=LogisticObjFn, restart=TRUE, 
                                                          X=XX, Xty=Xty, lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols)))
    
    ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, family="binomial",
                                                              nesterov.init=FALSE, control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.daarem.time <- system.time(ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, family="binomial",
                                                              nesterov.init=FALSE, control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.nidaaram.time <- system.time(ans.nidaaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, family="binomial",
                                                                  nesterov.init=TRUE, control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.time <- system.time(ans.nidaarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, family="binomial",
                                                                  nesterov.init=TRUE, control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.sq.time <- system.time(ans.sq <- squarem(par=beta.init, fixptfn=LogisticUpdate, objfn=NegLogisticObjFn, X=XX, Xty=Xty, 
                                                 lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols)))
    
    ans.mpe.time <- system.time(ans.mpe <- squarem(par=beta.init, fixptfn=LogisticUpdate, objfn=NegLogisticObjFn, X=XX, Xty=Xty, 
                                                   lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols, K=3, method="MPE")))
    
    DAAREMNI[k,j] <- ans.daarem$fpevals
    DAARAMNI[k,j] <- ans.daaram$fpevals
    NIDAAREMNI[k,j] <- ans.nidaarem$fpevals
    NIDAARAMNI[k,j] <- ans.nidaaram$fpevals
    SQNI[k,j] <- ans.sq$fpevals
    RENESTNI[k,j] <- ans.renest$fpevals
    NESTNI[k,j] <- ans.nest$fpevals
    PGDNI[k,j] <- ans.pgd$fpevals
    MPENI[k,j] <- ans.mpe$fpevals
    
    DAAREMObj[k,j] <- ans.daarem$value.objfn
    DAARAMObj[k,j] <- ans.daaram$value.objfn
    NIDAAREMObj[k,j] <- ans.nidaarem$value.objfn
    NIDAARAMObj[k,j] <- ans.nidaaram$value.objfn
    SQObj[k,j] <- (-1)*ans.sq$value.objfn
    RENESTObj[k,j] <- ans.renest$value.objfn
    NESTObj[k,j] <- ans.nest$value.objfn
    PGDObj[k,j] <- ans.pgd$value.objfn
    MPEObj[k,j] <- (-1)*ans.mpe$value.objfn
    
    DAAREMConv[k,j] <- ans.daarem$convergence
    DAARAMConv[k,j] <- ans.daaram$convergence
    NIDAAREMConv[k,j] <- ans.nidaarem$convergence
    NIDAARAMConv[k,j] <- ans.nidaaram$convergence
    SQConv[k,j] <- ans.sq$convergence
    RENESTConv[k,j] <- ans.renest$convergence
    NESTConv[k,j] <- ans.nest$convergence
    PGDConv[k,j] <- ans.pgd$convergence
    MPEConv[k,j] <- ans.mpe$convergence
    
    DAAREMTime[k,j] <- ans.daarem.time[3]
    DAARAMTime[k,j] <- ans.daaram.time[3]
    NIDAAREMTime[k,j] <- ans.nidaarem.time[3]
    NIDAARAMTime[k,j] <- ans.nidaaram.time[3]
    SQTime[k,j] <- ans.sq.time[3]
    RENESTTime[k,j] <- ans.renest.time[3]
    NESTTime[k,j] <- ans.nest.time[3]
    PGDTime[k,j] <- ans.pgd.time[3]
    MPETime[k,j] <- ans.mpe.time[3]
    
    print(c(k,j))
  }
  #if(k==1) {
  ##  fname <- "SimulationResults/LogisticLasso/MadelonFixedLambda1.RData"
  #} else if(k==2){
  #  fname <- "SimulationResults/LogisticLasso/MadelonFixedLambda2.RData"
  #}
  #save(DAAREMNI, DAARAMNI, NIDAAREMNI, NIDAARAMNI, SQNI, RENESTNI, NESTNI, PGDNI,
  #    DAAREMObj, DAARAMObj, NIDAAREMObj, NIDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, SIObj,
  #    DAAREMTime, DAARAMTime, NIDAAREMTime, NIDAARAMTime, SQTime, RENESTTime, NESTTime, PGDTime, SITime,
  #   DAAREMConv, DAARAMConv, NIDAAREMConv, NIDAARAMConv, SQConv, RENESTConv, NESTConv, PGDConv,
  #   file=fname)
}

## Save results:
fname <- "SimulationResults/LogisticLasso/MadelonFullPath.RData"
save(DAAREMNI, DAARAMNI, NIDAAREMNI, NIDAARAMNI, SQNI, RENESTNI, NESTNI, PGDNI, MPENI,
     DAAREMObj, DAARAMObj, NIDAAREMObj, NIDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, MPEObj,
     DAAREMTime, DAARAMTime, NIDAAREMTime, NIDAARAMTime, SQTime, RENESTTime, NESTTime, PGDTime, MPETime,
     DAAREMConv, DAARAMConv, NIDAAREMConv, NIDAARAMConv, SQConv, RENESTConv, NESTConv, PGDConv, MPEConv,
     file=fname)

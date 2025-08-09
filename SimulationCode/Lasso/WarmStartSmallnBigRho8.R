
### Simulations with l1-penalized regression with small n and big p
### n = 100 and p = 10000
### Covariates X_{ij} are generated so that cor(X_{ij}, X_{ik}) = rho^|j-k|.
###   (in these simulations rho = 0.95)
### Simulations are performed for a sequence of penalty parameters with "warm starts".
###
### Methods compared: (1) Regular Proximal Gradient (PG), (2) FISTA, (3) SQUAREM
###                   (4) DAAREM, (5) NIDAAREM, (6) Subsetted NIDAAREM,
###                   (7) FISTA with restarts

library(nidaarem)
library(SQUAREM)
setwd("~/Documents/nidaarem_reproduce")
source("SimulationCode/Lasso/LassoFunctions.R")


n <- 100
p <- 10000

ni <- 500000  ## maximum of 500,000 iterations.
tols <- 1e-7

rho <- 0.95
####################################################################
## Do a preliminary simulation to get a sensible lambda sequence.
true.beta <- rt(p, df=3) * rbinom(p, 1, 0.2)
XX <- matrix(0.0, nrow=n, ncol=p)
for(k in 1:n) {
  XX[k,] <- arima.sim(n = p, list(ar = rho), innov=rnorm(p))
}
XX <- scale(XX)
eps <- 0.01
K <- 20
n.lambda <- K


#####################################################################
#####################################################################
nreps <- 5

NIDAARAMNI <- NIDAAREMNI <- RNIDAARAMNI <- DAARAMNI <- RDAARAMNI <- SQNI <- RENESTNI <- matrix(0, nrow=nreps, ncol=n.lambda)
DAAREMNI <- NESTNI <- PGDNI <- MPENI <- SNIDAARAMNI <- SDAARAMNI <- matrix(0, nrow=nreps, ncol=n.lambda)
NIDAARAMObj <- NIDAAREMObj <- RNIDAARAMObj <- DAARAMObj <- RDAARAMObj <- SQObj <- RENESTObj <- matrix(0, nrow=nreps, ncol=n.lambda)
DAAREMObj <- NESTObj <- PGDObj <- MPEObj <- SNIDAARAMObj <- SDAARAMObj <- matrix(0, nrow=nreps, ncol=n.lambda)
RNIDAARAMConv <- NIDAAREMConv <- NIDAARAMConv <- RDAARAMConv <- DAARAMConv <- SQConv <- RENESTConv <- matrix(0, nrow=nreps, ncol=n.lambda)
DAAREMConv <- NESTConv <- PGDConv <- MPEConv <- SNIDAARAMConv <- SDAARAMConv <- matrix(0, nrow=nreps, ncol=n.lambda)
NIDAARAMTime <- NIDAAREMTime <- RNIDAARAMTime <- DAARAMTime <- RDAARAMTime <- SQTime <- RENESTTime <- matrix(0, nrow=nreps, ncol=n.lambda)
DAAREMTime <- NESTTime <- PGDTime <- MPETime <- SNIDAARAMTime <- SDAARAMTime <- matrix(0, nrow=nreps, ncol=n.lambda)
for(j in 1:nreps) {
  yy <- XX %*% true.beta + rnorm(n)
  yy <- yy - mean(yy)
  
  Xty <- crossprod(XX, yy)
  log.lambda.max <- log(max(abs(Xty)))
  log.lambda.min <- log(eps) + log.lambda.max
  log.lambdas <- seq(log.lambda.max, log.lambda.min, length.out=K)
  lambda.seq <- exp(log.lambdas)
  beta.init <- rep(0.0, p)
  Lmax <- norm(XX, "2")^2
  stp <- 1/Lmax
  for(k in 1:K) {
    
    lam <- lambda.seq[k]
    ans.pgd.time <- system.time(ans.pgd <- fpiter(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                  control=list(maxiter=ni, tol=tols)))
    
    ans.nest.time <- system.time(ans.nest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                      restart=FALSE, control=list(maxiter=ni, tol=tols)))
    
    ans.renest.time <- system.time(ans.renest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                          restart=TRUE, control=list(maxiter=ni, tol=tols)))
    
    ans.daarem.time <- system.time(ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.resid.time <- system.time(ans.nidaarem.resid <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.nidaaram.time <- system.time(ans.nidaaram <- daarem.lasso(par=rep(0, p), X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                  control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.time <- system.time(ans.nidaarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                  control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.sq.time <- system.time(ans.sq <- squarem(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                 control=list(maxiter=ni, tol=tols)))
    
    ans.snidaarem.time <- system.time(ans.snidaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                     sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.sdaarem.time <- system.time(ans.sdaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                                 sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    beta.init <- ans.renest$par
    ## Need to add subsetted lasso as well
    
    DAAREMNI[j,k] <- ans.daarem$fpevals
    NIDAAREMNI[j,k] <- ans.nidaarem$fpevals
    DAARAMNI[j,k] <- ans.daaram$fpevals
    NIDAARAMNI[j,k] <- ans.nidaaram$fpevals
    RNIDAARAMNI[j,k] <- ans.nidaarem.resid$fpevals
    SQNI[j,k] <- ans.sq$fpevals
    RENESTNI[j,k] <- ans.renest$fpevals
    NESTNI[j,k] <- ans.nest$fpevals
    PGDNI[j,k] <- ans.pgd$fpevals
    SNIDAARAMNI[j,k] <- ans.snidaarem$fpevals
    SDAARAMNI[j,k] <- ans.sdaarem$fpevals
    
    DAAREMObj[j,k] <- ans.daarem$value.objfn
    NIDAAREMObj[j,k] <- ans.nidaarem$value.objfn
    DAARAMObj[j,k] <- ans.daaram$value.objfn
    NIDAARAMObj[j,k] <- ans.nidaaram$value.objfn
    RNIDAARAMObj[j,k] <- ans.nidaarem.resid$value.objfn
    SQObj[j,k] <- (-1)*ans.sq$value.objfn
    RENESTObj[j,k] <- ans.renest$value.objfn
    NESTObj[j,k] <- ans.nest$value.objfn
    PGDObj[j,k] <- ans.pgd$value.objfn
    SNIDAARAMObj[j,k] <- ans.snidaarem$value.objfn
    SDAARAMObj[j,k] <- ans.sdaarem$value.objfn
    
    DAAREMConv[j,k] <- ans.daarem$convergence
    NIDAAREMConv[j,k] <- ans.nidaarem$convergence
    DAARAMConv[j,k] <- ans.daaram$convergence
    NIDAARAMConv[j,k] <- ans.nidaaram$convergence
    RNIDAARAMConv[j,k] <- ans.nidaarem.resid$convergence
    SQConv[j,k] <- ans.sq$convergence
    RENESTConv[j,k] <- ans.renest$convergence
    NESTConv[j,k] <- ans.nest$convergence
    PGDConv[j,k] <- ans.pgd$convergence
    SNIDAARAMConv[j,k] <- ans.snidaarem$convergence
    SDAARAMConv[j,k] <- ans.sdaarem$convergence
    
    DAAREMTime[j,k] <- ans.daarem.time[3]
    NIDAAREMTime[j,k] <- ans.nidaarem.time[3]
    DAARAMTime[j,k] <- ans.daaram.time[3]
    NIDAARAMTime[j,k] <- ans.nidaaram.time[3]
    RNIDAARAMTime[j,k] <- ans.nidaarem.resid.time[3]
    SQTime[j,k] <- ans.sq.time[3]
    RENESTTime[j,k] <- ans.renest.time[3]
    NESTTime[j,k] <- ans.nest.time[3]
    PGDTime[j,k] <- ans.pgd.time[3]
    SNIDAARAMTime[j,k] <- ans.snidaarem.time[3]
    SDAARAMTime[j,k] <- ans.sdaarem.time[3]
    
    print(c(j,k))
  }
}

fname <- "SimulationResults/Lasso/LassoRho95WarmStarts.RData"
save(NIDAARAMNI, RNIDAARAMNI, DAAREMNI, DAARAMNI, NIDAAREMNI, SQNI, RENESTNI, 
     NESTNI, PGDNI, SNIDAARAMNI, SDAARAMNI, NIDAARAMObj, RNIDAARAMObj, DAAREMObj, DAARAMObj,
     NIDAAREMObj, SQObj, RENESTObj, NESTObj, PGDObj, SNIDAARAMObj, SDAARAMObj, 
     RNIDAARAMConv, NIDAARAMConv, DAAREMConv, DAARAMConv, SQConv, RENESTConv, NIDAAREMConv,
     NESTConv, PGDConv, SNIDAARAMConv, SDAARAMConv, NIDAARAMTime, RNIDAARAMTime, 
     DAAREMTime, DAARAMTime, NIDAAREMTime, SQTime, RENESTTime, NESTTime, PGDTime, SNIDAARAMTime,
     SDAARAMTime, file=fname)



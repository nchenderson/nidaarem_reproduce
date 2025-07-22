### Simulations with l1-penalized regression with small n and big p
### n = 100 and p = 10000
### Covariates X_{ij} are generated so that cor(X_{ij}, X_{ik}) = rho^|j-k|.
###   (in these simulations rho = 0.8)
### Simulations are performed for a sequence of penalty parameters with "warm starts".
###
### Methods compared: (1) Regular Proximal Gradient (PG), (2) FISTA, (3) SQUAREM
###                   (4) DAAREM, (5) NIDAAREM (2 versions), (6) Subsetted NIDAAREM,
###                   (7) FISTA with restarts

library(nidaarem)
library(SQUAREM)
setwd("~/Documents/nidaarem_reproduce")
source("SimulationCode/Lasso/LassoFunctions.R")

set.seed(314127)
n <- 10000
p <- 100

## Use 500000 in full simulations
ni <- 500000
tols <- 1e-7
rho <- 0.8

#####################################################################
nreps <- 50 ## Do nreps = 50 in full simulations

lambda.vals <- rep(c(500, 50), each=3)
n.vals <- rep(c(10000, 2000, 100), 2)
p.vals <- rep(c(100, 500, 10000), 2)
SimDesign <- cbind(lambda.vals, n.vals, p.vals)
colnames(SimDesign) <- c("lambda", "n", "p")

NIDAARAMNI <- RNIDAARAMNI <- DAAREMNI <- DAARAMNI <- SQNI <- RENESTNI <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAAREMNI <- NESTNI <- PGDNI <- SNIDAARAMNI <- SDAARAMNI <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAARAMObj <- RNIDAARAMObj <- DAAREMObj <- DAARAMObj <- SQObj <- RENESTObj <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAAREMObj <- NESTObj <- PGDObj <- SNIDAARAMObj <- SDAARAMObj <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
RNIDAARAMConv <- NIDAARAMConv <- DAAREMConv <- DAARAMConv <- SQConv <- RENESTConv <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAAREMConv <- NESTConv <- PGDConv <- SNIDAARAMConv <- SDAARAMConv <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAARAMTime <- RNIDAARAMTime <- DAAREMTime <- DAARAMTime <- SQTime <- RENESTTime <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
NIDAAREMTime <- NESTTime <- PGDTime <- SNIDAARAMTime <- SDAARAMTime <- matrix(0, nrow=nreps, ncol=nrow(SimDesign))
for(k in 1:nrow(SimDesign)) {
  ## Generate true.beta to be used for all nreps reps
  lam <- lambda.vals[k]
  n <- n.vals[k]
  p <- p.vals[k]
  true.beta <- rt(p, df=3) * rbinom(p, 1, 0.2)
  beta.init <- rep(0.0, p)
  for(j in 1:nreps) {
    
    XX <- matrix(0.0, nrow=n, ncol=p)
    for(h in 1:n) {
      XX[h,] <- arima.sim(n = p, list(ar = rho), innov=rnorm(p))
    }
    XX <- scale(XX)
    yy <- XX %*% true.beta + rnorm(n)
    yy <- yy - mean(yy)
    Lmax <- norm(XX, "2")^2
    stp <- 1/Lmax
 
    ans.pgd.time <- system.time(ans.pgd <- fpiter(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                  control=list(maxiter=ni, tol=tols)))
    
    ans.nest.time <- system.time(ans.nest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                      restart=FALSE, control=list(maxiter=ni, tol=tols)))
    
    ans.renest.time <- system.time(ans.renest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                          restart=TRUE, control=list(maxiter=ni, tol=tols)))
    
    #ans.daaram.resid.time <- system.time(ans.daaram.resid <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
    #                                                                      control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.daarem.time <- system.time(ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
  
    ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.resid.time <- system.time(ans.nidaarem.resid <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.nidaaram.time <- system.time(ans.nidaaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                  control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.time <- system.time(ans.nidaarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                  control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.sq.time <- system.time(ans.sq <- squarem(par=beta.init, fixptfn=GDLassoStep, objfn=NegLassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                 control=list(maxiter=ni, tol=tols)))
    
    ans.snidaarem.time <- system.time(ans.snidaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                     sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ans.sdaarem.time <- system.time(ans.sdaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                                 sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=1)))
    
    ## Need to add subsetted lasso as well
    
    DAAREMNI[j,k] <- ans.daarem$fpevals
    DAARAMNI[j,k] <- ans.daaram$fpevals
    NIDAAREMNI[j,k] <- ans.nidaarem$fpevals
    NIDAARAMNI[j,k] <- ans.nidaaram$fpevals
    RNIDAARAMNI[j,k] <- ans.nidaarem.resid$fpevals
    SQNI[j,k] <- ans.sq$fpevals
    RENESTNI[j,k] <- ans.renest$fpevals
    NESTNI[j,k] <- ans.nest$fpevals
    PGDNI[j,k] <- ans.pgd$fpevals
    SNIDAARAMNI[j,k] <- ans.snidaarem$fpevals
    SDAARAMNI[j,k] <- ans.sdaarem$fpevals
    
    DAAREMObj[j,k] <- ans.daarem$value.objfn
    DAARAMObj[j,k] <- ans.daarem$value.objfn
    NIDAAREMObj[j,k] <- ans.nidaarem$value.objfn
    NIDAARAMObj[j,k] <- ans.nidaaram$value.objfn
    RNIDAARAMObj[j,k] <- ans.nidaarem.resid$value.objfn
    SQObj[j,k] <- (-1)*ans.sq$value.objfn
    RENESTObj[j,k] <- ans.renest$value.objfn
    NESTObj[j,k] <- ans.nest$value.objfn
    PGDObj[j,k] <- ans.pgd$value.objfn
    SNIDAARAMObj[j,k] <- ans.snidaarem$value.objfn
    SDAARAMObj[j,k] <- ans.sdaarem$value.objfn
    
    
    DAAREMConv[j,k] <- ans.daarem$convergence
    DAARAMConv[j,k] <- ans.daarem$convergence
    NIDAAREMConv[j,k] <- ans.nidaarem$convergence
    NIDAARAMConv[j,k] <- ans.nidaaram$convergence
    RNIDAARAMConv[j,k] <- ans.nidaarem.resid$convergence
    SQConv[j,k] <- ans.sq$convergence
    RENESTConv[j,k] <- ans.renest$convergence
    NESTConv[j,k] <- ans.nest$convergence
    PGDConv[j,k] <- ans.pgd$convergence
    SNIDAARAMConv[j,k] <- ans.snidaarem$convergence
    SDAARAMConv[j,k] <- ans.sdaarem$convergence
    
    DAAREMTime[j,k] <- ans.daarem.time[3]
    DAARAMTime[j,k] <- ans.daarem.time[3]
    NIDAAREMTime[j,k] <- ans.nidaarem.time[3]
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

fname <- "SimulationResults/Lasso/LassoRho8.RData"
save(SimDesign, NIDAARAMNI, RNIDAARAMNI, DAAREMNI, DAARAMNI, NIDAAREMNI, SQNI, RENESTNI, 
     NESTNI, PGDNI, SNIDAARAMNI, SDAARAMNI, NIDAARAMObj, RNIDAARAMObj, DAAREMObj, DAARAMObj,
     NIDAAREMObj, SQObj, RENESTObj, NESTObj, PGDObj, SNIDAARAMObj, SDAARAMObj, 
     RNIDAARAMConv, NIDAARAMConv, DAAREMConv, DAARAMConv, SQConv, RENESTConv, NIDAAREMConv,
     NESTConv, PGDConv, SNIDAARAMConv, SDAARAMConv, NIDAARAMTime, RNIDAARAMTime, 
     DAAREMTime, DAARAMTime, NIDAAREMTime, SQTime, RENESTTime, NESTTime, PGDTime, SNIDAARAMTime,
     SDAARAMTime, file=fname)







### Simulations with l1-penalized regression with small n and big p
### n = 100 and p = 10000
### Covariates X_{ij} are generated so that cor(X_{ij}, X_{ik}) = rho^|j-k|.
###   (in these simulations rho = 0.6)
### Simulations are performed for a sequence of penalty parameters with "warm starts".
###
### Methods compared: (1) Regular Proximal Gradient (PG), (2) FISTA, (3) SQUAREM
###                   (4) DAAREM, (5) NIDAAREM, (6) Subsetted NIDAAREM,
###                   (7) FISTA with restarts
###   (should we also compare DAAREM with and without objective function)
###   (maybe for another simulation.)

library(nidaarem)
library(SQUAREM)

set.seed(51367) ## the 2nd and 19th runs are quite long for the DAAREM methods.
set.seed(314127)
n <- 2000
p <- 500

ni <- 3000000
#ni <- 100000
#ni <- 5000
tols <- 1e-7

rho <- 0.6
####################################################################
## Do a preliminary simulation to get a sensible lambda sequence.
true.beta <- rt(p, df=3) * rbinom(p, 1, 0.2)

#####################################################################
nreps <- 100
beta.init <- rep(0.0, p)

### If we are going to try 2 values of lambda, use lambda = 1,000 and lambda = 100
lambda.seq <- c(1000, 50) ## (1000,100) works well enough for rho=0.9; also (2000, 100); reasonable for (2000, 200); works well for (4000, 300)
n.lambda <- 2

NIDAARAMNI <- RNIDAARAMNI <- DAARAMNI <- RDAARAMNI <- SQNI <- RENESTNI <- matrix(0, nrow=nreps, ncol=n.lambda)
NESTNI <- PGDNI <- MPENI <- SNIDAARAMNI <- SDAARAMNI <- matrix(0, nrow=nreps, ncol=n.lambda)
NIDAARAMObj <- RNIDAARAMObj <- DAARAMObj <- RDAARAMObj <- SQObj <- RENESTObj <- matrix(0, nrow=nreps, ncol=n.lambda)
NESTObj <- PGDObj <- MPEObj <- SNIDAARAMObj <- SDAARAMObj <- matrix(0, nrow=nreps, ncol=n.lambda)
RNIDAARAMConv <- NIDAARAMConv <- RDAARAMConv <- DAARAMConv <- SQConv <- RENESTConv <- matrix(0, nrow=nreps, ncol=n.lambda)
NESTConv <- PGDConv <- MPEConv <- SNIDAARAMConv <- SDAARAMConv <- matrix(0, nrow=nreps, ncol=n.lambda)
NIDAARAMTime <- RNIDAARAMTime <- DAARAMTime <- RDAARAMTime <- SQTime <- RENESTTime <- matrix(0, nrow=nreps, ncol=n.lambda)
NESTTime <- PGDTime <- MPETime <- SNIDAARAMTime <- SDAARAMTime <- matrix(0, nrow=nreps, ncol=n.lambda)
for(j in 1:nreps) {
  for(k in 1:n.lambda) {
    beta.init <- rep(0.0, p)
    lam <- lambda.seq[k]
    XX <- matrix(0.0, nrow=n, ncol=p)
    for(h in 1:n) {
      XX[h,] <- arima.sim(n = p, list(ar = rho), innov=rnorm(p))
    }
    XX <- scale(XX)
    yy <- XX %*% true.beta + rnorm(n)
    yy <- yy - mean(yy)
    # stp <- FindStepLengthGaussian(XX, yy, lam)
    Lmax <- norm(XX, "2")^2
    stp <- 1/Lmax
    
    #eps <- 0.01
    #K <- 50
    #Xty <- crossprod(XX, yy)
    #log.lambda.max <- log(max(abs(Xty)))
    #log.lambda.min <- log(eps) + log.lambda.max
    #log.lambdas <- seq(log.lambda.max, log.lambda.min, length.out=K)
    #lambda.seq <- exp(log.lambdas)
    
    ans.pgd.time <- system.time(ans.pgd <- fpiter(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                  control=list(maxiter=ni, tol=tols)))
    
    ans.nest.time <- system.time(ans.nest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                      restart=FALSE, control=list(maxiter=ni, tol=tols)))
    
    ans.renest.time <- system.time(ans.renest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
                                                          restart=TRUE, control=list(maxiter=ni, tol=tols)))
    
    ans.daaram.resid.time <- system.time(ans.daaram.resid <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                                          control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.nidaarem.resid.time <- system.time(ans.nidaarem.resid <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                              control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(.95,0), objfn.check=FALSE)))
    
    ans.nidaarem.time <- system.time(ans.nidaarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                  control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.sq.time <- system.time(ans.sq <- squarem(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp, 
                                                 control=list(maxiter=ni, tol=tols)))
    
    #ans.mpe.time <- system.time(ans.mpe <- squarem(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=stp,
    #                                               control=list(maxiter=ni, tol=tols, K=3, method="MPE")))
    
    ans.snidaarem.time <- system.time(ans.snidaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=TRUE, 
                                                                     sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ans.sdaarem.time <- system.time(ans.sdaarem <- sdaarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=stp, nesterov.init=FALSE, 
                                                                 sub.type="threshold", control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))
    
    ## Need to add subsetted lasso as well
    
    DAARAMNI[j,k] <- ans.daaram$fpevals
    NIDAARAMNI[j,k] <- ans.nidaarem$fpevals
    RNIDAARAMNI[j,k] <- ans.nidaarem.resid$fpevals
    SQNI[j,k] <- ans.sq$fpevals
    RENESTNI[j,k] <- ans.renest$fpevals
    NESTNI[j,k] <- ans.nest$fpevals
    PGDNI[j,k] <- ans.pgd$fpevals
    #MPENI[j,k] <- ans.mpe$fpevals
    SNIDAARAMNI[j,k] <- ans.snidaarem$fpevals
    SDAARAMNI[j,k] <- ans.sdaarem$fpevals
    
    DAARAMObj[j,k] <- ans.daaram$value.objfn
    NIDAARAMObj[j,k] <- ans.nidaarem$value.objfn
    RNIDAARAMObj[j,k] <- ans.nidaarem.resid$value.objfn
    SQObj[j,k] <- (-1)*ans.sq$value.objfn
    RENESTObj[j,k] <- ans.renest$value.objfn
    NESTObj[j,k] <- ans.nest$value.objfn
    PGDObj[j,k] <- ans.pgd$value.objfn
    #MPEObj[j,k] <- (-1)*ans.mpe$value.objfn
    SNIDAARAMObj[j,k] <- ans.snidaarem$value.objfn
    SDAARAMObj[j,k] <- ans.sdaarem$value.objfn
    
    
    DAARAMConv[j,k] <- ans.daaram$convergence
    NIDAARAMConv[j,k] <- ans.nidaarem$convergence
    RNIDAARAMConv[j,k] <- ans.nidaarem.resid$convergence
    SQConv[j,k] <- ans.sq$convergence
    RENESTConv[j,k] <- ans.renest$convergence
    NESTConv[j,k] <- ans.nest$convergence
    PGDConv[j,k] <- ans.pgd$convergence
    #MPEConv[j,k] <- ans.mpe$convergence
    SNIDAARAMConv[j,k] <- ans.snidaarem$convergence
    SDAARAMConv[j,k] <- ans.sdaarem$convergence
    
    DAARAMTime[j,k] <- ans.daaram.time[3]
    NIDAARAMTime[j,k] <- ans.nidaarem.time[3]
    RNIDAARAMTime[j,k] <- ans.nidaarem.resid.time[3]
    SQTime[j,k] <- ans.sq.time[3]
    RENESTTime[j,k] <- ans.renest.time[3]
    NESTTime[j,k] <- ans.nest.time[3]
    PGDTime[j,k] <- ans.pgd.time[3]
    #MPETime[j,k] <- ans.mpe.time[3]
    SNIDAARAMTime[j,k] <- ans.snidaarem.time[3]
    SDAARAMTime[j,k] <- ans.sdaarem.time[3]
    
    print(c(j,k))
  }
}

save(NIDAARAMNI, RNIDAARAMNI, DAARAMNI, RDAARAMNI, SQNI, RENESTNI, 
     NESTNI, PGDNI, MPENI, SNIDAARAMNI, SDAARAMNI, NIDAARAMObj, RNIDAARAMObj, DAARAMObj, 
     RDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, MPEObj, SNIDAARAMObj, SDAARAMObj, 
     RNIDAARAMConv, NIDAARAMConv, RDAARAMConv, DAARAMConv, SQConv, RENESTConv, 
     NESTConv, PGDConv, MPEConv, SNIDAARAMConv, SDAARAMConv, 
     NIDAARAMTime, RNIDAARAMTime, DAARAMTime, RDAARAMTime, SQTime, RENESTTime,
     NESTTime, PGDTime, MPETime, SNIDAARAMTime, SDAARAMTime,
     file="~/Documents/NIDAAREM/SimulationResults/LassoMednMedp2.RData")






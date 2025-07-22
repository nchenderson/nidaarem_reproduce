

NegMatrixCompleteObj <- function(par, lambda, A, PA, ind) {
  B <- matrix(par, nrow=nrow(PA), ncol=ncol(PA))
  svdB <- svd(B)
  
  ans <- sum((par[ind] - A)*(par[ind] - A))/2 + lambda*sum(svdB$d)
  return(-ans)
}

### Look at softImpute R package

library(softImpute)
library(SQUAREM)
library(nidaarem)

load("~/Documents/NIDAAREM/Data/movie_data.RData")
rm(VV)
n <- nrow(XX)
p <- ncol(XX)

ind <- which(XX != 0)

Y <- XX[ind]
XNA <- XX
XNA[-ind] <- NA

## Perform matrix completion using the fixed-point iteration (may take about 10 minutes)
mu.B <- mean(Y)
sd.B <- sd(Y)

tols <- 1e-7
maxit <- 15000
#maxit.prime <- 12000


lambda.star <- 20
nreps <- 10

NIDAARAMNI <- RNIDAARAMNI <- DAARAMNI <- RDAARAMNI <- SQNI <- RENESTNI <- rep(0, nreps)
NESTNI <- PGDNI <- MPENI <- rep(0, nreps)
NIDAARAMObj <- RNIDAARAMObj <- DAARAMObj <- RDAARAMObj <- SQObj <- RENESTObj <- rep(0, nreps)
NESTObj <- PGDObj <- MPEObj <- rep(0, nreps)
RNIDAARAMConv <- NIDAARAMConv <- RDAARAMConv <- DAARAMConv <- SQConv <- RENESTConv <- rep(0, nreps)
NESTConv <- PGDConv <- MPEConv <- rep(0, nreps)
NIDAARAMTime <- RNIDAARAMTime <- DAARAMTime <- RDAARAMTime <- SQTime <- RENESTTime <- rep(0, nreps)
NESTTime <- PGDTime <- MPETime <- rep(0, nreps)
SIObj <- SITime <- rep(0, nreps)



for(h in 1:nreps) {
  B.init <- matrix(rnorm(n*p, mean=mu.B, sd=sd.B), nrow=n, ncol=p)
  
  ans.pgd.time <- system.time(ans.pgd <- fpiter(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                A=Y, PA=XX, ind=ind, lambda=lambda.star, control=list(maxiter=maxit, tol=tols)))
  
  ans.nest.time <- system.time(ans.nest <- Nesterov(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, A=Y, PA=XX, 
                                                    ind=ind, lambda=lambda.star, control=list(maxiter=maxit, tol=tols)))
  
  ans.renest.time <- system.time(ans.renest <- Nesterov(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, restart=TRUE, 
                                                        A=Y, PA=XX, ind=ind, lambda=lambda.star, control=list(maxiter=maxit, tol=tols)))
  
  ans.daaram.time <- system.time(ans.daaram <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, nesterov.init=FALSE, 
                                                        A=Y, PA=XX, ind=ind, lambda=lambda.star, control=list(maxiter=maxit, order=5, 
                                                                                                              mon.tol=c(1,0), tol=tols)))
  
  ans.nidaarem.time <- system.time(ans.nidaarem <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                            A=Y, PA=XX, ind=ind, lambda=lambda.star, control=list(maxiter=maxit, order=5,
                                                                                                                  mon.tol=c(1,0), tol=tols)))
  
  ans.nidaarem.resid.time <- system.time(ans.nidaarem.resid <- nidaarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                                                                        A=Y, PA=XX, ind=ind, lambda=lambda.star, control=list(tol=tols, 
                                                                                                                              order=5, maxiter=maxit, mon.tol=c(.95,0), objfn.check=FALSE)))
  
  ans.sq.time <- system.time(ans.sq <- squarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=NegMatrixCompleteObj, A=Y, PA=XX, 
                                               ind=ind, lambda=lambda.star, control=list(maxiter=maxit, tol=tols)))
  
  #ans.mpe.time <- system.time(ans.mpe <- squarem(par=c(B.init), fixptfn=MatrixCompleteUpdate, objfn=NegMatrixCompleteObj, A=Y, PA=XX, 
  #                                               ind=ind, lambda=lambda.star, control=list(maxiter=maxit, tol=tols, K=3, method="MPE")))
  
  ans.si.time <- system.time(ans.si <- softImpute(XNA, rank.max=942, lambda=lambda.star, thresh=1e-12, maxit=maxit, trace.it=TRUE))
  
  ld <- length(ans.si$d)
  Xsi <- matrix(0, nrow=nrow(XNA), ncol=ncol(XNA))
  if(ld > 1) {
    for(l in 1:ld) {
      Xsi <- Xsi + ans.si$d[l]*outer(ans.si$u[,l], ans.si$v[,l])
    }
  } else if (ld == 1) {
    Xsi <- ans.si$d*outer(ans.si$u, ans.si$v)
  }
  #Xsin <- ans.si$u%*%diag(ans.si$d)%*%t(ans.si$v)
  si.obj <- MatrixCompleteObj(c(Xsi), lambda=lambda.star, A=Y, PA=XX, ind=ind)
  
  DAARAMNI[h] <- ans.daaram$fpevals
  NIDAARAMNI[h] <- ans.nidaarem$fpevals
  RNIDAARAMNI[h] <- ans.nidaarem.resid$fpevals
  SQNI[h] <- ans.sq$fpevals
  RENESTNI[h] <- ans.renest$fpevals
  NESTNI[h] <- ans.nest$fpevals
  PGDNI[h] <- ans.pgd$fpevals
  #MPENI[h] <- ans.mpe$fpevals
  
  DAARAMObj[h] <- ans.daaram$value.objfn
  NIDAARAMObj[h] <- ans.nidaarem$value.objfn
  RNIDAARAMObj[h] <- ans.nidaarem.resid$value.objfn
  SQObj[h] <- (-1)*ans.sq$value.objfn
  RENESTObj[h] <- ans.renest$value.objfn
  NESTObj[h] <- ans.nest$value.objfn
  PGDObj[h] <- ans.pgd$value.objfn
  #MPEObj[h] <- (-1)*ans.mpe$value.objfn
  SIObj[h] <- si.obj
  
  DAARAMConv[h] <- ans.daaram$convergence
  NIDAARAMConv[h] <- ans.nidaarem$convergence
  RNIDAARAMConv[h] <- ans.nidaarem.resid$convergence
  SQConv[h] <- ans.sq$convergence
  RENESTConv[h] <- ans.renest$convergence
  NESTConv[h] <- ans.nest$convergence
  PGDConv[h] <- ans.pgd$convergence
 # MPEConv[h] <- ans.mpe$convergence
  
  DAARAMTime[h] <- ans.daaram.time[3]
  NIDAARAMTime[h] <- ans.nidaarem.time[3]
  RNIDAARAMTime[h] <- ans.nidaarem.resid.time[3]
  SQTime[h] <- ans.sq.time[3]
  RENESTTime[h] <- ans.renest.time[3]
  NESTTime[h] <- ans.nest.time[3]
  PGDTime[h] <- ans.pgd.time[3].
  #MPETime[h] <- ans.mpe.time[3]
  SITime[h] <- ans.si.time[3]
  
  print(h)
}

cbind(DAARAMNI, NIDAARAMNI, RNIDAARAMNI, SQNI, RENESTNI, NESTNI, PGDNI, MPENI)
cbind(DAARAMObj, NIDAARAMObj, RNIDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, MPEObj, SIObj)
TT <- cbind(DAARAMTime, NIDAARAMTime, RNIDAARAMTime, SQTime, RENESTTime, NESTTime, PGDTime, MPETime, SITime)

#save(DAARAMNI, NIDAARAMNI, RNIDAARAMNI, SQNI, RENESTNI, NESTNI, PGDNI, MPENI,
#     DAARAMObj, NIDAARAMObj, RNIDAARAMObj, SQObj, RENESTObj, NESTObj, PGDObj, MPEObj, SIObj,
#     DAARAMTime, NIDAARAMTime, RNIDAARAMTime, SQTime, RENESTTime, NESTTime, PGDTime, MPETime, SITime,
#     DAARAMConv, NIDAARAMConv, RNIDAARAMConv, SQConv, RENESTConv, NESTConv, PGDConv, MPEConv,
#     file="~/Documents/NIDAAREM/SimulationResults/MovieSims2.RData")








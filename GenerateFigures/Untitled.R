
#library(daaremtest)
library(daarem)
library(nidaarem)

#source("~/ProximalAnderson/Lasso_fns.R")
set.seed(5671)

n <- 1000
p <- 500

ni <- 100000  ## maximum of 200,000 iterations.
ni.prime <- ni
tols <- 1e-12

rho <- 0.9
true.beta <- rt(p, df=3) * rbinom(p, 1, 0.8)
XX <- matrix(0.0, nrow=n, ncol=p)
for(h in 1:n) {
  XX[h,] <- arima.sim(n = p, list(ar = rho), innov=rnorm(p))
}
XX <- scale(XX)
yy <- XX %*% true.beta + rnorm(n)
yy <- yy - mean(yy)
##

lambda.seq <- c(17000, 10000, 6400, 4000, 2300, 1400, 800, 500, 300, 150)

n.lambda <- length(lambda.seq)


## Run 1

lam <- lambda.seq[5]

##
XtX <- crossprod(XX, XX)
VV <- eigen(XtX)
LL <- max(VV$values)
beta.init <- rep(0, p)

ni <- 800
ni.star <- ni*50

gd_vals <- rep(NA, ni+1)
gd_vals[1] <- LassoObjFn(beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL)
beta.new <- beta.init
for(k in 1:ni) {
  beta.new <- GDLassoStep(beta.new, X=XX, y=yy, lambda=lam, stplngth=1/LL) 
  gd_vals[k+1] <- LassoObjFn(beta.new, X=XX, y=yy, lambda=lam, stplngth=1/LL)
}
ans.pgd <- list()
ans.pgd$objfn.track <- gd_vals

ans.nest <- Nesterov(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, 
                     lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols))

ans.renest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=1/LL,
                       restart=TRUE, control=list(maxiter=ni, tol=tols))


ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=FALSE, 
                                                          control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))

ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=FALSE, 
                           control=list(tol=tols, order=5, maxiter=ni, mon.tol=1))


tmp <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=TRUE,
                    control=list(tol=0, order=5, maxiter=ni.star, mon.tol=1))
f.star <- min((-1)*tmp$value.objfn)


pgd.diff1 <- (-1)*ans.pgd$objfn.track - f.star
nest.diff1 <- (-1)*ans.nest$objfn.track - f.star
daaram.diff1 <- (-1)*ans.daaram$objfn.track - f.star
renest.diff1 <- (-1)*ans.renest$objfn.track - f.star

daaram.stop1 <- ans.daaram$fpevals + 1
nest.stop1 <- length(nest.diff1)
pgd.stop1 <- length(pgd.diff1)
renest.stop1 <- length(renest.diff1)

###################################################
### Run 2

lam <- lambda.seq[7]
XtX <- crossprod(XX, XX)
VV <- eigen(XtX)
LL <- max(VV$values)
beta.init <- rep(0, p)

ni <- 800
ni.star <- ni*50

gd_vals <- rep(NA, ni+1)
gd_vals[1] <- LassoObjFn(beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL)
beta.new <- beta.init
for(k in 1:ni) {
  beta.new <- GDLassoStep(beta.new, X=XX, y=yy, lambda=lam, stplngth=1/LL) 
  gd_vals[k+1] <- LassoObjFn(beta.new, X=XX, y=yy, lambda=lam, stplngth=1/LL)
}
ans.pgd <- list()
ans.pgd$objfn.track <- gd_vals

ans.nest <- Nesterov(par=beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, 
                     lambda=lam, stplngth=1/LL, control=list(maxiter=ni, tol=tols))

ans.renest <- Nesterov(beta.init, fixptfn=GDLassoStep, objfn=LassoObjFn, X=XX, y=yy, lambda=lam, stplngth=1/LL,
                       restart=TRUE, control=list(maxiter=ni, tol=tols))


ans.daaram.time <- system.time(ans.daaram <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=FALSE, 
                                                          control=list(tol=tols, order=5, maxiter=ni, mon.tol=c(1,0))))

ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=FALSE, 
                           control=list(tol=tols, order=5, maxiter=ni, mon.tol=1))


tmp <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL,
                    control=list(tol=0, order=5, maxiter=ni.star, mon.tol=1))
f.star <- min((-1)*tmp$value.objfn)


pgd.diff2 <- (-1)*ans.pgd$objfn.track - f.star
nest.diff2 <- (-1)*ans.nest$objfn.track - f.star
daaram.diff2 <- (-1)*ans.daaram$objfn.track - f.star
renest.diff2 <- (-1)*ans.renest$objfn.track - f.star

daaram.stop2 <- ans.daaram$fpevals + 1
nest.stop2 <- length(nest.diff2)
pgd.stop2 <- length(pgd.diff2)
renest.stop2 <- length(renest.diff2)



xl <- c(1,700)
renest.stop1 <- 420

#postscript(file="~/Documents/NIDAAREM/FinalPaperPlots/TracePlots1.eps", width=6, height=4.5, horizontal=FALSE)
par(mar=c(4.1,4.2, .5, .25), mfrow=c(1,1))
plot(daaram.diff1, type="n", xlab="Iteration Number k", las=1, log="y", xlim=xl,
     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), cex.axis=0.8, ylim=c(exp(-24), exp(16)), cex.lab=1.1)
lines(1:pgd.stop1, pgd.diff1, lwd=2)
lines(1:nest.stop1,nest.diff1, col="blue", lwd=2)
lines(1:renest.stop1, renest.diff1[1:renest.stop1], col="green", lwd=2)
lines(1:daaram.stop1, daaram.diff1, col="red", lwd=2)
legend('topright', legend=c("PGD", "Nesterov","Nesterov/restarts", "DAARAM"), col=c("black","blue","green","red"), 
       lwd=3, bty='n', cex=1.1)
#dev.off()

#plot(daaram.diff2, type="n", xlab="Iteration Number k", las=1, log="y", xlim=xl,
#     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), main="(b)")
#lines(1:pgd.stop2, pgd.diff2, lwd=2)
#lines(1:nest.stop2,nest.diff2, col="blue", lwd=2)
#lines(1:renest.stop2, renest.diff2, col="green", lwd=2)
#lines(1:daaram.stop2, daaram.diff2, col="red", lwd=2)
#legend('topright', legend=c("GD", "Nesterov","Nesterov w restarts", "DAARAM"), col=c("black","blue","green","red"), lwd=3, bty='n', cex=1.1)

#dev.off()

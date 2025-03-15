
setwd("~/Documents/nidaarem_reproduce")
library(nidaarem)

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

#stop.ind <- which(diff(ans.nest$objfn.track) < 0)[1]
#ne.tmp <- Nesterov(par=beta.init, fixptfn = GDLassoStep, objfn = LassoObjFn, X=XX, y=yy,
#                   lambda=lam, stplngth=1/LL, control=list(maxiter=stop.ind, tols=1.e-10))
#ne.init <- ne.tmp$par

ans.nidaarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=TRUE,
                             control=list(tol=tols, order=5, maxiter=ni, mon.tol=1))

ans.daarem <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=FALSE, 
                           control=list(tol=tols, order=5, maxiter=ni, mon.tol=1))


tmp <- daarem.lasso(par=beta.init, X=XX, y=yy, lambda=lam, stplngth=1/LL, nesterov.init=TRUE,
                    control=list(tol=0, order=5, maxiter=ni.star, mon.tol=1))
f.star <- min((-1)*tmp$value.objfn)


switch.ind <- which(diff(ans.nest$objfn.track) < 0)[1]
#switch.y <- ans.nest$objfn.track[switch.ind]
pgd.diff1 <- (-1)*ans.pgd$objfn.track - f.star
nest.diff1 <- (-1)*ans.nest$objfn.track - f.star
daarem.diff1 <- (-1)*ans.daarem$objfn.track - f.star
renest.diff1 <- (-1)*ans.renest$objfn.track - f.star
nidaarem.diff1 <- (-1)*ans.nidaarem$objfn.track - f.star
switch.y <- nest.diff1[switch.ind]

daarem.stop1 <- ans.daarem$fpevals + 1
nest.stop1 <- length(nest.diff1)
pgd.stop1 <- length(pgd.diff1)
renest.stop1 <- length(renest.diff1)
nidaarem.stop1 <- length(nidaarem.diff1)

xl <- c(1,700)
renest.stop1 <- 420
graph_cols <- c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
    
postscript(file="Figures/Figure1.eps", width=6, height=4.2, horizontal=FALSE)
par(mar=c(2.7,4.3, .5, .25), mfrow=c(1,2))
plot(daarem.diff1, type="n", xlab="", las=1, log="y", xlim=xl,
     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), cex.axis=0.7, ylim=c(exp(-24), exp(16)), 
     cex.lab=0.8, xaxt='n')
lines(1:pgd.stop1, pgd.diff1, col=graph_cols[1], lwd=2)
lines(1:nest.stop1,nest.diff1, col=graph_cols[2], lwd=2)
lines(1:renest.stop1, renest.diff1[1:renest.stop1], col=graph_cols[3], lwd=2)
lines(1:daarem.stop1, daarem.diff1, col=graph_cols[4], lwd=2)
legend('topright', legend=c("PGD", "Nesterov","Nesterov with restarts", "DAAREM"), col=graph_cols[1:4], 
       lwd=3, bty='n', cex=0.7)
axis(1, at=100*(0:7), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.6)
text(x=350, y=2e-13, labels="Iteration Number k", xpd=NA, cex=0.8)
    
plot(daarem.diff1, type="n", xlab="", las=1, log="y", xlim=xl,
     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), cex.axis=0.7, 
     ylim=c(exp(-24), exp(16)), cex.lab=0.8, xaxt='n')
lines(1:pgd.stop1, pgd.diff1, col=graph_cols[1], lwd=2)
lines(1:nest.stop1,nest.diff1, col=graph_cols[2], lwd=2)
lines(1:renest.stop1, renest.diff1[1:renest.stop1], col=graph_cols[3], lwd=2)
lines(1:daarem.stop1, daarem.diff1, col=graph_cols[4], lwd=2)
lines(1:nidaarem.stop1, nidaarem.diff1, col=graph_cols[5], lwd=2)
legend('topright', legend=c("PGD", "Nesterov","Nesterov with restarts", "DAAREM", "NIDAAREM"), col=graph_cols, 
       lwd=3, bty='n', cex=0.7)
arrows(x0=100, y0=2e4, x1=switch.ind, y1=switch.y, length=0.10, lwd=2)
text(50, 1.5e5, labels="switch", adj=0, cex=0.6)
text(50, 4e4, labels="iteration", adj=0, cex=0.6)
axis(1, at=100*(0:7), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.6)
text(x=350, y=2e-13, labels="Iteration Number k", xpd=NA, cex=0.8)
dev.off()


### Plot panels separately in Figures 1a and 1 b
postscript(file="Figures/Figure1a.eps", width=6, height=4.2, horizontal=FALSE)
par(mar=c(2.7,4.3, .5, .25), mfrow=c(1,1))
plot(daarem.diff1, type="n", xlab="", las=1, log="y", xlim=xl,
     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), cex.axis=0.7, ylim=c(exp(-24), exp(16)), 
     cex.lab=0.8, xaxt='n')
lines(1:pgd.stop1, pgd.diff1, col=graph_cols[1], lwd=2)
lines(1:nest.stop1,nest.diff1, col=graph_cols[2], lwd=2)
lines(1:renest.stop1, renest.diff1[1:renest.stop1], col=graph_cols[3], lwd=2)
lines(1:daarem.stop1, daarem.diff1, col=graph_cols[4], lwd=2)
legend('topright', legend=c("PGD", "Nesterov","Nesterov with restarts", "DAAREM"), col=graph_cols[1:4], 
       lwd=3, bty='n', cex=0.7)
axis(1, at=100*(0:7), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.6)
text(x=350, y=2e-13, labels="Iteration Number k", xpd=NA, cex=0.8)
dev.off()


postscript(file="Figures/Figure1b.eps", width=6, height=4.2, horizontal=FALSE)
par(mar=c(2.7,4.3, .5, .25), mfrow=c(1,1))
plot(daarem.diff1, type="n", xlab="", las=1, log="y", xlim=xl,
     ylab=expression( paste(varphi(x[k]), " - ", varphi^"*")), cex.axis=0.7, 
     ylim=c(exp(-24), exp(16)), cex.lab=0.8, xaxt='n')
lines(1:pgd.stop1, pgd.diff1, col=graph_cols[1], lwd=2)
lines(1:nest.stop1,nest.diff1, col=graph_cols[2], lwd=2)
lines(1:renest.stop1, renest.diff1[1:renest.stop1], col=graph_cols[3], lwd=2)
lines(1:daarem.stop1, daarem.diff1, col=graph_cols[4], lwd=2)
lines(1:nidaarem.stop1, nidaarem.diff1, col=graph_cols[5], lwd=2)
legend('topright', legend=c("PGD", "Nesterov","Nesterov with restarts", "DAAREM", "NIDAAREM"), col=graph_cols, 
       lwd=3, bty='n', cex=0.7)
arrows(x0=100, y0=2e4, x1=switch.ind, y1=switch.y, length=0.10, lwd=2)
text(50, 1.5e5, labels="switch", adj=0, cex=0.6)
text(50, 4e4, labels="iteration", adj=0, cex=0.6)
axis(1, at=100*(0:7), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.6)
text(x=350, y=2e-13, labels="Iteration Number k", xpd=NA, cex=0.8)
dev.off()

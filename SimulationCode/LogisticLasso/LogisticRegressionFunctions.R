NegLogisticObjFn <- function(par, X, Xty, lambda, stplngth) {
  X.beta <- as.vector(X%*%par)
  p1 <- sum(par*Xty)
  
  ## Do log-sum-exp trick for p2
  tmp <- X.beta > 0
  p2.pos <- sum(X.beta[tmp] + log1p(exp(-X.beta[tmp])))
  p2 <- sum(log1p(exp(X.beta[!tmp]))) + p2.pos
  
  p3 <- lambda*sum(abs(par))
  ans <- p3 + p2 - p1
  return(ans)
}


LogisticUpdate <- function(par, X, Xty, lambda, stplngth) {
    phat <- expit(X%*%par)
    xnew <- SoftThresh(par + stplngth*(Xty - crossprod(X, phat)), lambda=lambda*stplngth)
    return(xnew)
}

NegLassoObjFn <- function(bbeta, X, y, lambda, stplngth) {
  X.beta <- as.vector(X%*%bbeta)
  ans <- sum((y - X.beta)*(y - X.beta))/2 + lambda*sum(abs(bbeta))
  return(ans)
}
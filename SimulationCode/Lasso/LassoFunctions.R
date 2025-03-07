GDLassoStep <- function(par, X, y, lambda, stplngth) {
  beta.vec <- par
  X.beta <- X%*%beta.vec
  beta.vec <- SoftThresh(beta.vec + stplngth*crossprod(X, y - X.beta), 
                         lambda=lambda*stplngth)
  return(beta.vec)
}

LassoObjFn <- function(bbeta, X, y, lambda, stplngth) {
  X.beta <- as.vector(X%*%bbeta)
  ans <- sum((y - X.beta)*(y - X.beta))/2 + lambda*sum(abs(bbeta))
  return(-ans)
}

FindStepLengthGaussian <- function(X, y, lam, pmax=10000, nmin=1000) {
  ### Function for finding the steplength for Anderson
  p <- ncol(X)
  n <- nrow(X)
  if(n <= nmin) {
    Lmax <- norm(X, "2")^2
    step <- 1/Lmax
    accurate <- TRUE
  } else if (n > nmin & p <= pmax){
    XtX <- crossprod(X)
    Lmax <- eigs_sym(XtX, 1, retvec=FALSE)$values
    accurate <- TRUE
  } else if (p > pmax & p <= n) {
    XtX <- crossprod(X)
    Lmax <- max(rowSums(abs(XtXt)))
    accurate <- FALSE
  } else if (p > pmax & p > n) {
    XtXt <- tcrossprod(X)
    Lmax <- max(rowSums(abs(XtXt)))
    accurate <- FALSE
  }
  if(!accurate) {
    beta.init <- rep(0, p)
    step <- 1/Lmax
    k <- 1
    lchg <- 1
    while (lchg > 0) {
      k <- 2*k
      lold <- LassoObj(beta.init, X, y, lam, k*step)
      lnew <- LassoObj(GDLassostep(beta.init, X, y, lam, k*step), X, y, lam, k*step)
      lchg <- lnew - lold
    }
    k <- max(2, k/4)
    step <- k/Lmax
  }
  return(step)
}

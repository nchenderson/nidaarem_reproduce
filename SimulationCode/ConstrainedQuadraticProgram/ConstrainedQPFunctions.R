## Functions to perform the proximal gradient descent update, compute
## the objective function, and generate simulated data in the 
## constrained quadratic programming example.
ConstrQPFixpt <- function(par, Q, q, a, b, stplength) {
  tmp <- par - stplength*Q%*%par - stplength*q
  new.par <- tmp
  ind1 <- tmp < a
  ind2 <- tmp > b
  
  new.par[ind1] <- a
  new.par[ind2] <- b
  return(new.par)
}

ConstrQPObj <- function(par, Q, q, a, b, stplength) {
  ind1 <- par < a
  ind2 <- par > b
  num.violations <- sum(ind1) + sum(ind2)
  if(num.violations > 0) {
    ans <- Inf
  } else {
    ans <- sum(par*Q%*%par)/2 + sum(q*par)
  }
  return(ans)
}

ConstrQPNegObj <- function(par, Q, q, a, b, stplength) {
  ind1 <- par < a
  ind2 <- par > b
  num.violations <- sum(ind1) + sum(ind2)
  if(num.violations > 0) {
    ans <- -Inf
  } else {
    ans <- -sum(par*Q%*%par)/2 - sum(q*par)
  }
  return(ans)
}

SimulateConstrQP <- function(n, cond.number) {
   B <- matrix(rnorm(n*n), nrow=n, ncol=n)
   BtB <- crossprod(B)
   U <- tcrossprod(B, chol(solve(BtB)))
   D <- diag(exp(seq(log(cond.number), 0, length.out=n)))
   Q <- U%*%D%*%t(U)
   q <- colMeans(Q) + apply(Q, 2, sd)*rnorm(n)/5
   return(list(Q=Q, q=q))
}

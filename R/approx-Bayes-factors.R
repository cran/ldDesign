calc.B.ABF <- function(Z,n,a=1){
# calculate B from Z
  sqrt(a/(n+a))*exp(n^2*Z^2/(2*(n+a)))
}

calc.Zc.ABF <- function(B,n,a=1){
# calculate critical value of Z from B, n, a
  Zc <- sqrt(2*(n+a)/n^2 * log(sqrt(n+a)*B/sqrt(a)))
  Zc
}

calc.Zalpha.ABF <- function(alpha,n){
    sqrt(1/n*qchisq(1-alpha,1))
}

calc.Balpha.ABF <- function(alpha, n, a){
# find the Bayes factor corresponding to alpha,n,a
  calc.B.ABF(calc.Zalpha.ABF(alpha, n = n), n = n, a=a)
}

calc.alphaB.ABF <- function(B,n,a=1, alpha.start=1e-5, reduction.factor=1.5, niter=20,
	verbose=FALSE, show.progress=FALSE){
# find the alpha value corresp. to B,n,a
# alpha.start: upper bound initial guess
#browser()
alphas <- numeric(niter)
Bs <- numeric(niter)
alphas[1] <- alpha.start
Bs[1] <- calc.B.ABF(calc.Zalpha.ABF(alphas[1],n=n),n=n,a=a)
for(ii in 2:niter){
  if(show.progress)cat(".")
  alphas[ii] <- alphas[ii-1]/reduction.factor
  Bs[ii] <- calc.B.ABF(calc.Zalpha.ABF(alphas[ii],n=n),n=n,a=a)
}
if(show.progress)cat("\n")
if(verbose) print(rbind(alphas,Bs))
alphac <- exp(approx(log(Bs),log(alphas),xout=log(B))$y)

alphac
}


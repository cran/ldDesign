cc.design <- function(B, OR, D,p,q, power, baseline.risk, Dprime=NULL, R=NULL, prevalence=NULL,n.cases,
                      n.controls, model=c("additive","dominant","recessive","general"),a=1,
                      sigma2.eta=NULL, verbose=FALSE, amalgamate.cells=FALSE, 
                      pmin=0.1,pmax=0.99,
                      ninterp=20,print.power.curve=TRUE){
# give n.cases, n.controls as starting values
# print.it: if TRUE, print all calculated values 
# pmin, pmax values for the power curve

#browser()
args <- as.list(sys.call())[-1]
args$power <- NULL; args$print.power.curve <- NULL
args$nmin <- NULL; args$nmax <- NULL; args$ninterp <- NULL;
args$pmin <- NULL; args$pmax <- NULL
power1 <- do.call(cc.power,args)
#args$OR <- NULL
if(model=="general"){
  args$R <- attr(power1,"relative.risk(s)")
}else{
  args$R <- attr(power1,"relative.risk")
}
anti.logit <- function (q){
    u <- exp(q)
    u/(1 + u)
}
logit <- function(p, tol = 1e-06) {
        p1 <- ifelse(p < tol, tol, ifelse(p > 1 - tol, 1 - tol, 
            p))
        log(p1/(1 - p1))
}
power.func1 <- function(x){
#browser()
#  args1 <- eval(args,sys.parent(1))
  args1 <- args
  args1$power <- NULL
  args1$n.cases <- x*n.cases
  args1$n.controls <- x*n.controls
  do.call(cc.power,args1)
}
bracket.power <- function(power){
#browser()
x1 <- 1;
x.factor1 <- 1.5;
x.factor2 <-  1.2
p1 <- power.func1(x1)
delta.logit <- 0.2
plow <- anti.logit(logit(power)-delta.logit)
phigh <- anti.logit(logit(power)+delta.logit)
if(is.nan(p1))return(list(xlow=NaN,power.low=NaN,xhigh=NaN,power.high=NaN))
while(p1 < plow){
  x1 <- x1*x.factor1
  p1 <- power.func1(x1)
}
x2 <- x1
while(p1 > plow){
  x1 <- x1/x.factor2
  p1 <- power.func1(x1)
}
p2 <- power.func1(x2)
while(p2 > phigh){
  x2 <- x2/x.factor1
  p2 <- power.func1(x2)
}
while(p2 < phigh){
  x2 <- x2 * x.factor2
  p2 <- power.func1(x2)
}
list(xlow=x1,power.low=p1,xhigh=x2,power.high=p2)
}
#
#browser()
if(print.power.curve){
  bkt.lower <- bracket.power(pmin)
  xsoln.low <- optimize(function(x){(power.func1(x) - pmin)^2}, lower= bkt.lower$xlow,
                      upper=bkt.lower$xhigh)
  xlow <- xsoln.low$minimum
  bkt.upper <- bracket.power(pmax)
  xsoln.high <- optimize(function(x){(power.func1(x) - pmax)^2}, lower= bkt.upper$xlow,
                      upper=bkt.upper$xhigh)
  xhigh <- xsoln.high$minimum
  xs <- exp(seq(log(xlow),log(xhigh),length=ninterp+1))
  powers <- mapply(power.func1,xs)
  attributes(powers) <- NULL
  n.casess <- floor(xs*n.cases+0.999)
  n.controlss <- floor(xs*n.controls+0.999)
  power.curve <- cbind(n.controls = n.controlss, n.cases = n.casess, power = powers)
  cat("Power curve:\n")
  print(power.curve)
}

  bkt.soln <- bracket.power(power)
  xopt <- optimize(function(x){(power.func1(x) - power)^2}, lower= bkt.soln$xlow,
                      upper=bkt.soln$xhigh)
  xsoln <- xopt$minimum
  n.cases <- floor(xsoln*n.cases+0.999)
  n.controls <- floor(xsoln*n.controls+0.999)
  cat(paste(n.controls,"controls and", n.cases,"cases for power", power, "\n"))
  c(n=n.cases+n.controls, n.controls=n.controls, n.cases=n.cases)
}



cc.power <- function(B, OR=NULL, D, p, q, baseline.risk, Dprime=NULL, R=NULL, prevalence=NULL,n.cases,
                     n.controls, model=c("additive","dominant","recessive","general"),a=1,sigma2.eta=NULL,
                     verbose=FALSE, amalgamate.cells=FALSE, show.attributes=FALSE){
.Call <- sys.call()
args <- as.list(.Call)[-1]
if(missing(sigma2.eta))args$a <- a
#browser()       
#names(args)
args$show.attributes <- NULL
res <- switch(model,
       "additive"=  do.call(ccdesign2x3.ld.power,args),
       "dominant" = {args$model <- NULL;do.call(ccdesign2x2.ld.power.dominant,args)},
       "recessive" = {args$model <- NULL;do.call(ccdesign2x2.ld.power.recessive,args)},
       "general" = {args$model <- NULL;do.call(ccdesign2x3.general.ld.power,args)}  )
if(show.attributes){
  res
}else{
  res0 <- res
  attributes(res0) <- NULL
  #print(res0)
  #invisible(res)
  res0
}
}


ccdesign2x3.ld.power <- function(B,OR,D,p,q, baseline.risk, Dprime=NULL, R=NULL,
                                 prevalence=NULL, n.cases, n.controls, model="additive",
                                 a=1,sigma2.eta=NULL, verbose=FALSE,
                                 amalgamate.cells=FALSE, optim.trace=FALSE){
# p: marker allele frequency
# q: trait locus allele frequency
# OR: odds ratio  
#browser()
  Dmax <- min((1-p)*q,(1-q)*p)
  Dmin <- max(-p*q,-(1-p)*(1-q))

if(missing(D) && !missing(Dprime)){
  if(Dprime > 0) D <- Dprime*Dmax else D <-  (-Dprime)*Dmin
}  
#
if(missing(D) && missing(Dprime))Dprime <- 1
if(missing(D)){
  D <- Dmax
}else{
  if(D > Dmax || D < Dmin)
    stop(paste("must have  max(-p*q,-(1-p)*(1-q)) = Dmin =", Dmin,
               "< D < min((1-p)*q,(1-q)*p) = Dmax =", Dmax))
}
#
if(missing(baseline.risk) && missing(prevalence))stop("must give baseline.risk or prevalence")
if(missing(OR) && missing(R))stop("must give odds ratio or relative risk (R)")
#
if(!missing(baseline.risk))q0 <- baseline.risk else q0 <- NULL
if(!missing(prevalence))nu <- prevalence else nu <- NULL
#
#browser()
if(!missing(R)){
  if(missing(prevalence)){
    if(missing(baseline.risk))stop("must give baseline.risk or prevalence")
    nu <- prevalence <- q0*(q^2*R^2 + (2*q - 2*q^2)*R + q^2 - 2*q + 1)
  }else{
    q0 <- nu/(q^2*R^2 + (2*q - 2*q^2)*R + q^2 - 2*q + 1)
  }
  OR <- (q0-1)*R/(q0*R-1)
}
#
eta <- log(OR)
#
# [R = OR/(q0*OR-q0+1)]$
if(missing(R)){
if(missing(baseline.risk)){
  if(missing(prevalence))stop("need to give baseline.risk or prevalence")
  nu <- prevalence
#  browser()
  R0 <- 1.25
  q0 <- nu/(q^2*R0^2 + (2*q - 2*q^2)*R0 + q^2 - 2*q + 1)
  opt.fn <- function(x){
    q0 <- exp(x[1]); R <- exp(x[2])
    e1 <- OR - ((q0-1)*R)/(q0*R-1)
    e2 <- 1 - ( q^2*q0*R^2+2*(1-q)*q*q0*R+(1-q)^2*q0)/nu
    err <- e1^2+e2^2
    err
  }
  q0.soln <- NA
  for(ii in 1:10){
    opt1 <- optim(log(c(q0,R0)), opt.fn, method="Nelder-Mead",control=list(trace=optim.trace))
    opt1
    q0 <- q0.soln <- exp(opt1$par[1])
    R <- R0 <- exp(opt1$par[2])
    if(opt1$value < 1e-10)break;
  }
  if(opt1$convergence !=0){
    stop("failed to find solution for q0, R")
  }else{
    if(verbose)
      cat("found solution: q0 =", q0.soln,"R =", R,", squared error =", opt1$value,"\n")
  }
# problem theta used befor defined.
#  q0 <- baseline.risk <- nu/(q^2*exp(2*theta) + (2*q - 2*q^2)*exp(theta) + q^2 - 2*q + 1)
  q0 <- q0.soln 
  R <- OR/(q0*OR-q0+1)
}else{
  if(missing(prevalence) && missing(R)){
    q0 <- baseline.risk
    R <- OR/(q0*OR-q0+1)
    theta <- log(R)
    nu <-  prevalence <- q^2*q0*R^2 + 2*(1-q)*q*q0*R + (1-q)^2*q0
  }
}
}
#
#else{ # missing(OR)
#
#cat("handled above, should not get here\n")
#  if(missing(baseline.risk)){
#    nu <- prevalence
#    q0 <- baseline.risk <- nu/(q^2*R^2 + 2*(1-q)*q*R + (1-q)^2)
#  }
#  if(missing(prevalence)){
#    q0 <- baseline.risk
#    nu <-  prevalence <- q^2*q0*R^2 + 2*(1-q)*q*q0*R + (1-q)^2*q0
#  }
#  OR <- (q0-1)*R/(q0*R-1)
#}
if(missing(R)){
  R <- OR/(q0*OR-q0+1)
  theta <- log(R)
}
#browser()
#
q11 <-  (1-q)^2*(1-q0)/(1-nu)
q12 <- 2*(1-q)*q*(1-q0*R)/(1-nu)
q13 <- q^2*(1-q0*R^2)/(1-nu)
q21 <- (1-q)^2*q0/nu
q22 <- 2*(1-q)*q*q0*R/nu
q23 <- q^2*q0*R^2/nu
#
if(missing(D)){# non-LD case
p11 <- q11; p12 <- q12; p13 <- q13; p21 <- q21; p22 <- q22; p23 <- q23
}else{# LD case
p11 <- q11*(D+(1-p)*(1-q))^2/(1-q)^2  +  q12*((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)  +
  q13*((1-p)*q-D)^2/q^2;
p13 <- q11*(p*(1-q)-D)^2/(1-q)^2  +  q12*(p*(1-q)-D)*(D+p*q)/((1-q)*q)+  q13*(D+p*q)^2/q^2
p21 <- q21*(D+(1-p)*(1-q))^2/(1-q)^2  +  q22*((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)  +
  q23*((1-p)*q-D)^2/q^2;
p23 <- q21*(p*(1-q)-D)^2/(1-q)^2 +  q22*(p*(1-q)-D)*(D+p*q)/((1-q)*q) + q23*(D+p*q)^2/q^2
# NEW: add equations for p12, p22
p12 <- q12*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))/((1-q)*q) +
    2*q13*((1-p)*q-D)*(D+p*q)/q^2+2*q11*(p*(1-q)-D)*(D+(1-p)*(1-q))/(1-q)^2
p22 <- q22*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))/((1-q)*q) +
    2*q23*((1-p)*q-D)*(D+p*q)/q^2+2*q21*(p*(1-q)-D)*(D+(1-p)*(1-q))/(1-q)^2
#browser()
}
#browser()
#
n1 <- n.controls
n2 <- n.cases
n <- n1+n2
etax <- 1/2*(log(p11) - log(p13) - log(p21) + log(p23))
sigma1 <- 1/2*sqrt(n*(1/(n1*p11) + 1/(n1*p13) + 1/(n2*p21) + 1/(n2*p23)))

#browser()

calc.V1 <- function(n,c1,c2,p11,p12,p13){
  V1 <- (c1^2*(p12*p13+p11*p13)+(1-c1)^2*(p11*p13+p11*p12) - 
         2*(1-c1)*c1*p11*p13)/(n1*p11*p12*p13)
  return(V1);
}
calc.V2 <- function(n,c1,c2,p21,p22,p23){
  V2 <- (c1^2*(p22*p23+p21*p23)+(1-c1)^2*(p21*p23+p21*p22) -
         2*(1-c1)*c1*p21*p23)/(n1*p21*p22*p23)
  return(V2);
}

if(F){
c1.optimum <- (((2*n2*p11*p13+n2*p11*p12)*p21*p22+2*n1*p11*p12*p13*p21)*p23 + n1*p11*p12*p13*p21*p22) /
      (((((n2*p12+4*n2*p11)*p13+n2*p11*p12)*p21+n1*p11*p12*p13)*p22 +4*n1*p11*p12*p13*p21)*p23
          +n1*p11*p12*p13*p21*p22)
c1 <- c1.optimum
c2 <- 1-c1
etax.opt <- c1*(log(p11) - log(p12) - log(p21) + log(p22)) +
            (1-c1)*(log(p12) - log(p13) - log(p22) + log(p23))

V1 <- calc.V1(n1,c1,c2,p11,p12,p13)
V2 <- calc.V2(n2,c1,c2,p21,p22,p23)

sigma1.opt <- sqrt(V1+V2)*sqrt(n)

}
#

#
#
calc.ncp <- function(p11,p12,p13,p21,p22,p23,n1,n2,c1){
ncp <- ((1-c1)*(log(p23)-log(p22)-log(p13)+log(p12))
        +  c1*(log(p22)-log(p21)-log(p12)+log(p11)))^2 /
    ((c1^2*(p22*p23+p21*p23)+(1-c1)^2*(p21*p23+p21*p22)-2*(1-c1)*c1*p21*p23)
     /(n2*p21*p22*p23)
      +  (c1^2*(p12*p13+p11*p13)+(1-c1)^2*(p11*p13+p11*p12)-2*(1-c1)*c1*p11*p13)
      /(n1*p11*p12*p13))
return(ncp)
}
#
#
c1.opt.fn <- function(c1){
  calc.ncp(p11,p12,p13,p21,p22,p23,n1[1],n2[1],c1)
}
#
c1.opt.fn(c1=0.5)
ncp.values <- c1.opt.fn(c1 = seq(0,1,by=0.01))
if(verbose){
  plot(seq(0,1,by=0.01),ncp.values,xlab="c1 ", ylab="NCP")
}
#
c1.opt <- optimize(c1.opt.fn, interval=c(0,1), maximum=TRUE)
c1.optimum <- c1.opt$maximum
c1 <- c1.optimum
c2 <- 1 - c1
etax.opt <- c1*(log(p11) - log(p12) - log(p21) + log(p22)) +
            (1-c1)*(log(p12) - log(p13) - log(p22) + log(p23))
V1 <- calc.V1(n1,c1,c2,p11,p12,p13)
V2 <- calc.V2(n2,c1,c2,p21,p22,p23)
sigma1.opt <- sqrt(V1+V2)*sqrt(n)

#browser()


# calculate `a' (equivalent sample size for prior) if missing
if(missing(D)){
    detax.deta <- 1
}else{
    detax.deta <- cc2x3.calc.detax.deta(D,p,q,p11,p13,p21,p23,q0,R,verbose=verbose)
}
if(missing(D)){
    detax.deta.opt <- 1
}else{
#browser()
    detax.deta.opt <- cc2x3.calc.detax.deta.opt(D,p,q,p11,p12,p13,p21,p22,p23,q0,R,nu,
                                                c1, verbose=verbose)
}
#verbose <- TRUE
if(verbose){
cat(paste("D =", D,", Dmin =",Dmin, "Dmax =", Dmax,"\n"))
cat(paste("q0 =",q0,", R =",R, ", nu =", nu,", OR =", OR, "\n"))
cat(paste("q11 =", q11, ", q12 =", q12, ", q13 =", q13, ", q21 =", q21, ", q22 =", q22,
          ", q23 =", q23, "\n"))
cat(paste("p11 =", p11, ", p12 =",p12,", p13 =", p13, ", p21 =", p21, ", p22 =",p22,
          ", p23 =", p23,"\n"))
}

if(missing(a)){
  if(missing(sigma2.eta))stop("must give \"a\" or \"sigma2.eta\"")
  sigma2.etax <- detax.deta.opt^2*sigma2.eta
  ax <- a <- sigma1^2/sigma2.etax
}else{
  ax <- a/detax.deta.opt^2  
}

Zc <- calc.Zc.ABF(B,n=n.cases+n.controls,a=ax)
n <- n.cases+n.controls

# corresponding alpha value
# alphac <- 2*(1-pnorm(sqrt(n)*Zc))
alphac <- 1 - pchisq(n*Zc^2,1)
alphac

# equivalent non-centrality parameter in chi-squared test
ncp <- n*(etax/sigma1)^2
ncp

ncp.opt <- n*(etax.opt/sigma1.opt)^2
ncp.opt

if(verbose)
  cat(paste("after setting ncp.opt, ncp =", ncp, "ncp.opt =", ncp.opt,"\n"))

# use these values when allele frequency is low (amalgamate cells with 2 risk alleles,
#  cf dominant model)
sigma1.2x2 <- sqrt(n*(1/(n1*p11) + 1/(n1*(p12+p13)) + 1/(n2*p21) + 1/(n2*(p22+p23))))
etax.2x2 <- (log(p11) - log(p12+p13) - log(p21) + log(p22+p23))
detax.deta2x2 <- cc2x2.calc.detax.deta.dominant(D,p,q,p11,p12+p13,p21,p22+p23,q0,R,nu,
                                                verbose=verbose)
ax2 <- a/detax.deta2x2^2
Zc2x2 <- calc.Zc.ABF(B,n=n.cases+n.controls,a=ax2)
power.2x2 <- 1 - pnorm(sqrt(n)*(Zc2x2 - etax.2x2/sigma1.2x2))

## strictly speaking
if(F){
  power <- 1 - (pnorm(sqrt(n)*(Zc - etax/sigma1)) - pnorm(sqrt(n)*(-Zc - eta/sigma1)))
}
## neglecting contribution where statistic and value of opposite signs

power.non.opt <- 1 - pnorm(sqrt(n)*(Zc - etax/sigma1))
power <- 1 - pnorm(sqrt(n)*(Zc - etax.opt/sigma1.opt))
    if(verbose){
      cat(paste("power (not optimised) = ", power.non.opt," power (optimised, c1 =",c1,") =",power,"\n"))
      cat("pij:\n")
      # p12 <- 1-p11-p13; p22 <- 1-p21-p23;
      print(matrix(c(p11,p12,p13,p21,p22,p23),nr=2,byr=T))
      cat("qij:\n")
      q12 <- 1-q11-q13; q22 <- 1-q21-q23;
      print(matrix(c(q11,q12,q13,q21,q22,q23),nr=2,byr=T))
      print(c(n=n,q0=q0,R=R,nu=nu,ax=ax,etax=etax,sigma1=sigma1,
              etaxoversigma1=etax.opt/sigma1.opt,
              etax.2x2=etax.2x2, sigma1.2x2=sigma1.2x2,
              etax2x2oversigma12x2=etax.2x2/sigma1.2x2,
              OR=OR,B=B,Zc=Zc,power=power))
    }

if(etax.2x2/sigma1.2x2 > 1.05*etax.opt/sigma1.opt && ! amalgamate.cells){
  cat("low allele frequencies: power greater with amalgamated cells\n")
  cat("power (2x3) =", power, ", power (2x2) =", power.2x2,"\n")
}
if(amalgamate.cells){ 
 attr(power.2x2,"model") <- "additive (2x2)"
  attr(power.2x2,"n.cases") <- n.cases
  attr(power.2x2,"n.controls") <- n.controls
  attr(power.2x2,"baseline.risk") <- q0
  attr(power.2x2,"relative.risk") <- R
  attr(power.2x2,"odds-ratio") <- OR

  attr(power.2x2,"Dprime") <- Dprime
  attr(power.2x2,"p") <- p
  attr(power.2x2,"q") <- q
  attr(power.2x2,"ncp") <- ncp
  attr(power.2x2,"B") <- B
  attr(power.2x2,"alphac") <- alphac

  return(power.2x2)
}

attr(power,"model") <- "additive"
attr(power,"n.cases") <- n.cases
attr(power,"n.controls") <- n.controls
attr(power,"prevalence") <- nu
attr(power,"baseline.risk") <- q0
attr(power,"relative.risk") <- R
attr(power,"odds-ratio") <- OR
attr(power,"Dprime") <- Dprime
attr(power,"p") <- p
attr(power,"q") <- q
attr(power,"ncp") <- ncp.opt
attr(power,"B") <- B
attr(power,"alphac") <- alphac
attr(power,"c1.opt") <- c1.optimum
ps <- rbind(c(p11,p12,p13),c(p21,p22,p23))
dimnames(ps) <- list(c("control","case"),c("bb","Bb","BB"))
attr(power,"ps") <- ps

power
}


ccdesign2x3.general.ld.power <- function(B,OR,D,p,q, baseline.risk, Dprime=NULL, R=NULL,
                                         prevalence=NULL, n.cases, n.controls, model="general",
                                         a=1,sigma2.eta=NULL, verbose=FALSE,
                                 amalgamate.cells=FALSE, optim.trace=FALSE){
# p: marker allele frequency
# q: trait locus allele frequency
# OR: odds ratio  
# Note: does not vectorise!

  Dmax <- min((1-p)*q,(1-q)*p)
  Dmin <- max(-p*q,-(1-p)*(1-q))
#
if(missing(D) && !missing(Dprime)){
  if(Dprime > 0) D <- Dprime*Dmax else D <-  (-Dprime)*Dmin
}  
#
if(missing(D) && missing(Dprime))Dprime <- 1
if(missing(D)){
  D <- Dmax
}else{
  if(D > Dmax || D < Dmin)
    stop(paste("must have  max(-p*q,-(1-p)*(1-q)) = Dmin =", Dmin,
               "< D < min((1-p)*q,(1-q)*p) = Dmax =", Dmax))
}
#
if(missing(baseline.risk) && missing(prevalence))stop("must give baseline.risk or prevalence")
if(missing(OR) && missing(R))stop("must give odds ratio or relative risk (R)")
#
if(!missing(OR) && length(OR)!=2)stop("OR should be a vector of length 2")
if(!missing(R) && length(R)!=2)stop("R should be a vector of length 2")
#
if(!missing(baseline.risk))q0 <- baseline.risk else q0 <- NULL
if(!missing(prevalence))nu <- prevalence else nu <- NULL
#
if(missing(prevalence) && missing(baseline.risk))
  stop("need to give baseline.risk or prevalence")
#
#browser()
if(missing(OR) && !missing(R)){
  OR <- numeric(2)
  if(!missing(prevalence)) q0 <- nu/(q^2*R[1]*R[2] + (2*q - 2*q^2)*R[1] + q^2 - 2*q + 1)
  if(!missing(baseline.risk))nu <- q0*(q^2*R[1]*R[2] + (2*q - 2*q^2)*R[1] + q^2 - 2*q + 1)
  OR[1] <- (q0-1)*R[1]/(q0*R[1]-1)
  q1 <- q0*R[1]
  OR[2] <- (q1-1)*R[2]/(q1*R[2]-1)
}
#
eta <- log(OR)
#
# [R = OR/(q0*OR-q0+1)]$
if(missing(R)){
if(missing(baseline.risk)){
  if(missing(prevalence))stop("need to give baseline.risk or prevalence")
  nu <- prevalence
#  browser()
  R0 <- 1.25
  q0 <- nu/(q^2*R0^2 + (2*q - 2*q^2)*R0 + q^2 - 2*q + 1)
  opt.fn <- function(x){
    q0 <- exp(x[1]); # R1 <- exp(x[2]); R2 <- exp(x[3])
    OR1 <- OR[1]
    OR2 <- OR[2]
    nuhat <-(  ((((q^2-2*q+1)*q0^3+(2*q-q^2)*q0^2)*OR1^2 +
                 ((-q^2+2*q-1)*q0^3+(1-2*q)*q0^2+q^2*q0)*OR1)*OR2  +
                ((-q^2+2*q-1)*q0^3+(3*q^2-4*q+1)*q0^2+(2*q-2*q^2)*q0)*OR1+(q^2-2*q+1)*q0^3 +
                (-2*q^2+4*q-2)*q0^2+(q^2-2*q+1)*q0) /
             ((q0^2*OR1^2+(q0-q0^2)*OR1)*OR2+(q0-q0^2)*OR1+q0^2-2*q0+1) )
    e1 <- (nu-nuhat)/nu
    e1^2
  }
  # TODO try guesses to narrow interval
#  
  opt1 <- optimize(opt.fn, lower=log(1e-4), upper=log(1), tol=1e-6)
  q0 <- exp(opt1$minimum)
  R[1] <- OR[1]/(q0*OR[1]-q0+1)
  q1 <- q0*R[1]
  R[2] <- OR[2]/(q1*OR[2]-q1+1)
  if(verbose)
      cat("found solution: q0 =", q0,"R1 =", R[1],"R2 =", R[2],", squared error =",
          opt1$objective,"\n")
  if(F){
    if(opt1$convergence !=0){
      stop("failed to find solution for q0, R")
    }else{
      if(verbose)
        cat("found solution: q0 =", q0.soln,"R =", R,", squared error =", opt1$value,"\n")
    }
  }
}else{
  if(missing(prevalence) && missing(R)){
     q0 <- baseline.risk
     R[1] <- OR[1]/(q0*OR[1]-q0+1)
     q1 <- q0*R[1]
     R[2] <- OR[2]/(q1*OR[2]-q1+1)
     nu <-  prevalence <- q^2*q0*R[1]*R[2] + 2*(1-q)*q*q0*R[1] + (1-q)^2*q0
  }
}
}

if(!missing(baseline.risk)){
  q0 <- baseline.risk
  nu <-  prevalence <- q^2*q0*R[1]*R[2] + 2*(1-q)*q*q0*R[1] + (1-q)^2*q0
}
if(!missing(prevalence)){
  nu <- prevalence
  q0 <- baseline.risk <- nu/(q^2*R[1]*R[2] + 2*(1-q)*q*R[1] + (1-q)^2)
}
#R <- OR/(q0*OR-q0+1)
#theta <- log(R)
q11 <-  (1-q)^2*(1-q0)/(1-nu)
q12 <- 2*(1-q)*q*(1-q0*R[1])/(1-nu)
q13 <- q^2*(1-q0*R[1]*R[2])/(1-nu)
q21 <- (1-q)^2*q0/nu
q22 <- 2*(1-q)*q*q0*R[1]/nu
q23 <- q^2*q0*R[1]*R[2]/nu
#

if(missing(D)){# non-LD case
p11 <- q11; p12 <- q12; p13 <- q13; p21 <- q21; p22 <- q22; p23 <- q23
}else{# LD case
p11 <- q11*(D+(1-p)*(1-q))^2/(1-q)^2  +  q12*((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)  +
  q13*((1-p)*q-D)^2/q^2;
p13 <- q11*(p*(1-q)-D)^2/(1-q)^2  +  q12*(p*(1-q)-D)*(D+p*q)/((1-q)*q)+  q13*(D+p*q)^2/q^2
p21 <- q21*(D+(1-p)*(1-q))^2/(1-q)^2  +  q22*((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)  +
  q23*((1-p)*q-D)^2/q^2;
p23 <- q21*(p*(1-q)-D)^2/(1-q)^2 +  q22*(p*(1-q)-D)*(D+p*q)/((1-q)*q) + q23*(D+p*q)^2/q^2
# NEW: add equations for p12, p22
p12 <- q12*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))/((1-q)*q) +
    2*q13*((1-p)*q-D)*(D+p*q)/q^2+2*q11*(p*(1-q)-D)*(D+(1-p)*(1-q))/(1-q)^2
p22 <- q22*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))/((1-q)*q) +
    2*q23*((1-p)*q-D)*(D+p*q)/q^2+2*q21*(p*(1-q)-D)*(D+(1-p)*(1-q))/(1-q)^2
#browser()
}
#browser()
n1 <- n.controls
n2 <- n.cases
n <- n1+n2
#
etax1 <- log(p11) - log(p12) - log(p21) + log(p22)
etax2 <- log(p12) - log(p13) - log(p22) + log(p23)
etax <- c(etax1,etax2)
#
calc.V <- function(n1,n2){
  1/n1*matrix(c( (p11+p12)/(p11*p12), -1/p12, -1/p12, (p12+p13)/(p12*p13)), nr=2,nc=2) +
  1/n2*matrix(c( (p21+p22)/(p21*p22), -1/p22, -1/p22, (p22+p23)/(p22*p23)), nr=2,nc=2)
}
#
if(F){
  if(length(n.cases)>1){
  cat("n.cases = ", n.cases,"\n")
  cat("ccdesign2x3.general.ld.power, before calc.V\n")
  browser()
  }  
}
#
Vn1n2 <- calc.V(n1,n2)
c1 <- n1/n
# V1 <- calc.V(c1,1-c1)
V1 <- Vn1n2*n
# square  root geometric mean variance 
sigma1 <- sqrt(sqrt(det(V1)))
# need correct ax here
# calculate `a' (equivalent sample size for prior) if missing
if(missing(D)){
    detax.deta <- diag(rep(1,2))
}else{
    detax.deta <- cc2x3.calc.detax.deta.general(D,p,q,p11,p12,p13,p21,p22,p23,q0,R,nu,
                                                verbose=verbose)
}
detax.deta.scale <- sqrt(abs(det(detax.deta)))
if(missing(a)){
  if(missing(sigma2.eta))stop("must give \"a\" or \"sigma2.eta\"")
  sigma2.etax <- detax.deta.scale^2*sigma2.eta
  ax <- sigma1^2/sigma2.etax
}else{
  ax <- a/(detax.deta.scale^2)
}
#
H <- (1+ax/n)*solve(Vn1n2)
# ncp <- ax^2/(n*(n+a)) * t(etax) %*% H %*% etax
ncp <- n/(n+ax) * t(etax) %*% H %*% etax
#
if(F){
# compare terms inside the log for the power 
Ha <- ax/n*solve(Vn1n2)
sqrt(det(H)/det(Ha))
(n+ax)/ax
X2c <- 2*(n+ax)/n * log(sqrt(det(H)/det(Ha)) * B)
}
# a-->ax!
# X2c <- 2*(n+ax)/n * log((n+ax)/a * B)
X2c <- 2*(n+ax)/n * log((n+ax)/ax * B)
alphac <- 1 - pchisq(X2c,df=2)
power <- 1 - pchisq(X2c,df=2,ncp=ncp)
#
attr(power,"model") <- "general"
attr(power,"n.cases") <- n.cases
attr(power,"n.controls") <- n.controls
attr(power,"prevalence") <- nu
attr(power,"baseline.risk") <- q0
attr(power,"relative.risk(s)") <- R
attr(power,"odds-ratio(s)") <- OR
attr(power,"Dprime") <- Dprime
attr(power,"p") <- p
attr(power,"q") <- q
attr(power,"ncp") <- ncp
attr(power,"B") <- B
attr(power,"alphac") <- alphac
ps <- rbind(c(p11,p12,p13),c(p21,p22,p23))
dimnames(ps) <- list(c("control","case"),c("bb","Bb","BB"))
attr(power,"ps") <- ps
#attr(power,"H") <- H

power
}

cc2x3.calc.detax.deta.general <- function(D,p,q,p11,p12,p13,p21,p22,p23,q0,R,nu,
                                                verbose=FALSE){
R1 <- R[1]
R2 <- R[2]
# nu <-  q0*(q^2*R[1]*R[2]+2*(1-q)*q*R[1]+(1-q)^2)
# deta1x.deta1
#
#['diff(q0,R1,1) = -nu*(q^2*R2-2*q^2+2*q)
#                   /(q^2*R1*R2+(2*q-2*q^2)*R1+q^2-2*q+1)^2]$
#['diff(q0,R2,1) = -q^2*q0*R1/(q^2*R1*R2+(2*q-2*q^2)*R1+q^2-2*q+1)]$
#
if(F){ #check equal
dq0.dR1 <- -nu*(q^2*R2-2*q^2+2*q)/(q^2*R1*R2+(2*q-2*q^2)*R1+q^2-2*q+1)^2
dq0.dR2 <- -q^2*q0*R1/(q^2*R1*R2+(2*q-2*q^2)*R1+q^2-2*q+1)
}
dq0.dR1 <- -q0^2/nu *(q^2*R2 + 2*q*(1-q))
dq0.dR2 <- -q0^2/nu* q^2*R1
#['diff(OR1,R1,1)/OR1 = ('diff(q0,R1,1)*R1^2-'diff(q0,R1,1)*R1-q0+1)
#                     /((q0^2-q0)*R1^2+(1-q0)*R1)]$
# ['diff(OR1,R2,1)/OR1 = ('diff(q0,R2,1)*R1-'diff(q0,R2,1))/((q0^2-q0)*R1-q0+1)]$
deta1.dR1 <- (dq0.dR1*(R1^2-R1)-q0+1)/((q0^2-q0)*R1^2+(1-q0)*R1)
deta1.dR2 <- (dq0.dR2*(R1-1))/(-(1-q0)*q0*R1 + 1-q0)
#['diff(OR2,R1,1)/OR2 = (('diff(q0,R1,1)*R1+q0)*R2-'diff(q0,R1,1)*R1-q0)
#                     /((q0^2*R1^2-q0*R1)*R2-q0*R1+1)]$
#['diff(OR2,R2,1)/OR2 = ('diff(q0,R2,1)*R1*R2^2-'diff(q0,R2,1)*R1*R2-q0*R1+1)
#                     /((q0^2*R1^2-q0*R1)*R2^2+(1-q0*R1)*R2)]$
deta2.dR1 <- ((dq0.dR1*R1+q0)*R2-dq0.dR1*R1-q0)/((q0^2*R1^2-q0*R1)*R2-q0*R1+1)
deta2.dR2 <- (dq0.dR2*R1*R2^2-dq0.dR2*R1*R2-q0*R1+1)/((q0^2*R1^2-q0*R1)*R2^2+(1-q0*R1)*R2)
# 'diff(p[11],q11,1) = (D+(1-p)*(1-q))^2/(1-q)^2$
dp11dq11 <- (D+(1-p)*(1-q))^2/(1-q)^2
# 'diff(p[11],q12,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)$
dp11dq12 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
# 'diff(p[11],q13,1) = ((1-p)*q-D)^2/q^2$
dp11dq13 <-  ((1-p)*q-D)^2/q^2
dp11dqs <- c(dp11dq11,dp11dq12,dp11dq13)
#['diff(p12,q11,1) = (D+(1-p)*(1-q))^2/(1-q)^2,
# 'diff(p12,q12,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q),
# 'diff(p12,q13,1) = ((1-p)*q-D)^2/q^2]$
dp12dq11 <- (D+(1-p)*(1-q))^2/(1-q)^2
dp12dq12 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
dp12dq13 <- ((1-p)*q-D)^2/q^2
dp12dqs <- c(dp12dq11,dp12dq12,dp12dq13)
# 'diff(p[13],q11,1) = (D+p*q)^2/(1-q)^2$
#!! dp13dq11 <- (D+p*q)^2/(1-q)^2
dp13dq11 <- (p*(1-q)-D)^2/(1-q)^2
# 'diff(p[13],q12,1) = (p*(1-q)-D)*(D+p*q)/((1-q)*q)$
dp13dq12 <- (p*(1-q)-D)*(D+p*q)/((1-q)*q)
# 'diff(p[13],q13,1) = (D+p*q)^2/q^2$
dp13dq13 <- (D+p*q)^2/q^2
dp13dqs <- c(dp13dq11,dp13dq12,dp13dq13)
# 'diff(p[21],q21,1) = (D+(1-p)*(1-q))^2/(1-q)^2$
dp21dq21 <-  (D+(1-p)*(1-q))^2/(1-q)^2
# 'diff(p[21],q22,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)$
dp21dq22 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
# 'diff(p[21],q23,1) = ((1-p)*q-D)^2/q^2$
dp21dq23 <- ((1-p)*q-D)^2/q^2
dp21dqs <- c(dp21dq21,dp21dq22,dp21dq23)
# ['diff(p22,q21,1) = (D+(1-p)*(1-q))^2/(1-q)^2,
#  'diff(p22,q22,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q),
#  'diff(p22,q23,1) = ((1-p)*q-D)^2/q^2]$
dp22dq21 <- (D+(1-p)*(1-q))^2/(1-q)^2
dp22dq22 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
dp22dq23 <- ((1-p)*q-D)^2/q^2
dp22dqs <- c(dp22dq21,dp22dq22,dp22dq23)
#! 'diff(p[23],q21,1) = (D+p*q)^2/(1-q)^2$ 
#! dp23dq21 <- (D+p*q)^2/(1-q)^2
dp23dq21 <- (p*(1-q)-D)^2/(1-q)^2
# 'diff(p[23],q22,1) = (p*(1-q)-D)*(D+p*q)/((1-q)*q)$
dp23dq22 <- (p*(1-q)-D)*(D+p*q)/((1-q)*q)
# 'diff(p[23],q23,1) = (D+p*q)^2/q^2$
dp23dq23 <- (D+p*q)^2/q^2
dp23dqs <- c(dp23dq21,dp23dq22,dp23dq23)
#['diff(q11,R1,1) = -(1-q)^2*'diff(q0,R1,1)/(1-nu),
# 'diff(q12,R1,1) = 2*(1-q)*q*(-'diff(q0,R1,1)*R1-q0)/(1-nu),
# 'diff(q13,R1,1) = q^2*(-'diff(q0,R1,1)*R1*R2-q0*R2)/(1-nu)]$
#['diff(q11,R2,1) = -(1-q)^2*'diff(q0,R2,1)/(1-nu),
# 'diff(q12,R2,1) = -2*(1-q)*q*'diff(q0,R2,1)*R1/(1-nu),
# 'diff(q13,R2,1) = q^2*(-'diff(q0,R2,1)*R1*R2-q0*R1)/(1-nu)]$
dq11.dR1 <- -(1-q)^2/(1-nu)*dq0.dR1
dq11.dR2 <- -(1-q)^2/(1-nu)*dq0.dR2
dq12.dR1 <- -2*q*(1-q)/(1-nu)*(q0 + R1*dq0.dR1)
dq12.dR2 <- -2*q*(1-q)/(1-nu)*R1*dq0.dR2
dq13.dR1 <- -q^2/(1-nu)*(dq0.dR1*R1*R2+q0*R2)
dq13.dR2 <- -q^2/(1-nu)*(dq0.dR2*R1*R2 + q0*R1)
#
dq1s.dR1 <- c(dq11.dR1,dq12.dR1,dq13.dR1)
dq1s.dR2 <- c(dq11.dR2,dq12.dR2,dq13.dR2)
#['diff(q21,R1,1) = (1-q)^2*'diff(q0,R1,1)/nu,
# 'diff(q22,R1,1) = 2*(1-q)*q*'diff(q0,R1,1)*R1/nu+2*(1-q)*q*q0/nu,
# 'diff(q23,R1,1) = q^2*'diff(q0,R1,1)*R1*R2/nu+q^2*q0*R2/nu]$
#['diff(q21,R2,1) = (1-q)^2*'diff(q0,R2,1)/nu,
# 'diff(q22,R2,1) = 2*(1-q)*q*'diff(q0,R2,1)*R1/nu,
# 'diff(q23,R2,1) = q^2*'diff(q0,R2,1)*R1*R2/nu+q^2*q0*R1/nu]$
dq21.dR1 <- (1-q)^2/nu * dq0.dR1
dq21.dR2 <- (1-q)^2/nu * dq0.dR2
dq22.dR1 <- 2*q*(1-q)/nu * (q0+R1*dq0.dR1)
dq22.dR2 <- 2*q*(1-q)/nu * R1*dq0.dR2
dq23.dR1 <- q^2/nu * (dq0.dR1*R1*R2 + q0*R2)
dq23.dR2 <- q^2/nu * (dq0.dR2*R1*R2 + q0*R1)
#
dq2s.dR1 <- c(dq21.dR1,dq22.dR1,dq23.dR1)
dq2s.dR2 <- c(dq21.dR2,dq22.dR2,dq23.dR2)
#
detax1.dR1 <- (1/p11*dotprod(dp11dqs,dq1s.dR1) - 1/p12*dotprod(dp12dqs,dq1s.dR1) -
                 1/p21*dotprod(dp21dqs,dq2s.dR1)   + 1/p22*dotprod(dp22dqs,dq2s.dR1))
detax1.dR2 <- (1/p11*dotprod(dp11dqs,dq1s.dR2) - 1/p12*dotprod(dp12dqs,dq1s.dR2) -
                 1/p21*dotprod(dp21dqs,dq2s.dR2)   + 1/p22*dotprod(dp22dqs,dq2s.dR2))
detax2.dR1 <- (1/p12*dotprod(dp12dqs, dq1s.dR1) - 1/p13*dotprod(dp13dqs,dq1s.dR1) -
                 1/p22*dotprod(dp22dqs, dq2s.dR1) + 1/p23*dotprod(dp23dqs,dq2s.dR1))
detax2.dR2 <- (1/p12*dotprod(dp12dqs, dq1s.dR2) - 1/p13*dotprod(dp13dqs,dq1s.dR2) -
                 1/p22*dotprod(dp22dqs, dq2s.dR2) + 1/p23*dotprod(dp23dqs,dq2s.dR2))
# detax = Nx N^{-1} deta
Nx <- matrix(c(detax1.dR1,detax1.dR2,detax2.dR1,detax2.dR2),nc=2,byr=T)
N <- matrix(c(deta1.dR1,deta1.dR2,deta2.dR1,deta2.dR2),nc=2,byr=T)
detax.deta <- Nx %*% solve(N)

detax.deta
}


cc2x3.calc.detax.deta.opt <- function(D,p,q,p11,p12,p13,p21,p22,p23,q0,R,nu,c1.opt,
                                      verbose=FALSE){
# derivatives are taken with prevalence, nu, fixed
# version 12/7/2011 when using c1.opt*etax_1 + (1-c1.opt)*etax_2 as the estimator
#browser()
dq0.dR <- -2*q*q0/(q*R-q+1)
deta.dR <- (dq0.dR*(R^2-R)-q0+1)/((q0^2-q0)*R^2+(1-q0)*R)
    dp11dq11 <- (D + (1 - p) * (1 - q))^2/(1 - q)^2
    dp11dq12 <- ((1 - p) * q - D) * (D + (1 - p) * (1 - q))/((1 - 
        q) * q)
    dp11dq13 <- ((1 - p) * q - D)^2/q^2
    dp11dqs <- c(dp11dq11, dp11dq12, dp11dq13)
    dp12dq11 <- (D + (1 - p) * (1 - q))^2/(1 - q)^2
    dp12dq12 <- ((1 - p) * q - D) * (D + (1 - p) * (1 - q))/((1 - 
        q) * q)
    dp12dq13 <- ((1 - p) * q - D)^2/q^2
    dp12dqs <- c(dp12dq11, dp12dq12, dp12dq13)
    dp13dq11 <- (p * (1 - q) - D)^2/(1 - q)^2
    dp13dq12 <- (p * (1 - q) - D) * (D + p * q)/((1 - q) * q)
    dp13dq13 <- (D + p * q)^2/q^2
    dp13dqs <- c(dp13dq11, dp13dq12, dp13dq13)
    dp21dq21 <- (D + (1 - p) * (1 - q))^2/(1 - q)^2
    dp21dq22 <- ((1 - p) * q - D) * (D + (1 - p) * (1 - q))/((1 - 
        q) * q)
    dp21dq23 <- ((1 - p) * q - D)^2/q^2
    dp21dqs <- c(dp21dq21, dp21dq22, dp21dq23)
    dp22dq21 <- (D + (1 - p) * (1 - q))^2/(1 - q)^2
    dp22dq22 <- ((1 - p) * q - D) * (D + (1 - p) * (1 - q))/((1 - 
        q) * q)
    dp22dq23 <- ((1 - p) * q - D)^2/q^2
    dp22dqs <- c(dp22dq21, dp22dq22, dp22dq23)
    dp23dq21 <- (p * (1 - q) - D)^2/(1 - q)^2
    dp23dq22 <- (p * (1 - q) - D) * (D + p * q)/((1 - q) * q)
    dp23dq23 <- (D + p * q)^2/q^2
    dp23dqs <- c(dp23dq21, dp23dq22, dp23dq23)
#    nu <- q^2 * q0 * R^2 + 2 * (1 - q) * q * q0 * R + (1 - q)^2 * 
#        q0

dq11.dR <- -(1-q)^2*dq0.dR/(1-nu)
dq12.dR <- -2*(1-q)*q*(dq0.dR*R + q0)/(1-nu)
dq13.dR <- -q^2*(dq0.dR*R^2 + 2*q0*R)/(1-nu)
dq21.dR <- (1-q)^2*dq0.dR/nu
dq22.dR <- 2*(1-q)*q*(dq0.dR*R+q0)/nu
dq23.dR <- q^2*(dq0.dR*R^2+2*q0*R)/nu

dq1s.dR <- c(dq11.dR,dq12.dR,dq13.dR)
dq2s.dR <- c(dq21.dR,dq22.dR,dq23.dR)

#detax.deta <- 1/2*(1/p11*dotprod(dp11dqs, dq1s.dR) - 1/p13*dotprod(dp13dqs,dq1s.dR) -
#                   1/p21*dotprod(dp21dqs,dq2s.dR)   + 1/p23*dotprod(dp23dqs,dq2s.dR))

detax1.deta1 <- (1/p11*dotprod(dp11dqs,dq1s.dR) - 1/p12*dotprod(dp12dqs,dq1s.dR) -
                 1/p21*dotprod(dp21dqs,dq2s.dR)   + 1/p22*dotprod(dp22dqs,dq2s.dR))
detax2.deta2 <- (1/p12*dotprod(dp12dqs, dq1s.dR) - 1/p13*dotprod(dp13dqs,dq1s.dR) -
                 1/p22*dotprod(dp22dqs, dq2s.dR) + 1/p23*dotprod(dp23dqs,dq2s.dR))
detax.deta <-  c1.opt*detax1.deta1 + (1-c1.opt)*detax2.deta2
detax.deta <- detax.deta/deta.dR

detax.deta
}



cc2x3.calc.detax.deta <- function(D,p,q,p11,p13,p21,p23,q0,R,nu,verbose=FALSE){
# derivatives are taken with prevalence, nu, fixed
#browser()

dq0.dR <- -2*q*q0/(q*R-q+1)

deta.dR <- (dq0.dR*(R^2-R)-q0+1)/((q0^2-q0)*R^2+(1-q0)*R)


# 'diff(p[11],q11,1) = (D+(1-p)*(1-q))^2/(1-q)^2$
dp11dq11 <- (D+(1-p)*(1-q))^2/(1-q)^2
# 'diff(p[11],q12,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)$
dp11dq12 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
# 'diff(p[11],q13,1) = ((1-p)*q-D)^2/q^2$
dp11dq13 <-  ((1-p)*q-D)^2/q^2
dp11dqs <- c(dp11dq11,dp11dq12,dp11dq13)

# 'diff(p[13],q11,1) = (D+p*q)^2/(1-q)^2$
#!! dp13dq11 <- (D+p*q)^2/(1-q)^2
dp13dq11 <- (p*(1-q)-D)^2/(1-q)^2
# 'diff(p[13],q12,1) = (p*(1-q)-D)*(D+p*q)/((1-q)*q)$
dp13dq12 <- (p*(1-q)-D)*(D+p*q)/((1-q)*q)
# 'diff(p[13],q13,1) = (D+p*q)^2/q^2$
dp13dq13 <- (D+p*q)^2/q^2
dp13dqs <- c(dp13dq11,dp13dq12,dp13dq13)

# 'diff(p[21],q21,1) = (D+(1-p)*(1-q))^2/(1-q)^2$
dp21dq21 <-  (D+(1-p)*(1-q))^2/(1-q)^2
# 'diff(p[21],q22,1) = ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)$
dp21dq22 <- ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q)
# 'diff(p[21],q23,1) = ((1-p)*q-D)^2/q^2$
dp21dq23 <- ((1-p)*q-D)^2/q^2
dp21dqs <- c(dp21dq21,dp21dq22,dp21dq23)

#! 'diff(p[23],q21,1) = (D+p*q)^2/(1-q)^2$ 
#! dp23dq21 <- (D+p*q)^2/(1-q)^2
dp23dq21 <- (p*(1-q)-D)^2/(1-q)^2
# 'diff(p[23],q22,1) = (p*(1-q)-D)*(D+p*q)/((1-q)*q)$
dp23dq22 <- (p*(1-q)-D)*(D+p*q)/((1-q)*q)
# 'diff(p[23],q23,1) = (D+p*q)^2/q^2$
dp23dq23 <- (D+p*q)^2/q^2
dp23dqs <- c(dp23dq21,dp23dq22,dp23dq23)


if(F){ # values from  latest maxima output
dp11dqs.new <- c( (D+(1-p)*(1-q))^2/(1-q)^2,
                 ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q),
                 ((1-p)*q-D)^2/q^2 )
dp13dqs.new <- c( (p*(1-q)-D)^2/(1-q)^2,
                  (p*(1-q)-D)*(D+p*q)/((1-q)*q),
                  (p*(1-q)-D)*(D+p*q)/((1-q)*q) )
dp21dqs.new <- c( (D+(1-p)*(1-q))^2/(1-q)^2,
                  ((1-p)*q-D)*(D+(1-p)*(1-q))/((1-q)*q),
                  ((1-p)*q-D)^2/q^2 )
dp23dqs.new <- c( (p*(1-q)-D)^2/(1-q)^2,
                  (p*(1-q)-D)*(D+p*q)/((1-q)*q),
                  (p*(1-q)-D)*(D+p*q)/((1-q)*q) )
}

nu <-  q^2*q0*R^2+2*(1-q)*q*q0*R+(1-q)^2*q0

dq11.dR <- -(1-q)^2*dq0.dR/(1-nu)
dq12.dR <- -2*(1-q)*q*(dq0.dR*R + q0)/(1-nu)
dq13.dR <- -q^2*(dq0.dR*R^2 + 2*q0*R)/(1-nu)
dq21.dR <- (1-q)^2*dq0.dR/nu
dq22.dR <- 2*(1-q)*q*(dq0.dR*R+q0)/nu
dq23.dR <- q^2*(dq0.dR*R^2+2*q0*R)/nu

dq1s.dR <- c(dq11.dR,dq12.dR,dq13.dR)
dq2s.dR <- c(dq21.dR,dq22.dR,dq23.dR)

detax.deta <- 1/2*(1/p11* dotprod(dp11dqs, dq1s.dR) - 1/p13*dotprod(dp13dqs,dq1s.dR) -
                   1/p21*dotprod(dp21dqs,dq2s.dR)   + 1/p23*dotprod(dp23dqs,dq2s.dR))
detax.deta <- detax.deta/deta.dR

detax.deta
}

dotprod <- function(x,y){
#  if(length(x) !=length(y))stop("incompatible lengths")
  c(crossprod(x,y))
}


# trace(ccdesign2x2.ld.power.dominant,browser)
if(F){
cc.design(B=1e6,Dprime=0.5,p=0.3,q=0.05,OR=1.2,prevalence=1e-5,model="dominant",
          power=0.8,n.cases=1000,n.controls=1000)


}


ccdesign2x2.ld.power.dominant <- function(B,OR,D,p,q, baseline.risk, Dprime=NULL, R=NULL,
           prevalence=NULL,n.cases,n.controls,type="dominant",a=1,sigma2.eta=NULL, verbose=FALSE){
Dmax <- min((1-p)*q,(1-q)*p)
Dmin <- max(-p*q,-(1-p)*(1-q))
#
if(missing(D) && !missing(Dprime)){
  if(Dprime > 0) D <- Dprime*Dmax else D <-  (-Dprime)*Dmin
}  
#
if(missing(D) && missing(Dprime))Dprime <- 1
if(missing(D)){
  D <- Dmax
}else{
  if(D > Dmax || D < Dmin)
    stop("must have  max(-p*q,-(1-p)*(1-q)) < D < min((1-p)*q,(1-q)*p)")
}
#
if(missing(baseline.risk) && missing(prevalence))stop("must give baseline.risk or prevalence")
if(missing(OR) && missing(R))stop("must give odds ratio or relative risk (R)")
#
if(!missing(baseline.risk))q0 <- baseline.risk else q0 <- NULL
if(!missing(prevalence))nu <- prevalence else nu <- NULL
#
if(!missing(R)){
  if(missing(prevalence)){
    if(missing(baseline.risk))stop("must give baseline.risk or prevalence")
    nu <- q0*( (q^2+2*q*(1-q))*R + (1-q)^2)
  }else{
    q0 <- nu/( (q^2+2*q*(1-q))*R + (1-q)^2)
  }
  OR <- (q0-1)*R/(q0*R-1)
}
#
eta <- log(OR)
#
# [RR = OR/(q0*OR-q0+1)]$
if(missing(R)){
if(missing(baseline.risk)){
  if(missing(prevalence))stop("need to give baseline.risk or prevalence")
  nu <- prevalence
#########
#
# maxima solution for q0, dominant model:
# 
# [q0 = -(sqrt((q^4-4*q^3+(2*nu+4)*q^2-4*nu*q+nu^2)*OR^2
#               +(-2*q^4+8*q^3-10*q^2+4*q-2*nu^2+2*nu)*OR+q^4-4*q^3+(6-2*nu)*q^2
#               +(4*nu-4)*q+nu^2-2*nu+1)
#     +(-q^2+2*q-nu)*OR+q^2-2*q+nu+1)
#     /((2*q^2-4*q+2)*OR-2*q^2+4*q-2),
#  q0 = (sqrt((q^4-4*q^3+(2*nu+4)*q^2-4*nu*q+nu^2)*OR^2
#              +(-2*q^4+8*q^3-10*q^2+4*q-2*nu^2+2*nu)*OR+q^4-4*q^3+(6-2*nu)*q^2
#              +(4*nu-4)*q+nu^2-2*nu+1)
#     +(q^2-2*q+nu)*OR-q^2+2*q-nu-1)
#     /((2*q^2-4*q+2)*OR-2*q^2+4*q-2)]$
# 
#########

q0.soln1 <- -(sqrt((q^4-4*q^3+(2*nu+4)*q^2-4*nu*q+nu^2)*OR^2
               +(-2*q^4+8*q^3-10*q^2+4*q-2*nu^2+2*nu)*OR+q^4-4*q^3+(6-2*nu)*q^2
               +(4*nu-4)*q+nu^2-2*nu+1)
                +(-q^2+2*q-nu)*OR+q^2-2*q+nu+1) / ((2*q^2-4*q+2)*OR-2*q^2+4*q-2)
q0.soln2 <- (sqrt((q^4-4*q^3+(2*nu+4)*q^2-4*nu*q+nu^2)*OR^2
              +(-2*q^4+8*q^3-10*q^2+4*q-2*nu^2+2*nu)*OR+q^4-4*q^3+(6-2*nu)*q^2
              +(4*nu-4)*q+nu^2-2*nu+1)
             +(q^2-2*q+nu)*OR-q^2+2*q-nu-1)   /((2*q^2-4*q+2)*OR-2*q^2+4*q-2)
  if(q0.soln1 < 0 &&   q0.soln2 < 0) stop("no positive solution for baseline risk (q0)")
  if(q0.soln1 > 0 &&   q0.soln2 > 0){
    cat(paste("1st solution: ",q0.soln1, "2nd solution: ", q0.soln2,"\n"))
    warning("two positive solutions for baseline risk (q0)")
  }    
  if(q0.soln1 > 0){
    q0 <- q0.soln1
  }else{
    if(q0.soln2 > 0){
      q0 <- q0.soln2
    }
  }
#
# [R = OR/(q0*OR-q0+1)]$
#
    R <- OR/(q0*OR-q0+1)
}else{
  if(missing(prevalence)){
    q0 <- baseline.risk
    R <- OR/(q0*OR-q0+1)
#
# NB not eq (A10) for dominant model 
#
# nu = (q  + 2 (1 - q) q) q0 R + (1 - q)^2 q0 
    nu <- q0*( (q^2+2*q*(1-q))*R + (1-q)^2)
  }
}
}

# [R = OR/(q0*OR-q0+1)]$
R <- OR/(q0*OR-q0+1)

verify <- TRUE
if(verify){
  nu1 <- q0*( (q^2+2*q*(1-q))*R + (1-q)^2)
  err <- nu - nu1
  if(abs(err)>1e-6){
    cat(paste("nu =",nu, "nu1 = q0*(q^2*R + 1-q^2)",nu1, "error = nu - nu1 =", err, "\n"))
    stop("error in calculating q0,R")
  }
}

#
# [q11 = (1-q)^2*(1-q0)/(1-nu),q12 = (q^2+2*(1-q)*q)*(1-q0*R)/(1-nu),
#  q21 = (1-q)^2*q0/nu,q22 = (q^2+2*(1-q)*q)*q0*R/nu]$
#


q11 <- (1-q)^2*(1-q0)/(1-nu)
q12 <- (2*q-q^2)*(1-q0*R)/(1-nu)
q21 <- (1-q)^2*q0/nu
q22 <- (2*q-q^2)*q0*R/nu


if(missing(D)){
  p11 <- q11; p12 <- q12; p21 <- q21; p22 <- q22
}else{


# 
# subst(2*q - q^2, 1 - (1-q)^2, [p11simp2,p12simp2,p21simp2,p22simp2]);
# grind(%);
# [q11*(D+(1-p)*(1-q))^2/(1-q)^2-q12*((D+(1-p)*(1-q))^2-(1-p)^2)/(2*q-q^2),
#  q12*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(2*q-q^2)
#   -q11*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2,
#  q21*(D+(1-p)*(1-q))^2/(1-q)^2-q22*((D+(1-p)*(1-q))^2-(1-p)^2)/(2*q-q^2),
#  q22*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(2*q-q^2)
#   -q21*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2]$
# 
# 

p11 <-   q11*(D+(1-p)*(1-q))^2/(1-q)^2            - q12*((D+(1-p)*(1-q))^2-(1-p)^2)/(2*q-q^2)
p12 <-  -q11*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2  + q12*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(2*q-q^2)
p21 <-   q21*(D+(1-p)*(1-q))^2/(1-q)^2            - q22*((D+(1-p)*(1-q))^2-(1-p)^2)/(2*q-q^2)
p22 <-  -q21*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2  + q22*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(2*q-q^2) 

#browser()

if(F){ # values from latest maxima output.
p11.new <- q11*(D+(1-p)*(1-q))^2/(1-q)^2-q12*((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2)
p12.new <- q12*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2) -q11*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2
p21.new <- q21*(D+(1-p)*(1-q))^2/(1-q)^2-q22*((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2)
p22.new <- q22*((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2) -q21*((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2
matrix(c(p11,p12,p21,p22,p11.new,p12.new,p21.new,p22.new),nc=4,byr=T)
apply(matrix(c(p11,p12,p21,p22),nc=2,byr=T),1,sum)
apply(matrix(c(q11,q12,q21,q22),nc=2,byr=T),1,sum)
}
} 

if(verbose){
cat(paste("D =", D,", Dmin =",Dmin, "Dmax =", Dmax,"\n"))
cat(paste("q0 =",q0,", R =",R, ", nu =", nu,", OR =", OR, "\n"))
cat(paste("q11 =", q11, ", q12 =", q12, ", q21 =", q21, ", q22 =", q22,"\n"))
cat(paste("p11 =", p11, ", p12 =", p12, ", p21 =", p21, ", p22 =", p22,"\n"))
}

#browser()

n1 <- n.controls
n2 <- n.cases
n <- n1+n2
sigma1 <- sqrt(n*(1/(n1*p11) + 1/(n1*p12) + 1/(n2*p21) + 1/(n2*p22)))
# calculate `a' (equivalent sample size for prior) if missing
if(missing(D)){
    detax.deta <- 1
}else{
    detax.deta <- cc2x2.calc.detax.deta.dominant(D,p,q,p11,p12,p21,p22,q0,R,nu)
}

if(missing(a) && !missing(sigma2.eta)){
  if(missing(sigma2.eta))stop("must give \"a\" or \"sigma2.eta\"")
#  sigma2.etax <- detax.deta^2*sigma2.eta
  a <- sigma1^2/sigma2.eta
}
ax <- a/detax.deta^2  

# critical value of Z-statistic
etax <- log(p11) - log(p12) - log(p21) + log(p22)
Zc <- calc.Zc.ABF(B,n=n.cases+n.controls,a=ax)

# corresponding alpha value
# alphac <- 2*(1-pnorm(sqrt(n)*Zc))
alphac <- 1 - pchisq(n*Zc^2,1)
alphac

# equivalent non-centrality parameter in chi-squared test
ncp <- n*(etax/sigma1)^2
ncp

if(F){ 
power <- 1 - pchisq(qchisq(1-alphac,1),df=1,ncp=ncp)
power
}


## neglecting contribution where statistic and value of opposite signs
power <- 1 - pnorm(sqrt(n)*(Zc - etax/sigma1))

if(verbose){
  print(c(n=n,q0=q0,R=R,nu=nu,OR=OR,ax=ax,etax=etax,sigma1=sigma1,
          etaxoversigma1=etax/sigma1, ncp=ncp, alphac=alphac,
          OR=OR,B=B,Zc=Zc,power=power))
}

attr(power,"model") <- "dominant"
attr(power,"n.cases") <- n.cases
attr(power,"n.controls") <- n.controls
attr(power,"baseline.risk") <- q0
attr(power,"prevalence") <- nu
attr(power,"relative.risk") <- R
attr(power,"odds-ratio") <- OR
attr(power,"Dprime") <- Dprime
attr(power,"p") <- p
attr(power,"q") <- q
attr(power,"ncp") <- ncp
attr(power,"B") <- B
attr(power,"alphac") <- alphac
ps <- rbind(c(p11,p12),c(p21,p22))
dimnames(ps) <- list(c("control","case"),c("bb","Bb||BB"))
attr(power,"ps") <- ps

power
}


cc2x2.calc.detax.deta.dominant <- function(D,p,q,p11,p12,p21,p22,q0,R,nu,verbose=FALSE){
#
# ['diff(q0,R,1) = -q^2*q0/(q^2*R-q^2+1)]$
#
# eq (A37)
dq0.dR <- -(q^2-2*q)*q0/((q^2-2*q)*R-q^2+2*q-1)
dq0.dR <- (2*q - q^2)*q0/((2*q - q^2)*R + q^2 - 2*q + 1)

#
#['diff(eta,R,1) = ('diff(q0,R,1)*R^2-'diff(q0,R,1)*R-q0+1)/((q0^2-q0)*R^2+(1-q0)*R)]$
#
# eq (A13), all models
deta.dR <- (dq0.dR*(R^2-R) - q0 + 1)/((q0^2 - q0)*R^2 + (1-q0)*R)

#
#
# ['diff(p11,q11,1) = (D+(1-p)*(1-q))^2/(1-q)^2,
#  'diff(p11,q12,1) = -((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2),
#  'diff(p12,q11,1) = -((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2,
#  'diff(p12,q12,1) = ((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2)]$
# ['diff(p21,q21,1) = (D+(1-p)*(1-q))^2/(1-q)^2,
#  'diff(p21,q22,1) = -((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2),
#  'diff(p22,q21,1) = -((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2,
#  'diff(p22,q22,1) = ((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2)]$
# 
# 

dp11.dq11 <- (D+(1-p)*(1-q))^2/(1-q)^2
dp11.dq12 <- -((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2)
dp12.dq11 <- -((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2
dp12.dq12 <- ((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2)

dp21.dq21 <- (D+(1-p)*(1-q))^2/(1-q)^2
dp21.dq22 <- -((D+(1-p)*(1-q))^2-(1-p)^2)/(1-(1-q)^2)
dp22.dq21 <- -((D+(1-p)*(1-q))^2-(1-q)^2)/(1-q)^2
dp22.dq22 <- ((D+(1-p)*(1-q))^2-(1-q)^2-(1-p)^2+1)/(1-(1-q)^2)


#
# ['diff(q11,R,1) = -(1-q)^2*'diff(q0,R,1)/(1-nu),
#  'diff(q12,R,1) = (q^2+2*(1-q)*q)*(-'diff(q0,R,1)*R-q0)/(1-nu),
#  'diff(q21,R,1) = (1-q)^2*'diff(q0,R,1)/nu,
#  'diff(q22,R,1) = (q^2+2*(1-q)*q)*'diff(q0,R,1)*R/nu+(q^2+2*(1-q)*q)*q0/nu]$

dq11.dR <- -(1-q)^2*dq0.dR/(1-nu)
dq12.dR <- (q^2+2*(1-q)*q)*(-dq0.dR*R-q0)/(1-nu)
dq21.dR <- (1-q)^2*dq0.dR/nu
dq22.dR <- ((q^2+2*(1-q)*q)*dq0.dR*R+(q^2+2*(1-q)*q)*q0)/nu

# eq (37)
A <- 1/p11*(dp11.dq11*dq11.dR + dp11.dq12*dq12.dR) - 1/p12*(dp12.dq11*dq11.dR + dp12.dq12*dq12.dR) -
     1/p21*(dp21.dq21*dq21.dR + dp21.dq22*dq22.dR) + 1/p22*(dp22.dq21*dq21.dR + dp22.dq22*dq22.dR) 

detax.deta <- A/deta.dR

if(verbose){
cat(paste("p11 =",p11, ", p12 =",p12,", p21 = ",p21,", p22 =",p22,"\n"))
cat(paste("dp11.dq11 =",dp11.dq11,", dp11.dq12 =",dp11.dq12,", dp12.dq11 =",dp12.dq11,", dp12.dq12 =",
            dp12.dq12,"\n"))
cat(paste("dp21.dq21 =",dp21.dq21,", dp21.dq22 =",dp21.dq22,", dp22.dq21 =",dp22.dq21,", dp22.dq22 =",
            dp22.dq22,"\n"))
cat(paste("A =",A,", deta.dR =", deta.dR,", detax.deta =", detax.deta,"\n"))

}

detax.deta
}


ccdesign2x2.ld.power.recessive <- function(B,OR,D,p,q, baseline.risk, Dprime=NULL, R=NULL,
                 prevalence=NULL,n.cases, n.controls,model="recessive",a=1,sigma2.eta=NULL, verbose=FALSE){
#browser()
Dmax <- min((1-p)*q,(1-q)*p)
Dmin <- max(-p*q,-(1-p)*(1-q))
if(missing(D) && !missing(Dprime)){
  if(Dprime > 0) D <- Dprime*Dmax else D <-  (-Dprime)*Dmin
}  
#
if(missing(D) && missing(Dprime))Dprime <- 1
if(missing(D)){
  D <- Dmax
}else{
  if(D > Dmax || D < Dmin)
    stop("must have  max(-p*q,-(1-p)*(1-q)) < D < min((1-p)*q,(1-q)*p)")
}
#
if(missing(baseline.risk) && missing(prevalence))stop("must give baseline.risk or prevalence")
if(missing(OR) && missing(R))stop("must give odds ratio or relative risk (R)")
#
if(!missing(baseline.risk))q0 <- baseline.risk else q0 <- NULL
if(!missing(prevalence))nu <- prevalence else nu <- NULL
#
if(!missing(R)){
  if(missing(prevalence)){
    if(missing(baseline.risk))stop("must give baseline.risk or prevalence")
    nu <- q0*(q^2*R + 1-q^2)
  }else{
    q0 <- nu/(q^2*R + 1-q^2)
  }
  OR <- (q0-1)*R/(q0*R-1)
}
#
eta <- log(OR)
# [RR = OR/(q0*OR-q0+1)]$
if(missing(R)){
if(missing(baseline.risk)){
  if(missing(prevalence))stop("need to give baseline.risk or prevalence")
  nu <- prevalence
#########
#
# maxima solution for q0, recessive model:
# [q0 = -(sqrt((q^4-2*nu*q^2+nu^2)*OR^2+(-2*q^4+2*q^2-2*nu^2+2*nu)*OR+q^4
#                                      +(2*nu-2)*q^2+nu^2-2*nu+1)
#     +(nu-q^2)*OR+q^2-nu-1)
#     /((2*q^2-2)*OR-2*q^2+2),
#  q0 = (sqrt((q^4-2*nu*q^2+nu^2)*OR^2+(-2*q^4+2*q^2-2*nu^2+2*nu)*OR+q^4
#                                     +(2*nu-2)*q^2+nu^2-2*nu+1)
#     +(q^2-nu)*OR-q^2+nu+1)
#     /((2*q^2-2)*OR-2*q^2+2)]$
# 
#########
  q0.soln1 <- -(sqrt((q^4-2*nu*q^2+nu^2)*OR^2+(-2*q^4+2*q^2-2*nu^2+2*nu)*OR+q^4
               +(2*nu-2)*q^2+nu^2-2*nu+1) +(nu-q^2)*OR+q^2-nu-1)/((2*q^2-2)*OR-2*q^2+2)
  q0.soln2 <- (sqrt((q^4-2*nu*q^2+nu^2)*OR^2+(-2*q^4+2*q^2-2*nu^2+2*nu)*OR+q^4
               +(2*nu-2)*q^2+nu^2-2*nu+1) +(q^2-nu)*OR-q^2+nu+1)/((2*q^2-2)*OR-2*q^2+2)
  if(q0.soln1 < 0 &&   q0.soln2 < 0) stop("no positive solution for baseline risk (q0)")
  if(q0.soln1 > 0 &&   q0.soln2 > 0){
    cat(paste("1st solution: ",q0.soln1, "2nd solution: ", q0.soln2,"\n"))
    warning("two positive solutions for baseline risk (q0)")
  }    
  if(q0.soln1 > 0){
    q0 <- q0.soln1
  }else{
    if(q0.soln2 > 0){
      q0 <- q0.soln2
    }
  }
#
# [R = OR/(q0*OR-q0+1)]$
#
    R <- OR/(q0*OR-q0+1)
}else{
  if(missing(prevalence)){
    q0 <- baseline.risk
    R <- OR/(q0*OR-q0+1)
#
# NB not eq (A10) for recessive model 
#
# [nu = q^2*q0*R+(1-q^2)*q0]$
#
#  nu <- q^2*q0*R+(1-q^2)*q0
    nu <- q0*(q^2*R + 1-q^2)
  }
}
}
# [R = OR/(q0*OR-q0+1)]$
R <- OR/(q0*OR-q0+1)
verify <- TRUE
if(verify){
  nu1 <- q0*(q^2*R + 1-q^2)
  err <- nu - nu1
  if(abs(err)>1e-6){
    cat(paste("nu =",nu, "nu1 = q0*(q^2*R + 1-q^2)",nu1, "error = nu - nu1 =", err, "\n"))
    stop("error in calculating q0,R")
  }
}
#
# [q11 = (2*(1-q)*q+(1-q)^2)*(1-q0)/(1-nu),q12 = q^2*(1-q0*R)/(1-nu),
#  q21 = (2*(1-q)*q+(1-q)^2)*q0/nu,q22 = q^2*q0*R/nu]$
#
q11 <- (1-q^2)*(1-q0)/(1-nu)
q12 <- q^2*(1-q0*R)/(1-nu)
q21 <- (1-q^2)*q0/nu
q22 <- q^2*q0*R/nu
#
if(missing(D)){
  p11 <- q11; p12 <- q12; p21 <- q21; p22 <- q22
}else{
##############
#
# [p11 = q11*(((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))
#            /(2*(1-q)*q+(1-q)^2)
#            +(2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
#             +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))
#             /(2*(1-q)*q+(1-q)^2))
#      +q12*(2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2),
#  p12 = q12*(D+p*q)^2/q^2+q11*(2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)
#                          /(2*(1-q)*q+(1-q)^2),
#  p21 = q21*(((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))
#            /(2*(1-q)*q+(1-q)^2)
#            +(2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
#             +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))
#             /(2*(1-q)*q+(1-q)^2))
#      +q22*(2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2),
#  p22 = q22*(D+p*q)^2/q^2+q21*(2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)
#                          /(2*(1-q)*q+(1-q)^2)]$
# 
###############
#
# simplified equations for pij
#
p11 <-  q11/(1-q^2) * ((D+p*q)^2 + 1 - p^2 - q^2) -  q12/q^2 * ((D+p*q)^2 - q^2)
p12 <- -q11/(1-q^2) * ((D+p*q)^2 - p^2)           +  q12/q^2 * (D+p*q)^2
p21 <-  q21/(1-q^2) * ((D+p*q)^2 + 1 - p^2 - q^2) -  q22/q^2 * ((D+p*q)^2 - q^2)
p22 <- -q21/(1-q^2) * ((D+p*q)^2 - p^2)           +  q22/q^2 * (D+p*q)^2
#
#browser()
#
if(F){ # check against values from latest maxima output, OK.
p11.new <- q11*((D+p*q)^2-q^2-p^2+1)/(1-q^2)-q12*((D+p*q)^2-q^2)/q^2
p12.new <- q12*(D+p*q)^2/q^2-q11*((D+p*q)^2-p^2)/(1-q^2)
p21.new <- q21*((D+p*q)^2-q^2-p^2+1)/(1-q^2)-q22*((D+p*q)^2-q^2)/q^2
p22.new <- q22*(D+p*q)^2/q^2-q21*((D+p*q)^2-p^2)/(1-q^2)
matrix(c(p11,p12,p21,p22,p11.new,p12.new,p21.new,p22.new),nc=4,byr=T)
apply(matrix(c(p11,p12,p21,p22),nc=2,byr=T),1,sum)
apply(matrix(c(q11,q12,q21,q22),nc=2,byr=T),1,sum)
}
}

if(verbose){
cat(paste("D =", D,", Dmin =",Dmin, "Dmax =", Dmax,"\n"))
cat(paste("q0 =",q0,", R =",R, ", nu =", nu,", OR =", OR, "\n"))
cat(paste("q11 =", q11, ", q12 =", q12, "\nq21 =", q21, ", q22 =", q22,"\n"))
cat(paste("p11 =", p11, ", p12 =", p12, ", p21 =", p21, ", p22 =", p22,"\n"))
}


n1 <- n.controls
n2 <- n.cases
n <- n1+n2
sigma1 <- sqrt(n*(1/(n1*p11) + 1/(n1*p12) + 1/(n2*p21) + 1/(n2*p22)))
# calculate `a' (equivalent sample size for prior) if missing
if(missing(D)){
    detax.deta <- 1
  }else{
    detax.deta <- cc2x2.calc.detax.deta.recessive(D,p,q,p11,p12,p21,p22,q0,R,nu)
  }
if(missing(a) && !missing(sigma2.eta)){
  if(missing(sigma2.eta))stop("must give \"a\" or \"sigma2.eta\"")
#  sigma2.etax <- detax.deta^2*sigma2.eta
  a <- sigma1^2/sigma2.eta
}
ax <- a/detax.deta^2  
#
# critical value of Z-statistic
etax <- log(p11) - log(p12) - log(p21) + log(p22)
Zc <- calc.Zc.ABF(B,n=n.cases+n.controls,a=ax)

# corresponding alpha value
# alphac <- 2*(1-pnorm(sqrt(n)*Zc))
alphac <- 1 - pchisq(n*Zc^2,1)
alphac

# equivalent non-centrality parameter in chi-squared test
ncp <- n*(etax/sigma1)^2
ncp

if(F){ 
power <- 1 - pchisq(qchisq(1-alphac,1),df=1,ncp=ncp)
power
}

if(verbose){
  print(c(n=n,q0=q0,R=R,nu=nu,OR=OR,ax=ax,etax=etax,sigma1=sigma1,
          etaxoversigma1=etax/sigma1, ncp=ncp, alphac=alphac,
          OR=OR,B=B,Zc=Zc,power=power))
}

#
## neglecting contribution where statistic and value of opposite signs
if(F){
  power <- 1 - pnorm(sqrt(n)*(Zc - etax/sigma1))
}
# sqrt(n)Z | H1  ~= N(sqrt(n)*etax/sigma1,1)
# Or: n*Z^2 |H1 ~= chisquared(. ,df=1,ncp=n*(etax/sigma1)^2)
power <- 1 - pchisq(n*Zc^2, df=1,ncp=ncp)

attr(power,"model") <- "recessive"
attr(power,"n.cases") <- n.cases
attr(power,"n.controls") <- n.controls
attr(power,"prevalence") <- nu
attr(power,"baseline.risk") <- q0
attr(power,"relative.risk") <- R
attr(power,"odds-ratio") <- OR
attr(power,"Dprime") <- Dprime
attr(power,"p") <- p
attr(power,"q") <- q
attr(power,"ncp") <- ncp
attr(power,"B") <- B
attr(power,"alphac") <- alphac
ps <- rbind(c(p11,p12),c(p21,p22))
dimnames(ps) <- list(c("control","case"),c("bb||Bb","BB"))
attr(power,"ps") <- ps

power
}


cc2x2.calc.detax.deta.recessive <- function(D,p,q,p11,p12,p21,p22,q0,R,nu,verbose=FALSE){

#
# ['diff(q0,R,1) = -q^2*q0/(q^2*R-q^2+1)]$
#
# eq (A45)
dq0.dR <- -q^2*q0/(q^2*R-q^2+1)

#
#['diff(eta,R,1) = ('diff(q0,R,1)*R^2-'diff(q0,R,1)*R-q0+1)/((q0^2-q0)*R^2+(1-q0)*R)]$
#
# eq (A13), additive (!?) model
deta.dR <- (dq0.dR*(R^2-R) - q0 + 1)/((q0^2 - q0)*R^2 + (1-q0)*R)

###################
#
# ['diff(p11,q11,1) = ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))
#                   /(2*(1-q)*q+(1-q)^2)
#                   +(2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
#                    +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))
#                    /(2*(1-q)*q+(1-q)^2),
#  'diff(p11,q12,1) = 2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2,
#  'diff(p12,q11,1) = (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2),
#  'diff(p12,q12,1) = (D+p*q)^2/q^2]$
# ['diff(p21,q21,1) = ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))
#                   /(2*(1-q)*q+(1-q)^2)
#                   +(2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
#                    +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))
#                    /(2*(1-q)*q+(1-q)^2),
#  'diff(p21,q22,1) = 2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2,
#  'diff(p22,q21,1) = (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2),
#  'diff(p22,q22,1) = (D+p*q)^2/q^2]$
# 
###################


#
# ['diff(p11,q11,1) = ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2)
#             +((D+(1-p)*(1-q))*(D+p*q)+2*(p*(1-q)-D)*(D+(1-p)*(1-q))
#                                           +(p*(1-q)-D)*((1-p)*q-D))/
#                   (2*(1-q)*q+(1-q)^2)]
#
dp11.dq11.old <- ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2)+
                     ((D+(1-p)*(1-q))*(D+p*q)+2*(p*(1-q)-D)*(D+(1-p)*(1-q))
                         +(p*(1-q)-D)*((1-p)*q-D))/
                       (2*(1-q)*q+(1-q)^2)
# fix to Pr(AaIQq), change from eq (71)
dp11.dq11 <-     ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2) +
                     (2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
                       +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))/
                       (2*(1-q)*q+(1-q)^2)

#
# 'diff(p11,q12,1) = 2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2,
#
dp11.dq12 <-     2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2

#
# 'diff(p12,q11,1) = (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2),
#
dp12.dq11 <-  (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2)
#
# 'diff(p12,q12,1) = (D+p*q)^2/q^2]$
#
dp12.dq12 <- (D+p*q)^2/q^2

#
# ['diff(p21,q21,1) = ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))
#                  /(2*(1-q)*q+(1-q)^2)
#                  +((D+(1-p)*(1-q))*(D+p*q)+2*(p*(1-q)-D)*(D+(1-p)*(1-q))
#                                           +(p*(1-q)-D)*((1-p)*q-D))
#                   /(2*(1-q)*q+(1-q)^2),
#
dp21.dq21.old <- ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2) + 
             ((D+(1-p)*(1-q))*(D+p*q)+2*(p*(1-q)-D)*(D+(1-p)*(1-q))
                    +(p*(1-q)-D)*((1-p)*q-D)) / (2*(1-q)*q+(1-q)^2)

dp21.dq21.old <- ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)) +
              (D+(1-p)*(1-q))*(D+p*q)+2*(p*(1-q)-D)*(D+(1-p)*(1-q))
                    +(p*(1-q)-D)*((1-p)*q-D) ) / (2*(1-q)*q+(1-q)^2)

# fix to Pr(AaIQq) -> change from eq ()
dp21.dq21 <-     ((D+(1-p)*(1-q))^2+2*((1-p)*q-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2) +
                     (2*((D+(1-p)*(1-q))*(D+p*q)+(p*(1-q)-D)*((1-p)*q-D))
                       +2*(p*(1-q)-D)*(D+(1-p)*(1-q)))/(2*(1-q)*q+(1-q)^2)
#
#  'diff(p21,q22,1) = 2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2,
#
if(F){
dp21.dq22 <-   2*((1-p)*q-D)*(D+p*q)/q^2+((1-p)*q-D)^2/q^2
dp21.dq22
}
dp21.dq22 <-   (2*((1-p)*q-D)*(D+p*q) + ((1-p)*q-D)^2)/q^2
#
#  'diff(p22,q21,1) = (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2),
#
dp22.dq21 <- (2*(p*(1-q)-D)*(D+p*q)+(p*(1-q)-D)^2)/(2*(1-q)*q+(1-q)^2)

#
#  'diff(p22,q22,1) = (D+p*q)^2/q^2]$
#
dp22.dq22 <- (D+p*q)^2/q^2

#
# 'diff(q11,R,1) = -(2*(1-q)*q+(1-q)^2)*'diff(q0,R,1)/(1-nu),
#

dq11.dR <-  -(2*(1-q)*q+(1-q)^2)*dq0.dR/(1-nu)
#
#  'diff(q12,R,1) = q^2*(-'diff(q0,R,1)*R-q0)/(1-nu),
#
dq12.dR <-  q^2*(-dq0.dR*R-q0)/(1-nu)
#
# 'diff(q21,R,1) = (2*(1-q)*q+(1-q)^2)*dq0.dR/nu,
#
dq21.dR <-  (2*(1-q)*q+(1-q)^2)*dq0.dR/nu
#
# 'diff(q22,R,1) = q^2*dq0.dR*R/nu+q^2*q0/nu]$
# 
#dq22.dR <- q^2*dq0.dR*R/nu+q^2*q0/nu
dq22.dR <- q^2/nu * (dq0.dR*R + q0)

# eq (61)
A <- 1/p11*(dp11.dq11*dq11.dR + dp11.dq12*dq12.dR) - 1/p12*(dp12.dq11*dq11.dR + dp12.dq12*dq12.dR) -
     1/p21*(dp21.dq21*dq21.dR + dp21.dq22*dq22.dR) + 1/p22*(dp22.dq21*dq21.dR + dp22.dq22*dq22.dR) 

detax.deta <- A/deta.dR

if(verbose){
cat(paste("p11 =",p11, ", p12 =",p12,", p21 = ",p21,", p22 =",p22,"\n"))
cat(paste("dp11.dq11 =",dp11.dq11,", dp11.dq12 =",dp11.dq12,", dp12.dq11 =",dp12.dq11,", dp12.dq12 =",
            dp12.dq12,"\n"))
cat(paste("dp21.dq21 =",dp21.dq21,", dp21.dq22 =",dp21.dq22,", dp22.dq21 =",dp22.dq21,", dp22.dq22 =",
            dp22.dq22,"\n"))
cat(paste("A =",A,", deta.dR =", deta.dR,", detax.deta =", detax.deta,"\n"))
}



detax.deta
}


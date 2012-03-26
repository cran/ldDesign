## ldDesign R package for design of experiments for association studies
## for detection of linkage disequilibrium
#
# Copyright (C) 2003 Roderick D. Ball
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# Or see www.gnu.org/copyleft/gpl.htlm


ld.power <- function(n,p,q,D,h2,phi,Bf,missclass.rate=0){
# function to calculate power of experiment to detect
# linkage equilibrium between a biallelic marker and QTL                      
# with a given Bayes factor  
# n : vector of sample sizes
# Bf : desired Bayes factor
alphas <- numeric(length(n))
powers <- numeric(length(n))
for(ii in seq(along=n)){
  alphas[ii] <- oneway.bf.alpha.required(n=n[ii],group.sizes=c(p^2,2*p*(1-p),(1-p)^2)*n[ii],Bf=Bf)
  powers[ii] <- luo.ld.power(n=n[ii],p=p,q=q,D=D,h2=h2,phi=phi,alpha=alphas[ii],print.it=FALSE,
                             missclass.rate=missclass.rate)
}
res <- cbind(n=n,power=powers)
attr(res,"parms") <- c(p=p,q=q,D=D,h2=h2,phi=phi,Bf=Bf,missclass.rate=missclass.rate)
res
}


ld.design <- function(p,q,D,h2,phi,Bf,power,nmin=50,nmax=100000,ninterp=50,missclass.rate=0,
                      print.it=FALSE){
logit <- function(p,tol=1e-6){p1 <- ifelse(p<tol,tol,ifelse(p>1-tol,1-tol,p));log(p1/(1-p1))}
ns <- exp(seq(log(nmin),log(nmax),length=ninterp))
power.res <- ld.power(ns,p,q,D,h2,phi,Bf,missclass.rate)
if(print.it)print(power.res)
n <- exp(approx(logit(power.res[,"power"]),log(ns),xout=logit(power))$y)
#browser()
n
}

luo.ld.power <- function(n, p, q, D, h2, phi, Vp=100, alpha, print.it=TRUE,  missclass.rate=0){
# function to calculate power of experiment to detect
# linkage equilibrium between a biallelic marker and QTL                      
# Anova model power calculations from Z.W. Luo(1998):
# Detecting linkage disequilibrium between a polymorphic marker locus 
# and a trait locus in natural populations,  Heredity 80, 198--208
# n: census population size
# p: allele frequency of  marker type M (2 types: Mm)                      
# q: allele frequency of  QTL type A (2 types: Aa)
# D: linkage disequilibrium between marker and QTL.
# h2: QTL heritability, i.e. proportion of total variance explained by the QTL
# phi: dominance ratio, if d,h are QTL additive and dominance effects 
#     then AA,Aa,aa -> effects d, h, -d, on phenotype, and h = phi*d
# Vp: total phenotypic variance (arbitrary)
# alpha: type I error rate (comparison-wise)
# for missclassificationwith rate m use : beta' <- (1 - m)*beta
#                                         omega' <- omega + beta - beta'
#                                         f2' <- (1-m)*f2 + m*f2.1
# where  f2.1[i,j]  = f1[i]  * f1.q[j], where f1.q are the overall
# QTL genotype frequencies, i.e. (q^2,2q(1-q),(1-q)^2) for AA,Aa,aa  
# We then calculate EMS.omega using beta', omega',f2'  etc.
# version for CRAN
# Changes: `method' option removed
calc.sig2 <- function(ems.beta,ems.omega,n,f1){
# variance components from mean squares
  df.beta <- length(f1) - 1
  ni <- f1 * n
  sig2.beta <- (ems.beta - ems.omega)/((n - sum(ni^2)/n)/(df.beta))
  sig2.omega <- ems.omega
  list(sig2.beta=sig2.beta,sig2.omega=sig2.omega)
}
calc.ems <- function(sig2.beta,sig2.omega,n,f1){
# expected mean squares from variance components
  df.beta <- length(f1) - 1
  ni <- f1 * n
  ems.beta <- 1/df.beta*(n - sum(ni^2)/n)*sig2.beta + sig2.omega
  ems.omega <- sig2.omega
  list(ems.beta=ems.beta,ems.omega=ems.omega)
}

method1.calc.ems <- function(beta,omega,n,G,q,f1,f2,Vp,sig2.error,
                                missclass.rate=0){
# calculate 
beta0 <- beta
omega0 <- omega
beta <- (1-missclass.rate)*beta0
omega <- omega0 + beta0 - beta
f2.0 <- f2

#browser()

f2.1 <- f1*matrix(rep(c(q^2,2*q*(1-q),(1-q)^2),3),ncol=3,byrow=TRUE)

f2 <- (1-missclass.rate)*f2.0 + missclass.rate*f2.1
  
I <- row(f2)
J <- col(f2)
su <- I < J
IJ <- cbind(c(I),c(J))
ijk <- expand.grid(I=1:3,J=1:3,K=1:3)
sjk <- ijk$J < ijk$K

ems.omega <-  1/(n-G)*(
    n*sum(f2[IJ]*omega[IJ]^2)  
-   sum(1/f1[I]*f2[IJ]*(1+(n-1)*f2[IJ])*omega[IJ]^2)  
-   sum(1/f1[ijk$I[sjk]]*2*(n-1)*f2[cbind(ijk$I[sjk],ijk$J[sjk])]*
    f2[cbind(ijk$I[sjk],ijk$K[sjk])]*omega[cbind(ijk$I[sjk],ijk$J[sjk])]*
        omega[cbind(ijk$I[sjk],ijk$K[sjk])])
)  + sig2.error

ems.beta <- 
    1/(G-1)*( (n-1)*Vp - (n-G)*ems.omega)

list(ems.beta=ems.beta,ems.omega=ems.omega)
}


if(missclass.rate > 1 || missclass.rate <0)
  stop("invalid missclassification rate: not in [0,1]")


G <- 3 # number of marker genotypes                        
nu1 <- 2
nu2 <- n - 1

# QTL values d,h,-d for AA,Aa,aa respectively
if(h2 < 0 || h2 > 1)stop("heritability not in [0,1]")
Vq <- h2*Vp

# Vq = (2q - 2q^2)d^2 + (4q^5 -6q^4 -6q^3 + 6q^2 -2q)\phi d^2 +
#       (2q -6q^2 +8q^3 - 4q^4)\phi^2 d^2
# Vq = (2q - 2q^2)d^2 + (-16q^3 + 12q^2 -4q)\phi d^2 +
#       (2q -6q^2 +8q^3 - 4q^4)\phi^2 d^2

d <- sqrt(Vq/((2*q - 2*q^2) + (8*q^3-12*q^2 + 4*q)*phi  +
           (2*q -6*q^2 +8*q^3 - 4*q^4)*phi^2 ))
h <- phi*d
Vq 
sig2.error <- Vp - Vq

# self test for Vq calculation
do.test.Vq <- FALSE
if(do.test.Vq){
nsim.Vq <- 1000000
test.Vq <- var(sample(c(d,h,-d),size=nsim.Vq,replace=TRUE,prob=c(q^2,2*q*(1-q),(1-q)^2)))
se.test.Vq <- sqrt(2*Vq^2/(nsim.Vq - 2))
se.test.Vq
if(abs(test.Vq - Vq ) > 6*se.test.Vq){
  cat("discrepancy between Vq, and sample test for calculated d,h\n")
  cat("Vq =",Vq, "test.Vq =", test.Vq,"\n")
}
}

# Marker/QTL genotypes 1,2,3 as MM,Mm,mm and AA,Aa,aa respectively
# Cf Table 1.

# marker genotype frequencies for MM,Mm,mm (f <- i in paper)
# p^2,2*p*(1-p),(1-p)^2
f1 <- c(p^2,2*p*(1-p),(1-p)^2)

# Q, R are frequencies of A where marker = M,m resp.
Q  <-  q + D/p
R  <-  q - D/(1-p)

tol <- 1.0e-4
if(abs(Q - 1) < tol) Q <- 1
if(abs(Q) < tol) Q <- 0
if(abs(R - 1) < tol) R <- 1
if(abs(R) < tol) R <- 0

if(Q > 1 || Q < 0 || R > 1 || R < 0){
  cat(paste(
    "warning: invalid value for disequilibrium: p =",p,"q =",q, "D =",D, 
    "Dmax =",min(p*(1-q),q*(1-p)), "Dmin =", max(-p*q,(q-1)*(1-p)),"\n"))
#  browser()
}
# marker/qtl genotype frequencies
# f2[i,j] is frequency for ith marker genotype with jth QTL genotype (f_{ij} in paper)
# Cf Table 1 of Luo.
f2  <- matrix( 
      c(p^2*Q^2,2*p^2*Q*(1-Q), p^2*(1-Q)^2,
        2*p*(1-p)*Q*R, 2*p*(1-p)*(Q+R - 2*Q*R),2*p*(1-p)*(1-Q)*(1-R),
        (1-p)^2*R^2,2*(1-p)^2*R*(1-R),(1-p)^2*(1-R)^2),
           nrow=3,ncol=3,byrow=TRUE)


# beta[i] is the 'between marker genotype effects' for the ith marker genotype
beta <- c(2*D*(p*d - (D - p + 2*p*q)*h)/p^2,
       D*((1 - 2*p)*d + (2*D + (1 - 2*p)*(1 - 2*q))*h)/ (p*(1-p)),
       -2*D*((1 - p)*d + (D + (1-p)*(1-2*q))*h)/(1-p)^2 )
# omega[i,j] is the effect of the jth QTL genotype within the ith marker
# genotype.
omega <- matrix(NA,nrow=3,ncol=3)
omega[1,1] <- 2*(D - p*(1-q))*((D+p*q)*h -p*d)/p^2
omega[1,2] <-  -p*(2*D - p*(1 - 2*q))*d + (2*D^2 - 2*p*(1 - 2*q)*D +
             (1 - 2*q +2*q^2)*p^2)*h
omega[1,2] <- omega[1,2]/p^2
omega[1,3] <- 2*(D+p*q)*((D - p + p*q)*h - p*d)/p^2
omega[2,1] <- ((1-2*p)*D - 2*p*(1-p)*(1-q))*d + (2*D^2 + (1-2*p)*(1-2*q)*D + 2*p*q*(1-p)*(1-q))*h
omega[2,1] <-   -omega[2,1]/(p*(1-p))
omega[2,2]  <-  ((1-2*p)*D - p*(1-p)*(1-2*q))*d + (2*D^2 + (1-2*p)*(1-2*q)*D -p*(1-p)*(1-2*q+2*q^2))*h
omega[2,2] <-  -omega[2,2]/(p*(1-p))
omega[2,3] <- ((1-2*p)*D +2*p*q*(1-p))*d +(2*D^2 + (1-2*p)*(1-2*q)*D + 2*p*q*(1-p)*(1-q))*h
# apparent error (wrong sign) in algebra on p207
#omega[2,3] <- omega[2,3]/(p*(1-p))
omega[2,3] <-  - omega[2,3]/(p*(1-p))
omega[3,1] <- 2*(1+D-p-q+p*q)*((1-p)*d +(D-q+p*q)*h)/(1-p)^2
omega[3,2] <- (1-p)*(2*D + (1-p)*(1-2*q))*d + (2*D^2+2*(1-p)*(1-2*q)*D + (1-2*q+2*q^2)*(1-p)^2)*h
omega[3,2] <- omega[3,2]/(1-p)^2
omega[3,3] <- 2*(D-q+p*q)*((1-p)*d + (D + (1-p)*(1-q))*h)
omega[3,3] <- omega[3,3]/(1-p)^2

beta0 <- beta
omega0 <- omega
beta <- (1 - missclass.rate)*beta0
omega <-  beta0 - beta + omega0


I <- row(f2)
J <- col(f2)
su <- I < J
IJ <- cbind(c(I),c(J))


ijk <- expand.grid(I=1:3,J=1:3,K=1:3)
sjk <- ijk$J < ijk$K
ijkl <- expand.grid(I=1:3,J=1:3,K=1:3,L=1:3)
suu <- ijkl$I < ijkl$J & ijkl$K < ijkl$L                        

s1 <- n*sum(f1*beta^2) 
s1
s2 <- sum(f1*(1 + (n-1)*f1)*beta^2) 
s2
s3 <- 2*(n-1)*sum(f1[I[su]]*f1[J[su]]*beta[I[su]]*beta[J[su]]) 
s3
s4 <- sum(1/f1[I]*f2[IJ]*(1 + (n-1)*f2[IJ])*omega[IJ]^2) 
s4
#sum(1/f1[I]*f2[IJ]*(1 + (n-1)*f2[IJ])*omega[IJ]^2) 
s5 <- 0;for(i in 1:3)for(j in 1:3) for(k in 1:3) if (j < k) s5 <- s5 + 1/f1[i]*2*(n-1)*f2[i,j]*f2[i,k]*omega[i,j]*omega[i,k]
s5
s6 <- 0; for(i in 1:3)for(j in 1:3)s6 <- s6+f2[i,j]*(1+(n-1)*f2[i,j])*omega[i,j]^2
s6
s7 <- 0; for(i in 1:2)for(j in (i+1):3)for(k in 1:2)for(l in (k+1):3)
s7 <- s7+f2[i,j]*f2[k,l]*omega[i,j]*omega[k,l]
s7 <- 2*(n-1)*s7
s7

# pop1 Luo: 177.6, here 76.5
1/(G-1)*(s1-s2-s3+s4+s5-s6-s7) + sig2.error


ems.beta <- 1/(G - 1)*(
        n*sum(f1*beta^2) 
-       sum(f1*(1 + (n-1)*f1)*beta^2) 
-       2*(n-1)*sum(f1[I[su]]*f1[J[su]]*beta[I[su]]*beta[J[su]]) 
#+       sum(1/f1[I]*f2[I,J]*(1 + (n-1)*f2[I,J])*omega[I,J]^2) 
+       sum(1/f1[I]*f2[IJ]*(1 + (n-1)*f2[IJ])*omega[IJ]^2) 
+       sum(1/f1[ijk$I[sjk]]*2*(n-1)*f2[cbind(ijk$I[sjk],ijk$J[sjk])]*
            f2[cbind(ijk$I[sjk],ijk$K[sjk])]*
            omega[cbind(ijk$I[sjk],ijk$J[sjk])]*
            omega[cbind(ijk$I[sjk],ijk$K[sjk])])  
#-       sum(f2[I,J]*(1 + (n-1)*f2[I,J])*omega[I,J]^2) 
-       sum(f2[IJ]*(1 + (n-1)*f2[IJ])*omega[IJ]^2) 
-       2*(n-1)*sum(f2[cbind(ijkl$I[suu],ijkl$J[suu])] *
                    f2[cbind(ijkl$K[suu],ijkl$L[suu])] *
                    omega[cbind(ijkl$I[suu],ijkl$J[suu])] *
                    omega[cbind(ijkl$K[suu],ijkl$L[suu])])
 )  + sig2.error

s8 <- 0;for(i in 1:3)for(j in 1:3)s8 <- s8+f2[i,j]*omega[i,j]^2
s8 <- s8*n
s8
s9 <- s4
s9
s10 <- s5
s10
# pop1 Luo: 98.4, here 98.4
1/(n-G)*(s8 - s9 - s10) + sig2.error


ems.omega <-  1/(n-G)*(
    n*sum(f2[IJ]*omega[IJ]^2)  
-   sum(1/f1[I]*f2[IJ]*(1+(n-1)*f2[IJ])*omega[IJ]^2)  
-   sum(1/f1[ijk$I[sjk]]*2*(n-1)*f2[cbind(ijk$I[sjk],ijk$J[sjk])]*
    f2[cbind(ijk$I[sjk],ijk$K[sjk])]*omega[cbind(ijk$I[sjk],ijk$J[sjk])]*
        omega[cbind(ijk$I[sjk],ijk$K[sjk])])
)  + sig2.error


l1 <- method1.calc.ems(beta,omega,n,G,q,f1,f2,Vp,sig2.error,
                          missclass.rate=missclass.rate)
ems.beta <- l1$ems.beta
ems.omega <- l1$ems.omega






#Vq <- 0.5*d^2 + 0.25*h^2
#sig2.total <-  Vq + sig2.error

#browser()
# correct ems.beta, using total SS

if(FALSE){
if(correct.ems.beta)
  ems.beta <- 
    1/(G-1)*( (n-1)*Vp - (n-G)*ems.omega)
}

# expected variance components for Luo population 1
test.pop1 <- FALSE
if(test.pop1){
  calc.sig2(177.6,98.4,n,f1)
  calc.ems(2.53,98.4,n,f1)
}
# expected 
method <- "fixed"
if(method=="fixed"){
  l1 <- calc.sig2(ems.beta,ems.omega,n,f1)
  sig2.beta <- l1$sig2.beta
  sig2.omega <- l1$sig2.omega
}


debug <- FALSE
if(debug){
  ems.beta
  ems.omega
  calc.sig2(ems.beta,ems.omega,n,f1)
  calc.ems(sig2.beta,sig2.omega,n,f1)
  calc.ems(0,sig2.omega,n,f1)
}

#browser()                              
# in paper
ncp <- (ems.beta/ems.omega) * (nu1 * (nu2 - 1))/nu2  - nu1
# round off small negatives
if(D==0 || abs((ems.beta - ems.omega)/ems.omega) < 0.001) ncp <- 0
if(ncp < 0){
  cat("ld <- power: warning negative ncp set to zero: ncp =",ncp,"\n")
#  browser()
ncp <- 0
}

if(Q > 1 || Q < 0 || R > 1 || R < 0){
  power <- NA
}else{
  power <- 1 - pf(qf(1 - alpha, nu1,nu2), nu1,nu2,ncp)
}
if(print.it){
  cat("ems.beta =",round(ems.beta,2),
    ", ems.omega =",round(ems.omega,2),", power =",round(power,4),"\n")
  invisible(power)
}else power
}

oneway.bf.alpha <- 
function(n,group.sizes=c(0.25,0.5,0.25)*n,alpha=0.05){
  Fstat <- qf(1-alpha,2,n-3)
  SS.oneway.bf(group.sizes=group.sizes,Fstat=Fstat)
}

oneway.bf.alpha.required <-function(n,group.sizes=c(0.25,0.5,0.25)*n,Bf){
# function to calculate type one error rate required to achieve a given
# Bayes factor  
# n: population size
# group.sizes : expected number in each marker class
# Bf: Bayes factor(s)
# returns alpha, type 1 error rate
Fvals <-seq(0.1,30,by=0.1)  
Bfs <- SS.oneway.bf(group.sizes=group.sizes,Fstat=Fvals)
F.required <-approx(Bfs,Fvals,xout=Bf)$y
#browser()
alphas <- 1 - pf(F.required,2,n-3)
names(alphas) <- n
alphas
}

SS.oneway.bf <- function(group.sizes,Fstat)
# Bayes factor for oneway analysis of variance with improper prior
# a limit of the normal/inverse chi-squared conjugate prior
# From: Spiegelhalter, D.J.  and Smith, A.F.M. (1982)
# Bayes factors for linear and log-linear models with
# vague prior information JRSS B (1982) 44(3) 377-387
# Based on the idea that the Bayes factor for an
# imaginary minimal training sample should be approximately
# unity, to evaluate the constant.
# n.groups: number of groups
# group.sizes: size of each group
# Fstat: from analysis of variance  
# returns: Bayes factor B10 = p(y|M1)/p(y|M0)
{
if(any(group.sizes==0))group.sizes <- group.sizes[group.sizes > 0]
m <- length(group.sizes)
if(m<=1)stop("no models to compare")
n <- sum(group.sizes)
B01 <- sqrt(1/2*(m+1)*prod(group.sizes)/n)*(1 + (m-1)/(n-m)*Fstat)^(-n/2)
B10 <- 1/B01
B10
}  

ld.sim <- function(nsim, n, p, q, D, h2, Vp, phi, missclass.rate=0,
   sim.chunk=100, method=1, print.it=TRUE, data.only=FALSE){
# sim.chunk: number of nreps per call of ld.sim1()
# method 1 for calculations better if sim.chunk > 1
# nsim: number of simulations
# n: sample size
# p,q: marker and QTL allele frequencies
# D: disequilibrium
# h2: proportion of total variance explained by the QTL
# Vp: total or phenotypic variance
# phi: controls degree of dominance (phi=0 additive, phi=1 fully dominant)
# missclass.rate: miss-classification rate for marker genotyping
# method: Value 1 or 2. 
#   if method==1: simulate markers using p, then QTLs using Q,R.
#   if method==2: simulate markers and QTL directly from table of joint probs (f <- ij)
# print.it: if TRUE, print results
if(h2<0 || h2 > 1)stop("heritability not in [0,1]")
Vq <- h2*Vp
#d <- sqrt(Vq/(0.5 + 0.25*phi^2))  # only for q=0.5
d <- sqrt(Vq/( (2*q - 2*q^2) + (8*q^3-12*q^2 + 4*q)*phi + (2*q -6*q^2 +8*q^3 - 4*q^4)*phi^2))
h <- phi*d
sig2.error <- Vp - Vq
res <- list()
ncalls <- 1+trunc(nsim/sim.chunk)
for(ii in 1:ncalls){
  res[[ii]] <- ld.sim1(n, p, q, D, d, h, sig2.error, missclass.rate,
                    nreps=sim.chunk, method=method, data.only=data.only)
  cat(".")
}
cat("\n")
#browser()
if(data.only){ # return a dataframe with marker, trait and replicate info.
#browser()
nys <- sapply(res,ncol) - 1
res1 <- lapply(as.list(seq(along=res)),function(ii,res,nys){
  res1 <- res[[ii]]; ny1 <- nys[ii]
  marker <- rep(res1[,1],ny1)
  y <- c(res1[,-1])
  cbind(marker=c(marker),y=y,repl=rep(1:ny1,rep(length(res1[,1]),ny1)))
}, res=res,nys=nys)
#str(res1)

ns <- sapply(res1,nrow)

res2 <- data.frame(marker=c(sapply(res1,function(m)m[,1])),
      y=c(sapply(res1,function(m)m[,2])))
repl1 <- c(sapply(res1,function(m)m[,3]))
repl2 <- rep(1:length(res1),sapply(res1,nrow))
res2$replicate <- (repl2-1)*max(repl1) + repl1
#str(res2)

#browser()              
res2
}else{
  matrix(aperm(array(sapply(res,c),dim=c(sim.chunk,4,ncalls)),c(1,3,2)),ncol=4)
}
}

ld.sim1 <- function(n, p, q, D, d, h, sig2.error, missclass.rate=0, nreps=1,
                 method=2, print.it=TRUE, data.only=FALSE){
# simulate marker and QTL values for a population and test for linkage disequilibrium
# do a single simulation or nreps replicate simulations with the same set of marker genotypes
# n : sample size
# p : marker allele frequency
# q : QTL allele frequency
# D : linkage disequlibrium
# d,h : QTL parameters
# sig2.error : error variance when QTL genotype known and modelled
# missclass.rate :  proportion assumed to be missclassified
# nreps: do nreps  replicate simulations at once with same set of marker genotypes
#        using multivariate anova
# method: Value 1 or 2. 
#   if method==1: simulate markers using p, then QTLs using Q,R. Marker values common
#                 marker values common to all runs, => use manova
#   if method==2: simulate markers and QTL directly from table of joint probs (f <- ij)

Q  <-  q + D/p
R  <-  q - D/(1-p)

n0 <- trunc(n*(1 - missclass.rate))
n1 <- n - n0

# method 1: simulate markers using p, then QTLs using Q,R.
# simulate random marker genotypes, for 2 chromosomes

if(method==1){
mval1 <- rbinom(n,size=1,prob=p)
mval2 <- rbinom(n,size=1,prob=p)
# simulate random QTL genotypes, conditional on marker values
# n1 missclassified values sampled at random
r1 <- rbinom(n1*nreps,size=1,prob=q)
r2 <- rbinom(n1*nreps,size=1,prob=q)
qval11.M <- matrix(rbinom(n0*nreps,size=1,prob=Q),ncol=nreps)
qval12.M <- matrix(rbinom(n0*nreps,size=1,prob=R),ncol=nreps)
qval21.M <- matrix(rbinom(n0*nreps,size=1,prob=Q),ncol=nreps)
qval22.M <- matrix(rbinom(n0*nreps,size=1,prob=R),ncol=nreps)
qval1 <- ifelse(rep(mval1[1:n0]==1,nreps),qval11.M,qval12.M)
qval2 <- ifelse(rep(mval2[1:n0]==1,nreps),qval21.M,qval22.M)

#browser()
yq0 <- matrix(ifelse(qval1==1,ifelse(qval2==1,d,h),ifelse(qval2==1,h,-d)),ncol=nreps)
yq1 <- matrix(ifelse(r1==1,ifelse(r2==1,d,h),ifelse(r2==1,h,-d)),ncol=nreps)
yq <- rbind(yq0,yq1)

MM <- factor(ifelse(mval1==1,ifelse(mval2==1,"MM","Mm"),ifelse(mval2==1,
                                   "Mm","mm")))

} 
# end if(method==1)

if(method==2){
# method 2: simulate markers and QTL directly from table of joint probs (f <- ij)
f2 <- matrix(c(p^2*Q^2, 2*p^2*Q*(1-Q), p^2*(1-Q)^2,
     2*p*(1-p)*Q*R, 2*p*(1-p)*(Q+R-2*Q*R), 2*p*(1-p)*(1-Q)*(1-R),
     (1-p)^2*R^2, 2*(1-p)^2*R*(1-R), (1-p)^2*(1-R)^2), ncol=3,byrow=TRUE)

f2.0 <- f2
f1.q <- c(q^2,2*q*(1-q),(1-q)^2)
f2.1 <- matrix(rep(f1.q,3),ncol=3,byrow=TRUE)

#browser()
mq0 <- sample(1:9,size=n0*nreps,replace=TRUE, prob=c(f2))
mq1 <- sample(1:9,size=n1*nreps,replace=TRUE, prob=c(f2.1))
mq <- rbind(matrix(mq0,ncol=nreps),matrix(mq1,ncol=nreps))

yq <- c(d,h,-d)[col(f2)[mq]]
yq <- matrix(yq,ncol=nreps)
MM <- factor(c("MM","Mm","mm")[row(f2)[mq]])
} # end method 2

#browser()
mu <- 0
err <- rnorm(n*nreps)*sqrt(sig2.error)
err <- matrix(err,ncol=nreps)

y <- mu + yq + err

if(data.only){
#browser()
res <- cbind(MM,y)
res
}else{
 if(method==1){ # marker values common -> manova
 fit1 <- aov(y ~ MM)
 summ1 <- summary(fit1)
 MS <- sapply(summ1,function(u)u$"Mean Sq")
 summ1
 res <- cbind(MS.beta=MS[1,],MS.within=MS[2,],F.value=MS[1,]/MS[2,],
     P.value=sapply(summ1,function(u)u$"Pr(>F)"[1]))
 res
 }else{ # marker values different -> sep anova
  res <- matrix(nrow=nreps,ncol=4)
  repl <- col(y)
  for(ii in 1:nrow(res)){
    fit1 <- aov(y[,ii] ~ MM[repl==ii])
    summ1 <- summary(fit1)
#browser()
    MS <- summ1[[1]]$"Mean Sq"
    res[ii,] <- c(MS[1],MS[2],MS[1]/MS[2],summ1[[1]]$"Pr(>F)"[1])

#browser()
  }
  dimnames(res) <- list(NULL,c("MS.beta","MS.within","F.value","P.value"))
  res
 }
}

}




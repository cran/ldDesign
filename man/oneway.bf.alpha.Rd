\name{oneway.bf.alpha}
\alias{oneway.bf.alpha}
\alias{oneway.bf.alpha.required}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Correspondence between significance levels and Bayes factors for
  effects of marker genotype classes.}
\description{
Functions to calculate the correspondence between significance levels alpha 
  and the Bayes factor, for association between a bi-allelic marker and
  QTL, for given sample sizes and marker genotype frequencies for
  bi-allelic marker.
}
\usage{
  oneway.bf.alpha(n, group.sizes = c(0.25, 0.5, 0.25) * n, alpha = 0.05)
  oneway.bf.alpha.required(n, group.sizes = c(0.25, 0.5, 0.25) * n, Bf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size, i.e. number of individuals genotyped
    and phenotyped for the trait}
  \item{group.sizes}{Number in each of the 3 possible marker genotype
    classes MM, Mm, mm}
  \item{alpha}{Significance level, i.e. threshold for `detection'}
  \item{Bf}{Bayes factor, used as threshold for detection} }

\details{These functions implement the correspondence between the 
  significance levels and Bayes factors used in Ball (2003) to design
  experiments for detecting linkage disequilibrium with a given power 
  to achieve a given Bayes factor. The function \code{\link{SS.oneway.bf}} is used to calculate
  the Bayes factor corresponding to a given F statistic (Spiegelhalter
  and Smith 1982). This is combined with a call to \code{qf}, for \code{oneway.bf.alpha}
  or calls to \code{pf} and interpolation for
  \code{oneway.bf.alpha.required}, to calculate the Bayes factor
  corresponding to a given alpha or alpha values for a given Bayes factor.
 }
}
\value{
\code{oneway.bf.alpha} returns the Bayes factor
corresponding to a given significance level (alpha).
\code{oneway.bf.alpha.required} returns the significance level (alpha)
corresponding to a given Bayes factor.
}
\references{ 
  Ball, R.D. 2003
  Experimental designs for reliable detection of linkage disequilibrium
  in unstructured random population association studies.

  Spiegelhalter, D. and A.F.M. Smith 1982
  Bayes factors for linear and log-linear models with vague prior
  information J. Royal Statist Soc. B 44: 377--387.
}
\author{Rod Ball \email{rod.ball@forestresearch.co.nz} 
\url{www.forestresearch.co.nz}}
%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{\code{\link{SS.oneway.bf}},\code{\link{luo.ld.power}} }
\examples{
# calculations for Table 4 in the manuscript
data(luo.ld.populations)
Bs <- numeric(nrow(luo.ld.populations))
n.Bf20s <- numeric(nrow(luo.ld.populations))
ns <- c(seq(200,400,by=25),,450,seq(500,4000,by=100))
powers <- numeric(length(ns))
alphas <- numeric(length(ns))
P.Bf20s <- numeric(length(ns))
for(ii in 1:nrow(luo.ld.populations)){
   cat("ii=",ii,"\n")
   powers[ii] <- luo.ld.power(n=luo.ld.populations[ii,"n"],
                              p=luo.ld.populations[ii,"p"],
                              q=luo.ld.populations[ii,"q"],
                              D=luo.ld.populations[ii,"D"],
                              h2=luo.ld.populations[ii,"h2"],
                              phi=luo.ld.populations[ii,"phi"],
                              Vp=100,
                              alpha=0.05)

   p1 <- luo.ld.populations[ii,"p"]
   Bs[ii] <- oneway.bf.alpha(n=luo.ld.populations[ii,"n"],
                 group.sizes=c(p1^2,2*p1*(1-p1),(1-p1)^2)*
                                 luo.ld.populations[ii,"n"])
   for(jj in seq(along=ns)){
     alphas[jj] <- oneway.bf.alpha.required(ns[jj],
                       group.sizes=c(p1^2,2*p1*(1-p1),(1-p1)^2)*ns[jj],Bf=20)
     P.Bf20s[jj] <- luo.ld.power(n=ns[jj],
                                 p=luo.ld.populations[ii,"p"],
                                 luo.ld.populations[ii,"q"],
                                 D=luo.ld.populations[ii,"D"],
                                 h2=luo.ld.populations[ii,"h2"],
                                 phi=luo.ld.populations[ii,"phi"],
                                 Vp=100,
                                 alpha=alphas[jj],
                                 print.it=FALSE)
   }
   n.Bf20s[ii] <- approx(P.Bf20s,ns,xout=0.9)$y
   cat("n =",n.Bf20s[ii],"\n")
}
cbind(luo.ld.populations,powers,n.Bf20s)

}



\keyword{htest}
\keyword{models}

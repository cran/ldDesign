\name{luo.ld.power}
\alias{luo.ld.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classical deterministic power calculation for association studies to detect
  linkage disequilibrium}
\description{Classical deterministic power calculation for power to detect linkage
disequilibrium between a bi-allelic QTL and a bi-allelic marker, at a
given significance level in a population level association study.  
}
\usage{
luo.ld.power(n, p, q, D, h2, phi, Vp = 100, alpha, print.it = TRUE, 
             missclass.rate = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The sample size, i.e. number of individuals
    genotyped and tested for the trait of interest }
  \item{p}{Bi-allelic marker allele frequency}
  \item{q}{Bi-allelic QTL allele frequency}
  \item{D}{Linkage disequilibrium coefficient}
  \item{h2}{QTL `heritability', i.e. proportion of total or phenotypic
    variance explained by the QTL}
  \item{phi}{Dominance ratio: \code{phi} = 0 denotes purely additive, 
\code{phi} = 1 denotes purely dominant allele effects}
  \item{Vp}{Total or phenotypic variance: and arbitrary value may be used}
  \item{alpha}{Significance level for hypothesis tests}
  \item{print.it}{If TRUE print summary of results}
  \item{missclass.rate}{Proportion of marker values which are
    missclassified, i.e. incorrect}
}
\details{
This is based on the 'fixed model' power calculation from Luo (1998,
Heredity 80, 198--208), with corrections described in Ball (2003).  This
function, combined with \code{\link{oneway.bf.alpha}}, 
\code{\link{oneway.bf.alpha.required}}, is used in
Ball (2003) to design experiments for detecting linkage disequilibrium
with a given power to achieve a given Bayes factor.
}
\value{
  Returns the power, or probability of detecting an effect, with the
  given parameters, at the given significance level.
}

\references{
  Ball, R.D. 2003
  Experimental designs for reliable detection of linkage disequilibrium
  in unstructured random population association studies.

  Luo, Z.W. 1988
  Detecting linkage disequilibrium between a polymorphic marker locus and 
  a trait locus in natural populations. Heredity 80, 198--208}

\author{Rod Ball \email{rod.ball@forestresearch.co.nz} 
\url{www.forestresearch.co.nz}}
%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{\code{\link{ld.sim}},\code{\link{oneway.bf.alpha}},\code{\link{oneway.bf.alpha.required}},
\code{\link{SS.oneway.bf}}}
\examples{
data(luo.ld.populations)
options(digits=3)
powers <- numeric(nrow(luo.ld.populations))
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
}
cbind(luo.ld.populations,power=powers)
}
\keyword{design}% at least one, from doc/KEYWORDS


\name{ld.design}
\alias{ld.design}
\alias{ld.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for design of experiments to detect linkage disequilibrium}
\description{
  Find the sample size required to detect linkage disequilibrium with
  a given Bayes factor, with a given power, or find the power of
  experimental designs to detect linkage equilibrium with a given
  Bayes factor.
}
\usage{
ld.design(p, q, D, h2, phi, Bf, power, nmin = 50, nmax = 1e+05, ninterp = 50, 
          missclass.rate = 0, print.it = FALSE)
ld.power(n, p, q, D, h2, phi, Bf, missclass.rate = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{\code{ld.power:} vector of sample sizes}
  \item{p}{Bi-allelic marker allele frequency}
  \item{q}{Bi-allelic QTL allele frequency}
  \item{D}{Linkage disequilibrium coefficient}
  \item{h2}{QTL `heritability', i.e. proportion of total or phenotypic
    variance explained by the QTL}
  \item{phi}{Dominance ratio: \code{phi} = 0 denotes purely additive, 
    \code{phi} = 1 denotes purely dominant allele effects}
  \item{Bf}{Bayes factor}
  \item{power}{\code{ld.design:} Power, or probability of detecting an effect with Bayes
    factor greater than Bf}
  \item{nmin}{\code{ld.design:} Lower bound for sample size}
  \item{nmax}{\code{ld.design:} Upper bound for sample size}
  \item{ninterp}{\code{ld.design:} Number of sample sizes to try}
  \item{missclass.rate}{Proportion of marker values which are
    missclassified, i.e. incorrect (to allow for genotyping errors)}
  \item{print.it}{If TRUE print results for sample sizes tried}
}
\details{
These functions implement the method described in Ball (2005) for
obtaining the power of designs for detecting linkage disequilibrium with
a given Bayes factor. The F values, (and hence significance levels) 
corresponding to the given Bayes factors, sample sizes, and marker genotype 
frequecies, are calculated using the method of Spiegelhalter and Smith (1982) 
(R functions \code{\link{oneway.bf.alpha.required}},
\code{\link{SS.oneway.bf}}). The power is obtained using a corrected
version of the classical deterministic power calculation from 
Luo (1988) (R function \code{\link{luo.ld.power}}).
}

\value{
  For \code{ld.power}, a matrix with columns: 
  \item{n}{Sample sizes}
  \item{power}{Power of the design with the given sample sizes}
  Additionally the return value has attributes indicating the
  linkage disequilibrium parameters used.
  For \code{ld.design} the sample size is returned.
}
\references{
  Ball, R.D. 2005: Experimental designs for reliable detection of linkage
  disequilibrium in unstructured random population association studies.
  Genetics 170: 859--873.

  Ball, R.D. 2007:  Statistical analysis and experimental design.
  Chapter 8, pp133--196 In: Association Mapping in Plants, N.C. Oraguzie,
  E.H.A. Rikkerink, S.E. Gardiner, and H.N. DeSilva (Editors), Springer,
  New York.

  Luo, Z.W. 1988:
  Detecting linkage disequilibrium between a polymorphic marker locus and 
  a trait locus in natural populations. Heredity 80, 198--208

  Spiegelhalter, D. and A.F.M. Smith 1982:
  Bayes factors for linear and log-linear models with vague prior
  information J. Royal Statist Soc. B 44: 377--387.
}
\author{Rod Ball \email{rod.ball@scionresearch.com}
\url{www.scionresearch.com}}
%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{ \code{\link{luo.ld.power}}, \code{\link{ld.sim}}, \code{\link{oneway.bf.alpha}},
\code{\link{oneway.bf.alpha.required}}, \code{\link{SS.oneway.bf}}}
\examples{
ld.power(n=seq(100,1000,by=100),p=0.5,q=0.5,D=0.1,h2=0.1,phi=0,Bf=20)
ld.design(p=0.5,q=0.5,D=0.1,h2=0.1,phi=0,Bf=20,power=0.9,print.it=TRUE,nmin=600,nmax=4000)
ld.design(p=0.5,q=0.5,D=0.1,h2=0.1,phi=0,Bf=20,power=0.9,print.it=FALSE,nmin=1700,nmax=1900)
}
\keyword{design }% at least one, from doc/KEYWORDS


\name{calc.B.ABF}
\alias{calc.B.ABF}
\alias{calc.Zc.ABF}
\alias{calc.Zalpha.ABF}
\alias{calc.Balpha.ABF}
\alias{calc.alphaB.ABF}

\title{Calculate Z-values or alpha-values from approximate Bayes factors 
and vice versa.}

\description{These functions find the critical values of \code{Z}
and \code{alpha} corresponding to a given approximate Bayes factor, or
vice versa for given sample size (\code{n}) and prior equivalent sample
size (\code{a}), where {Z} is a statistic with sampling distribution
N(0, 1/n) for sample size \code{n}.  The critical value Zc is used by
\link{cc.power} to calculate the power, for an experiment to obtain a
given Bayes factor, and the corresponding alpha value alphac may be used
by the Genetic power calculator function \link{gpc.power}.}

\usage{
calc.B.ABF(Z,n,a=1)
calc.Zc.ABF(B,n,a=1)
calc.Zalpha.ABF(alpha,n)
calc.Balpha.ABF(alpha, n, a)
calc.alphaB.ABF(B,n,a=1,alpha.start=1e-5,reduction.factor=1.5,niter=20, verbose=FALSE, show.progress=FALSE)
}

\arguments{
  \item{B}{Bayes factor}
  \item{Z}{a sample statistic with sampling distribution N(0,1/n) estimating a 
	multiple of a quantity of interest}
  \item{n}{sample size}
  \item{a}{equivalent sample size for information in prior}
  \item{alpha}{significance level}
  \item{alpha.start}{initial guess (upper bound) for significance level}
  \item{reduction.factor}{factor to decrease \code{alpha.start} by each iteration
	 in the search for an alpha value corresponding to B}
  \item{niter}{number of iterations in search}
  \item{show.progress}{if TRUE print dots to show progress}
  \item{verbose}{if TRUE print array of alpha and B values calculated}
}

\details{
These functions calculate the approximate Bayes factor and threshold
calculations used in the power calculations for case-control studies
(\link{cc.power}, \link{cc.design}), but can be used more widely. For given 
\code{Z, n, a}, the approximate Bayes factor is given by the formula:

(1)	  B = sqrt(a/(n+a))*exp(n^2*Z^2/(2*(n+a)))

where \code{n} is the sample size and the sampling distribution of
\code{Z} is N(0,1/n), and the prior for \code{Z} is N(0,1/a), (so that
\code{a} is the equivalent sample size for information in the prior). In
applications a subjective prior for the quantity of interest can be
given and the corresponding value of a determined.}

\value{
For \code{calc.B.ABF} the Bayes factor from (1); for \code{calc.Zc.ABF} the value
of \code{Z} corresponding to \code{B} in (1); 
}

\references{
  Ball, R.D. 2005: Experimental designs for reliable detection of linkage
  disequilibrium in unstructured random population association studies.
  Genetics 170: 859--873.

  Ball, R.D. 2007 (Jan.): ``Statistical analysis and experimental
  design'', Chapter 8, In: Association mapping in plants. N.C. Oraguzie
  et al. editors, Springer Verlag, ISBN 0387358447. (69pp) (Approximate
  Bayes factors for S-TDT test: pp166-167)

  Wakefield, J. 2007 (Jul.): A Bayesian measure of the probability of
  false discovery in genetic epidemiology studies.  Am. J. Hum. Genet. 81:
  208--227.

  Ball, R.D. 2011: Experimental designs for robust detection of effects in genome-wide
  case-control studies (submitted).
}

\author{Rod Ball \email{rod.ball@scionresearch.com}
\url{www.scionresearch.com}}

\examples{
calc.Zc.ABF(B=1e6,n=10000,a=1) 
calc.B.ABF(Z=0.0607,n=10000,a=1)
calc.Balpha.ABF(alpha=0.01,n=100,a=1)
calc.alphaB.ABF(B=1e6,n=10000,a=1,alpha.start=1e-7)
calc.Zalpha.ABF(alpha=1.28e-9,n=10000)
}

\keyword{design, htest}

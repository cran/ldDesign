\name{SS.oneway.bf}
\alias{SS.oneway.bf}
\title{Bayes factors for one-way analysis of variance models.
}
\description{
Function to calculate the Bayes factor for a one-way analysis of
  variance layout with vague or improper priors.
}
\usage{
SS.oneway.bf(group.sizes, Fstat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group.sizes}{Sizes of groups in the one-way layout}
  \item{Fstat}{F statistic obtained}
}
\details{
  The function the Bayes factor corresponding to a given F statistic
  in a one-way analysis of variance model is calculated using the method 
  of Spiegelhalter and Smith 1982. With improper priors the marginal
  probabilities of the data under each of the models (corresponding to
  the NULL and alternative hypotheses) is indeterminate. This is
  resolved by updating each prior with a small imaginary training
  sample, which is equivalent to normalising the Bayes factor to be 1
  for the small training sample. Spiegelhalter and Smith obtain a
  formula for the Bayes factor in terms of the classical F value.
  
}
\value{
  Returns the Bayes factor corresponding to the given design and
  observed value of F statistic.
}
\references{
Spiegelhalter, D. and A.F.M. Smith 1982
  Bayes factors for linear and log-linear models with vague prior
  information J. Royal Statist Soc. B 44: 377--387.
}
\author{Rod Ball \email{rod.ball@scionresearch.com} 
\url{www.scionresearch.com}}
%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{\code{\link{oneway.bf.alpha}},\code{\link{oneway.bf.alpha.required}}}
\examples{
# Bayes factors corresponding to P-values 0.05,0.01,0.001,0.0001 for n=200
SS.oneway.bf(group.sizes=c(50,100,50),Fstat=qf(0.95,2,197))
SS.oneway.bf(group.sizes=c(50,100,50),Fstat=qf(0.99,2,197))
SS.oneway.bf(group.sizes=c(50,100,50),Fstat=qf(0.999,2,197))
SS.oneway.bf(group.sizes=c(50,100,50),Fstat=qf(0.9999,2,197))
}

\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line

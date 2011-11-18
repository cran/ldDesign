\name{gpc.power}
\alias{gpc.power}
\title{Power calculations for case-control studies using the Genetic
Power calculator}
\description{
  Find the power of tests and the sample sizes required for case-control for given
frequentist significance level (alphac). Uses the Genetic Power Calculator 
web site (Purcell 2008). Can be used to check the Bayesian calculations from
\link{cc.design}, \link{cc.power},  if called with the appropriate alpha level 
(alphac attribute from the output of \link{cc.power}.  }
\usage{
gpc.power(Dprime, p, q, prevalence, R1, R2, n.cases, n.controls, alpha, power,
	show.html=TRUE)
}
\arguments{
  \item{Dprime}{Dprime, i.e. linkage disequilibrium as a proportion of the maximum
	 (or minimum, if negative)}
  \item{p}{Bi-allelic marker allele frequency}
  \item{q}{Bi-allelic QTL allele frequency (for risk allele)}
  \item{prevalence}{disease prevalence, i.e. the probability of being a case in
        the population}
  \item{R1}{Relative risk for the genotype with one copy of the risk allele, compared
	with the genotype with no copies of the risk allele}
  \item{R2}{Relative risk for the genotype with two copies of the risk allele}
  \item{n.cases}{number of cases}
  \item{n.controls}{number of controls}
  \item{alpha}{significance level or threshold}
  \item{power}{Power, or probability of detecting an effect with the given relative risks,
	with significance level alpha (i.e. p < alpha)}
  \item{show.html}{if TRUE, try to show the html web page returned by the query}
}
\details{
This function generates a query to the Genetic Power Calculator website, and extracts the
power and sample sizes corresponding to the input values of \code{alpha} and
\code{power}. Requires wget, sed, and files in the ./bin subdirectory.
}
\value{A table showing the values of alphac, power, and number of cases. The power is
the power to obtain the specified significance level alphac, with the specified numbers
of cases and controls. The number of cases required is the number of cases required
to obtain the specified power with the specified significance level}
\references{
Purcell, S., Cherny, S.S., and Sham, P.C. 2003:
  Genetic Power Calculator: design of linkage and association genetic
  mapping studies of complex traits.  Bioinformatics 19: 149--150.
\url{http://pngu.mgh.harvard.edu/~purcell/gpc/cc2.html}

Purcell, S. 2008: Genetic Power Calculator: 
   \url{http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi}

}
\author{Rod Ball \email{rod.ball@scionresearch.com}
\url{www.scionresearch.com}}
\examples{
# Not run
power.res <- cc.power(B=1e6, R=1.2, Dprime=0.5, p=0.5, q=0.5, prevalence=0.13,
         n.cases=147546/2, n.controls=147546/2, model="dominant", a=1,
         verbose=TRUE)
n <- 147546
alphac <- attr(power.res,"alphac")
## use R1 = R2 for the dominant model, R1=1.0 for the recessive model
#gpc.res <- gpc.power(Dprime=0.5,p=0.5,q=0.5,prevalence=0.13,R1=1.2,R2=1.2,n.cases=n/2,
#                     n.controls=n/2,alpha=alphac,power=0.8)
## 
## gpc.res
##                 alphac   power ncases(pwr=80%)
## dominant      8.31e-10 0.80030           73752
## recessive     8.31e-10 0.02266          210022
## general (2df) 8.31e-10 0.80260           73580
## allelic       8.31e-10 0.74160           78035
}



  





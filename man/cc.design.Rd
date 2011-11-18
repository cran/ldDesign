\name{cc.design}
\alias{cc.design}
\alias{cc.power}
\title{Functions for design of experiments to detect linkage disequilibrium in genome-wide case-control studies}
\description{
  Find the sample sizes (number of cases and number of controls) required to detect linkage disequilibrium with a given Bayes factor, with a given power, or find the power of experimental designs to detect linkage equilibrium with a given Bayes factor.  }
\usage{
cc.design(B, OR, D, p, q, power, baseline.risk, Dprime=NULL, R=NULL, 
	  prevalence=NULL, n.cases, n.controls, model=c("additive",
	  "dominant","recessive","general"), a=1, sigma2.eta=NULL, verbose=FALSE, 
	  amalgamate.cells=FALSE, pmin=0.1, pmax=0.99, ninterp=20, print.power.curve=TRUE)
cc.power(B, OR, D, p, q, baseline.risk, Dprime=NULL, R=NULL, prevalence=NULL, 
	 n.cases, n.controls, model=c("additive","dominant","recessive","general"), 
	 a=1, sigma2.eta=NULL, verbose=FALSE, amalgamate.cells=FALSE,show.attributes=FALSE)
}
\arguments{
  \item{B}{Bayes factor}
  \item{OR}{Odds ratio}
  \item{D}{Linkage disequilibrium coefficient}
  \item{p}{Bi-allelic marker allele frequency}
  \item{q}{Bi-allelic QTL allele frequency (for risk allele)}
  \item{power}{\code{cc.design:} Power, or probability of detecting an effect with Bayes
    factor greater than B}
  \item{baseline.risk}{Baseline risk, i.e. the probability of being a case for a
	genotype with no risk alleles}
  \item{Dprime}{D': i.e. linkage disequilibrium as proportion of the maximum 
	(or minimum, if negative); need to give \code{D} _or_ \code{Dprime}}
  \item{R}{Relative risk: in the additive model: relative risk per copy of the
	risk allele, in dominant or recessive models relative risk for the high
	risk genotype(s) compared with the low risk genotype(s) in the 2x2
	contingency table; need to give \code{R} _or_ \code{OR}.}
  \item{prevalence}{disease prevalence, i.e. the probability of being a case in
	the population; need to give \code{prevalence} _or_ \code{baseline.risk}}
  \item{n.cases}{number of cases (can be vectorised), for \code{cc.design} only used for 
	determining the proportion of cases and controls}
  \item{n.controls}{number of controls (can be vectorised) for \code{cc.design} only used for 
        determining the proportion of cases and controls}
  \item{model}{genetic model assumed for the locus, one of: "additive", "dominant",
	"recessive" or "general"}
  \item{a}{equivalent sample size for information in prior, (alternative to specifying 
	sigma2.eta); a=1 is a conservative `default' value}
  \item{sigma2.eta}{prior variance for the log odds ratio}
  \item{verbose}{if TRUE display intermediate information}
  \item{amalgamate.cells}{if TRUE, amalgamate the genotypes with 1 or 2
        copies of the risk allele (additive model only); this may be advantageous when 
	the risk allele frequency (q) is low}
  \item{pmin}{\code{cc.design:} Lower bound for power --- power curve calculated from pmin to pmax}
  \item{pmax}{\code{cc.design:} Upper bound for power --- power curve calculated from pmin to pmax}
  \item{ninterp}{\code{cc.design:} Number of sample sizes to calculate for interpolation}
  \item{print.power.curve}{\code{cc.design:} If TRUE print power curve for power from pmin to pmax}
  \item{show.attributes}{\code{cc.power:} If TRUE print attributes of the returned object, otherwise
               just print the power.}
}
\details{
These functions implement the method described in Ball (2011) for
obtaining the power of designs for detecting linkage disequilibrium with
a given Bayes factor in genome-wide case-control studies. The method
uses an asymptotic approximation to the Bayes factor (Ball 2007, 2011;
Wakefield 2007) together with an estimator for the log-odds ratio and its
standard error.  This extends the method from Ball (2005) to case-control
studies used for studying human diseases.
}
\value{
  For \code{cc.power}, the power for the given sample size(s) (\code{n.cases}, 
  \code{n.controls}). The returned object also includes attributes including
  Bayes factor, model, sample sizes, odds ratios, relative risks, non-centrality parameter (\code{ncp}), 
  expected marker genotype frequencies for cases and controls (ps), and equivalent frequentist
  threshold (\code{alphac}).
  For \code{cc.design}, the sample sizes (\code{n}, \code{n.cases}, \code{n.controls}) 
  required to obtain the specified power.
}

\references{
  Ball, R.D. 2005: Experimental designs for reliable detection of linkage
  disequilibrium in unstructured random population association studies.
  Genetics 170: 859--873.

  Ball, R.D. 2007: ``Statistical analysis and experimental
  design'', Chapter 8, In: Association mapping in plants. N.C. Oraguzie
  et al. editors, Springer Verlag, ISBN 0387358447. (69pp)
 
  Ball, R.D. 2011: Experimental designs for robust detection of effects in genome-wide
  case-control studies (submitted).
}

\author{Rod Ball \email{rod.ball@scionresearch.com}
\url{www.scionresearch.com}}
\examples{
cc.power(B=1e6, OR=1.6, D=0.1, p=0.3, q=0.2, baseline.risk=0.1, 
         n.cases=1000*seq(2,12,by=2), n.controls=1000*seq(2,12,by=2),
	 model="additive", a=1)
cc.design(B=1e6, OR=2.0, D=0.1, p=0.3, q=0.2, power=0.9, baseline.risk=0.1, 
	  n.cases=2000, n.controls=3000, model="additive", a=1, 
	  pmin=0.1, pmax=0.99, ninterp=20, print.power.curve=TRUE)
cc.design(B=1e6, OR=c(2.0,2.0), D=0.1, p=0.3, q=0.2, power=0.9, baseline.risk=0.1, 
	  n.cases=2000, n.controls=3000, model="general", a=1, 
	  pmin=0.1, pmax=0.99, ninterp=20, print.power.curve=TRUE)
cc.design(B=1e6, OR=3.0, Dprime=0.5, p=0.03, q=0.02, power=0.9, baseline.risk=0.1,
          n.cases=5000, n.controls=5000, model="recessive",
	  pmin=0.1, pmax=0.9, ninterp=20, print.power.curve=TRUE)
cc.design(B=1e6,OR=c(1.2,1.2),D=0.25,p=0.5,q=0.5,prevalence=0.13,power=0.8,
                             n.cases=5000,n.controls=5000,model="general",a=1)
}

\keyword{design }% at least one, from doc/KEYWORDS

\name{ld.sim}
\alias{ld.sim}
\alias{ld.sim1}
\title{Functions to simulate populations with a bi-allelic marker and QTL in linkage 
disequilibrium and test for association.}
\description{
For a bi-allelic marker and QTL, with given allele frequencies, linkage disequilibrium, and QTL 
heritability, multiple replicate populations with marker, QTL, and trait values are simulated 
and tested for a marker-trait association. Results can be used to estimate the power of an experimental 
design for detecting linkage disequilibrium.
}
\usage{
ld.sim(nsim, n, p, q, D, h2, Vp, phi, missclass.rate = 0, sim.chunk = 100, 
       method = 1, print.it = TRUE, data.only=FALSE)
ld.sim1(n, p, q, D, d, h, sig2.error, missclass.rate = 0, nreps = 1, 
        method = 2, print.it = TRUE, data.only=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{Number of replicate simulations to do}
  \item{n}{The sample size, i.e. number of individuals
    genotyped and tested for the trait of interest }
  \item{p}{Bi-allelic marker allele frequency}
  \item{q}{Bi-allelic QTL allele frequency}
  \item{D}{Linkage disequilibrium coefficient}
  \item{h2}{QTL `heritability', i.e. proportion of total or phenotypic
    variance explained by the QTL}
  \item{Vp}{\code{ld.sim:} Total or phenotypic variance: an arbitrary value may be used}
  \item{phi}{\code{ld.sim:} Dominance ratio: \code{phi = 0} denotes purely additive
    gene action, \code{phi = 1} denotes completely dominant gene action}
  \item{d}{\code{ld.sim1:} Expected value for trait when QTL genotype is QQ,qq respectively is d,-d}
  \item{h}{\code{ld.sim1:} Expected value for trait when QTL genotype is Qq is h}
  \item{sig2.error}{\code{ld.sim1:} Error variance when QTL genotype known and modelled}
  \item{missclass.rate}{Proportion of marker values which are
    missclassified, i.e. incorrect}
  \item{sim.chunk}{\code{ld.sim:} Number of replicates to do in a
    `chunk' in each call to \code{ld.sim1}}
  \item{nreps}{\code{ld.sim1:} Number of replicate simulations to do for
    the given set of marker genotypes}
  \item{method}{If \code{method=1} simulate random QTL genotypes
    conditional on marker values in \code{ld.sim1}; if \code{method=2} 
    simulate markers and QTL directly from table of joint probabilities. 
    With \code{method=1}, a common set of marker values are used for each of 
    the \code{nreps} replicates per call to \code{ld.sim1}, enabling MANOVA to be used.}
  \item{print.it}{if \code{TRUE}, print results}
  \item{data.only}{if \code{TRUE}, just return the simulated trait and 
    marker genotype data}
}
\details{
  Marker, QTL, and trait values are simulated according to the genetic
  model with normal errors. In \code{ld.sim}, QTL parameters \code{d,h} are determined from the
  parameters \code{h2,q,phi}, and \code{Vq}, and the main simulation done for each
  chunk of replicates by a call to \code{ld.sim1}. Marker-trait  association is tested by 
  a one-way analysis of variance of trait values in terms of marker classes.
  The proportion of  results with P-value over a given threshold gives a stochastic estimate
  of the power calculated by \code{\link{luo.ld.power}}.
}
\value{
  If \code{data.only=FALSE}, an array with 1 row per simulation run, 
  and 4 columns with values for each run:
  \item{MS.beta}{Between marker classes mean square}
  \item{MS.within}{Within marker classes mean square}
  \item{F.value}{F value}
  \item{P.value}{P value}
  otherwise, if \code{data.only=TRUE}, an array nsim*n rows and  3 columns with 
   marker and trait values for simulated populations:
  \item{marker}{Marker genotype indicator with values \{1,2,3\}, corresponding to
    genotypes \{MM,Mm,mm\}}
  \item{y}{trait values}
  \item{replicate}{Replicate population indicator. Each blocks of rows with
  a given replicate number is a simulated population with the given parameters}
}
\references{Luo, Z.W. 1988 Detecting linkage disequilibrium between a
 polymorphic marker locus and a trait locus in natural populations.
 Heredity 80, 198--208.}
\author{Rod Ball \email{rod.ball@scionresearch.com} 
\url{www.scionresearch.com}}
%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{\code{\link{luo.ld.power}}}
\examples{
# Power from stochastic simulation for Luo's population 12
data(luo.ld.populations)
luo.pop12.sim <- ld.sim(nsim=3000,
                        n=luo.ld.populations[12,"n"],
                        p=luo.ld.populations[12,"p"],
                        q=luo.ld.populations[12,"q"],
                        D=luo.ld.populations[12,"D"],
                        h2=luo.ld.populations[12,"h2"],
                        phi=luo.ld.populations[12,"phi"],
                        Vp=100)
# power
table(luo.pop12.sim[,4] < 0.05)[2]/sum(table(luo.pop12.sim[,4] < 0.05))
# Cf power from deterministic calculation
luo.ld.power(n=luo.ld.populations[12,"n"],
             p=luo.ld.populations[12,"p"],
             q=luo.ld.populations[12,"q"],
             D=luo.ld.populations[12,"D"],
             h2=luo.ld.populations[12,"h2"],
             phi=luo.ld.populations[12,"phi"],
             Vp=100,
             alpha=0.05)
}
\keyword{design}% at least one, from doc/KEYWORDS

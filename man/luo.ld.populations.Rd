\name{luo.ld.populations}
\alias{luo.ld.populations}
\docType{data}
\title{Luo's Linkage disequilibrium example populations}
\description{
Matrix with rows containing parameters (population size, allele
frequencies, disequilibrium, dominance ratio) for example populations
with a bi-allelic marker and QTL in linkage disequilibrium, from Luo (1998).
}
\usage{data(luo.ld.populations)}
\format{
  The format is:
 num [1:12, 1:7] 1 2 3 4 5 6 7 8 9 10 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:7] "pop." "n" "p" "q" ...
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Luo, Z.W. 1988 Detecting linkage disequilibrium between a
 polymorphic marker locus and a trait locus in natural populations.
 Heredity 80, 198--208.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(luo.ld.populations)
luo.ld.populations
%## maybe str(luo.ld.populations) ; plot(luo.ld.populations) ...
}
\keyword{datasets}

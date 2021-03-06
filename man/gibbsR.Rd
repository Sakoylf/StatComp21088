% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp21088R.R
\name{gibbsR}
\alias{gibbsR}
\title{gibbs function}
\usage{
gibbsR(N, a, b, n)
}
\arguments{
\item{N}{sample size}

\item{a}{The first parameter of the distribution}

\item{b}{The second parameter of the distribution}

\item{n}{The third parameter of the distribution}
}
\value{
generated random numbers matrix \code{X}
}
\description{
a function of gibbs sample
}
\examples{
\dontrun{
C <- gibbsR(5000,1,1,25)
print(C)
}
}

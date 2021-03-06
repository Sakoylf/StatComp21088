% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp21088R.R
\name{SGD}
\alias{SGD}
\title{Stochastic Gradient Descent}
\usage{
SGD(input_data, real_result, alpha, theta)
}
\arguments{
\item{input_data}{Input_data  matrix after adding constant 1 column}

\item{real_result}{Real_result vector whose length is equal to the column number of data.}

\item{alpha}{Learning rate}

\item{theta}{The initial parameters of linear regression}
}
\value{
theta after iterations \code{theta}
}
\description{
Stochastic gradient descent.The idea is to use one sample per iteration to update the parameters.
}
\examples{
\dontrun{
x <- seq(0.1,10,0.002)
n <- length(x)
y <- 2*x+5+rnorm(n)
z <- as.matrix(data.frame(rep(1,n),x))
theta <- SGD(z, y,0.002,c(1,1))
print(theta)
}
}

\name{StatComp21088-package}
\alias{StatComp21088-package}
\alias{StatComp21088}
\docType{package}
\title{
  Achieve the Three kinds of gradient descent algorithm
}
\description{
  This package is developed to developed to achieve three kinds of gradient descent algorithm.Include Batch Gradient Descent,Stochastic Gradient Descent,Mini-batch gradient descent algorithm.And it will record the loss function value of the model after every iteration.
}
\details{
  Use train functions to train the model and get the result and the recorded cost.
}
\author{
Haifeng Chen

Maintainer: Haifeng Chen <cfhai@email.ustc.edu.cn>
}
\keyword{ package }
\examples{
  \dontrun{
    x <- seq(0.1,10,0.01)
    n <- length(x)
    z <- rnorm(n)
    y <- 2*x+5+z
    theta <- train_MBGD(x,y,500,200,0.01,c(1,1))$theta
    print(theta)
    cost <- train_MBGD(x,y,500,200,0.01,c(1,1))$cost
    print(cost)
  }
}


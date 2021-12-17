\name{StatComp21088-package}
\alias{StatComp21088-package}
\alias{StatComp21088}
\docType{package}
\title{
  Achieve the Mini-batch Gradient Descent algorithm
}
\description{
  This package is developed to add the first column which is full of 1 in input_data matrix,train the model one time in Mini-batch gradient descent algorithm and record the loss function value of the model after every iteration.
}
\details{
  Use train_MBGD function to train the model and get the result and the recorded cost.
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


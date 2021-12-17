#' @title gibbs function
#' @name gibbsR
#' @description a function of gibbs sample
#' @importFrom stats rbinom rbeta
#' @param N sample size
#' @param a The first parameter of the distribution
#' @param b The second parameter of the distribution
#' @param n The third parameter of the distribution
#' @return generated random numbers matrix \code{X}
#' @examples
#' \dontrun{
#' C <- gibbsR(5000,1,1,25)
#' print(C)
#' }
#' @export

gibbsR <- function(N, a, b, n){
  X <- matrix(0, N, 2)
  X[1,] <- c(0,0.5)
  for(i in 2:N){
    X2 <-  X[i-1, 2]
    X[i,1] <- rbinom(1,n,X2)
    X1 <- X[i,1]
    X[i,2] <- rbeta(1,X1+a,n-X1+b)
  }

  return(X)
}

#' @title Shuffle sequence of data
#' @name Shuffle_sequence
#' @description Before using the algorithm, the data needs to be randomly shuffled to eliminate the impact of the data order.
#' @param data Data matrix after adding constant 1 column.
#' @return New sequence  \code{random_sequence}
#' @examples
#' \dontrun{
#' A = matrix(
#' c(1,1,1,1,1,2,3,4),4,2)
#' x <- Shuffle_sequence(A)
#' print(x)
#' }
#' @export

Shuffle_sequence <- function(data){
  n <- nrow(data)
  random_sequence <- sample(1:n,n)
  return(random_sequence)
}


#' @title Mini-Batch Gradient Descent
#' @name MBGD
#' @description It is a compromise between batch gradient descent and stochastic gradient descent.The idea is to use batch_size samples per iteration to update the parameters.
#' @param input_data Input_data  matrix after adding constant 1 column
#' @param real_result Real_result vector whose length is equal to the column number of data.
#' @param batch_size Batch_size parameter constant
#' @param alpha Learning rate
#' @param theta The initial parameters of linear regression
#' @return theta after iterations \code{theta}
#' @examples
#' \dontrun{
#' x <- seq(0.1,10,0.002)
#' n <- length(x)
#' y <- 2*x+5+rnorm(n)
#' z <- as.matrix(data.frame(rep(1,n),x))
#' theta <- MBGD(z, y, 100,0.002,c(1,1))
#' print(theta)
#' }
#' @export

MBGD <- function(input_data,real_result,batch_size,alpha,theta){
  shuffle_sequence <- Shuffle_sequence(input_data)
  x <- input_data[shuffle_sequence,]
  y <- real_result[shuffle_sequence]
  n <- nrow(x)
  for(start in seq(1,n,batch_size)){

    end <- min(c(start+batch_size-1,n))
    Mini_train_data <- x[start:end,]
    Mini_train_result <- y[start:end]

    n0 <- nrow(Mini_train_data)
    p0 <- ncol(Mini_train_data)
    np <- rep(0,each = n0*p0)
    gradient_increasment <- matrix(np,nrow = n0, ncol = p0)

    for (i in 1:n0){
    g <- as.numeric(Mini_train_data[i,]*as.vector(Mini_train_result[i] - Mini_train_data[i,]%*%theta))

    gradient_increasment[i,] <- g
    }

    avg_g <- colMeans(gradient_increasment)
    theta <- theta + alpha * avg_g
  }
  return(theta)
}


#' @title Cost
#' @name Cost
#' @description loss function
#' @param input_data Input_data  matrix after adding constant 1 column
#' @param real_result Real_result vector whose length is equal to the column number of data.
#' @param theta The parameters of linear regression
#' @return cost \code{cost}
#' @examples
#' \dontrun{
#' x <- seq(0.1,10,0.002)
#' n <- length(x)
#' y <- 2*x+5+rnorm(n)
#' z <- as.matrix(data.frame(rep(1,n),x))
#' cost <- Cost(z, y, c(1,1))
#' print(cost)
#' }
#' @export

Cost <- function(input_data,real_result,theta){
  predict <- input_data%*%theta
  cost <- predict - real_result
  cost <- mean(cost^2)
  return(cost)
}


#' @title Using Mini-Batch Gradient Descent to train the model.
#' @name train_MBGD
#' @description Use Mini-Batch Gradient Descent to train the model.
#' @param input_data Input_data  matrix after adding constant 1 column
#' @param real_result Real_result vector whose length is equal to the column number of data.
#' @param iter iterations
#' @param batch_size Batch_size parameter constant
#' @param alpha Learning rate
#' @param theta The initial parameters of linear regression
#' @return theta after iterations and cost after every iteration \code{(theta ,cost)}
#' @examples
#' \dontrun{
#' x <- seq(0.1,10,0.01)
#' n <- length(x)
#' z <- rnorm(n)
#' y <- 2*x+5+z
#' theta <- train_MBGD(x,y,500,200,0.01,c(1,1))$theta
#' print(theta)
#' cost <- train_MBGD(x,y,500,200,0.01,c(1,1))$cost
#' print(cost)
#' }
#' @export

train_MBGD <- function(input_data,real_result,iter,batch_size,alpha,theta){
  input_data <- as.matrix(data.frame(rep(1,length(real_result)),input_data))
  cost <- numeric(nrow(input_data))
  for(i in 1:iter){
    theta <- MBGD(input_data,real_result,batch_size,alpha,theta)
    cost[i] <- Cost(input_data,real_result,theta)
  }
  return(list(theta = theta ,cost = cost))
}

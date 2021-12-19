## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#library(StatComp21088)
Shuffle_sequence <- function(data){
  n <- nrow(data)
  random_sequence <- sample(1:n,n)
  return(random_sequence)
}

n <- 10
A <- matrix(rep(1,n*n),n,n)
shuffle_sequence <- Shuffle_sequence(A)
print(shuffle_sequence)

## -----------------------------------------------------------------------------
Addo <- function(data){
    n <- nrow(as.matrix(data))
    c <- rep(1,n)
    new_data <- as.matrix(data.frame(c,data))
  
    return(new_data)
}

A = matrix(
c(1,1,1,1,1,2,3,4),4,2)
x <- Addo(A)
print(x)


## -----------------------------------------------------------------------------
BGD <- function(input_data,real_result,alpha,theta){
  n0 <- nrow(input_data)
  p0 <- ncol(input_data)
  gradient_increasment <- matrix(NA,nrow = n0, ncol = p0)
  for (i in 1:n0){
      g <- as.numeric(input_data[i,]*as.vector(real_result[i] - input_data[i,]%*%theta))
      gradient_increasment[i,] <- g
    }
    avg_g <- colMeans(gradient_increasment)
    theta <- theta + alpha * avg_g
  return(theta)
}

x <- seq(0.1,10,0.002)
n <- length(x)
y <- 2*x+5+rnorm(n)
z <- as.matrix(data.frame(rep(1,n),x))
theta <- BGD(z, y,0.002,c(1,1))
print(theta)

## -----------------------------------------------------------------------------
SGD <- function(input_data,real_result,alpha,theta){
  
  shuffle_sequence <- Shuffle_sequence(input_data)
  x <- input_data[shuffle_sequence,]
  y <- real_result[shuffle_sequence]
  n0 <- nrow(x)
  
  for (i in 1:n0){
    g <- as.numeric(x[i,]*as.vector(y[i] - x[i,]%*%theta))
    theta <- theta + alpha * g
  }
  return(theta)
}

x <- seq(0.1,10,0.002)
n <- length(x)
y <- 2*x+5+rnorm(n)
z <- as.matrix(data.frame(rep(1,n),x))
theta <- SGD(z, y,0.002,c(1,1))
print(theta)

## -----------------------------------------------------------------------------
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
    gradient_increasment <- matrix(NA,nrow = n0, ncol = p0)

    for (i in 1:n0){
    g <- as.numeric(Mini_train_data[i,]*as.vector(Mini_train_result[i] - Mini_train_data[i,]%*%theta))

    gradient_increasment[i,] <- g
    }

    avg_g <- colMeans(gradient_increasment)
    theta <- theta + alpha * avg_g
  }
  return(theta)
}

x <- seq(0.1,10,0.002)
n <- length(x)
y <- 2*x+5+rnorm(n)
z <- as.matrix(data.frame(rep(1,n),x))
theta <- MBGD(z, y, 100,0.002,c(1,1))
print(theta)


## -----------------------------------------------------------------------------

Cost <- function(input_data,real_result,theta){
  predict <- input_data%*%theta
  cost <- predict - real_result
  cost <- mean(cost^2)
  return(cost)
}


x <- seq(0.1,10,0.002)
n <- length(x)
y <- 2*x+5+rnorm(n)
z <- as.matrix(data.frame(rep(1,n),x))
cost <- Cost(z, y, c(1,1))
print(cost)

## -----------------------------------------------------------------------------
train_BGD <- function(input_data,real_result,iter,alpha,theta){
  input_data <- Addo(input_data)
  cost <- numeric(iter)
  for(i in 1:iter){
    theta <- BGD(input_data,real_result,alpha,theta)
    cost[i] <- Cost(input_data,real_result,theta)
  }
  return(list(theta = theta ,cost = cost))
}

x <- seq(0.1,10,0.01)
n <- length(x)
z <- rnorm(n)
y <- 2*x+5+z
theta <- train_BGD(x,y,500,0.02,c(1,1))
print(theta)


## -----------------------------------------------------------------------------
train_SGD <- function(input_data,real_result,iter,alpha,theta){
  input_data <- Addo(input_data)
  cost <- numeric(iter)
  for(i in 1:iter){
    theta <- SGD(input_data,real_result,alpha,theta)
    cost[i] <- Cost(input_data,real_result,theta)
  }
  return(list(theta = theta ,cost = cost))
}

x <- seq(0.1,10,0.01)
n <- length(x)
z <- rnorm(n)
y <- 2*x+5+z
theta <- train_SGD(x,y,200,0.02,c(1,1))
print(theta)


## -----------------------------------------------------------------------------
train_MBGD <- function(input_data,real_result,iter,batch_size,alpha,theta){
  input_data <- Addo(input_data)
  cost <- numeric(iter)
  for(i in 1:iter){
    theta <- MBGD(input_data,real_result,batch_size,alpha,theta)
    cost[i] <- Cost(input_data,real_result,theta)
  }
  return(list(theta = theta ,cost = cost))
}

x <- seq(0.1,10,0.01)
n <- length(x)
z <- rnorm(n)
y <- 2*x+5+z
theta <- train_MBGD(x,y,500,200,0.01,c(1,1))
print(theta)



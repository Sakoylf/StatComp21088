## -----------------------------------------------------------------------------
y <- runif(1000)
x <- -log(1-y) # 指数分布
hist(x,prob = TRUE)
z <- seq(0,8,.01)
lines(z,exp(-z))

## -----------------------------------------------------------------------------
sigma <- c(0.1,1,10)
for (i in 1:length(sigma)){
m <- 1000
u <- runif(m)
x <- sqrt(-2*sigma[i]^2*log(u))
{if(sigma[i]==0.1){
hist(x,freq=F,main=expression(sigma==0.1),col=rainbow(5))}
else if(sigma[i]==1)
{hist(x,freq=F,main=expression(sigma==1),col=rainbow(5))}
else
{hist(x,freq=F,main=expression(sigma==10),col=rainbow(5))}}
y <- seq(0,100,0.1)
lines(y,y/sigma[i]^2*exp(-y^2/(2*sigma[i]^2)))

}

## -----------------------------------------------------------------------------
m <- 1000
r1 <- rnorm(m)
r2 <- rnorm(m,mean=3)
p1 <- 0.75
p2 <- 1-p1
p0 <- sample(c(0,1),m,replace=TRUE,prob=c(p2,p1))
r <- p0*r1 + (1-p0)*r2
hist(r,freq=F,main=expression(p1==0.75),col=2)

p <- c(0.1,0.5,0.9)
for(i in 1:length(p)){
r01 <- rnorm(m)
r02 <- rnorm(m,3)
p00 <- sample(c(0,1),m,replace=TRUE,prob=c(1-p[i],p[i]))
r0 <- p00*r01 + (1-p00)*r02
if(p[i]==0.1){
  hist(r0,freq=F,main=expression(p1==0.1),col=i)
}else if(p[i]==0.5){
  hist(r0,freq=F,main=expression(p1==0.5),col=i)
}else{
    hist(r0,freq=F,main=expression(p1==0.9),col=i)
  }

}

## ----echo=FALSE---------------------------------------------------------------
lambda <- c(1,2,3,4,5)
t0 <- 10

r <-c(2,9,6,4,7)
beta <- c(0.5,3,2,0.8,7)
for (i in 1:length(r)){
  
  X <- numeric(10000)
for (j in 1:10000) {
N <- rpois(1,lambda[i]*t0)
Tn <- rgamma(1000,r[i],beta[i])
Xn <- sum(Tn[1:N])
X[j] <- Xn
}
thrye <- lambda[i]*t0*(r[i]/beta[i]) #由公式得出的x(t)的期望
thryv <- lambda[i]*t0*((r[i]/beta[i])^2+r[i]/beta[i]^2)# 由公式得到的x(t)的方差
if (i==1){
  print('理论期望    样本均值  理论方差   样本方差')
}
print(c(thrye,mean(X),thryv,var(X)))#输出的顺序从左到右为理论期望 样本均值 理论方差 样本方差
}

## ----echo=FALSE---------------------------------------------------------------
betacdf <-function(x,y,z){
  if (length(x)==1){
 if(0<x&x<1){
   u <- runif(10000)
theta <- (u^2)*(x^3)*(1-u*x)^2
beta0 <- mean(theta)/beta(y,z)
return(beta0)
 }
  else if(x<0){
   return(0)
 }else{
   return(1)
 }}
else{
  beta1 <- numeric(length(x))
  for (i in 1:length(x)){
    beta1[i]=betacdf(x[i],y,z)
  }
  return(beta1)
}
}
x <- seq(0.1,0.9,length = 9)
cdf <- betacdf(x,3,3)
pbeta <- pbeta(x,3,3)
print(round(rbind(x,cdf,pbeta),3))

## -----------------------------------------------------------------------------
sigma <- 0.1

m <- 100
u <- runif(m)
x <- sqrt(-2*sigma^2*log(u))
xp <- sqrt(-2*sigma^2*log(1-u))
print(round(rbind(x,xp),3))

MCa <- function(x, R = 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic){
  v <- runif(R/2) 
}else{
v <- 1 - u
}
u <- c(u, v)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
g <- u*x[i]*exp(-u*x[i]/(sigma^2))
cdf[i] <- mean(g) / (sigma^2)
}
cdf
}
n <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.95
for (i in 1:m) {
MC1[i] <- MCa(x, R = 1000, anti = FALSE)
MC2[i] <- MCa(x, R = 1000)
}
print((var(MC1) - var(MC2))/var(MC1))


## -----------------------------------------------------------------------------
x <- seq(1,20,by = 0.1)
g <- x^2/sqrt(2*pi)*exp((-x^2)/2)
f1 <- 1/sqrt(2*pi)*exp((-x^2)/2)
f2 <- x*exp((-x^2)/2)
plot(x,g)
lines(x,f1,col = 'red')#f1是红色
lines(x,f2,col = 'blue')#f2是蓝色


## -----------------------------------------------------------------------------
x <- abs(runif(10000))
y <- rexp(10000)
fg1 <- x^2
fg2 <- y/sqrt(2*pi)
theta1 <- mean(fg1)
theta2 <- mean(fg2)
se1 <- sd(fg1)
se2 <- sd(fg2)
theta <- c(theta1,theta2)
se <- c(se1,se2)
print(rbind(theta,se))

## -----------------------------------------------------------------------------
x <- abs(runif(10000))
fg1 <- x^2
theta1 <- mean(fg1)
se1 <- sd(fg1)
print(c(theta1,se1))

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05

set.seed(1)

# calculate the half length of CIs
CL <- replicate(1000, expr = {
  x <- rchisq(n, df = 2)
  mu.hat <- mean(x)
  sd.hat <- sd(x)
  sqrt(n) * abs(mu.hat - 2) / sd.hat
})

# calculate coverage proportion
mean(CL < qt(1 - alpha / 2, df = n - 1))

## -----------------------------------------------------------------------------
set.seed(1)
R <- 1000
n <- 20
alpha <- 0.05
cover <- numeric(R)
for (i in 1:R) {
  x <- rchisq(n, df = 2)
  UpperBound <- (n - 1) * var(x) / qchisq(alpha, df = n - 1)
  cover[i] <- (UpperBound > 4)
}
(ECP <- sum(cover) / R)

## ----echo=TRUE----------------------------------------------------------------
library(MASS)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)  #协方差矩阵
n <- c(10,20,30,50,100,200,500) #样本量
d <- 3 # 维度
cv <- qchisq(.95, d*(d+1)*(d+2)/6) #置信界

b1d <- function(x){
  k <- nrow(x)
  b1dm <- matrix(rep(0,k^2),nrow=k,ncol = k)
  x.lc <- matrix(rep(0,d*k),nrow=k,ncol = d)
  for (t in 1:d){
    x.lc[,t] <- x[,t]-mean(x[,t])
  } 
  sigma.hat <- t(x.lc)%*%x.lc/k
  
  b1dm <- (x.lc%*%solve(sigma.hat)%*%t(x.lc))^3
 
  return(sum(b1dm)/k^2)
}

p.reject <- numeric(length(n)) #储存结果
m <- 1000 #每个样本量重复m次实验
for (i in 1:length(n)) {
mst <- numeric(m) #储存每次实验的检验结果
for (j in 1:m) {
x <- mvrnorm(n[i],rep(0,3), sigma)
mst[j] <- as.integer(n[i]*b1d(x)/6 >= cv)#拒绝原假设为1，左侧拒绝域
}
p.reject[i] <- mean(mst)#拒绝原假设的频率
}
print(p.reject)

## ----echo=TRUE----------------------------------------------------------------


#sigma <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)  #协方差矩阵
n <- 30 #样本量
#d <- 3 # 维度
cv <- qchisq(.9, d*(d+1)*(d+2)/6) #置信界



p <- seq(0,1,0.02)
pow <- numeric(length(p)) #储存结果
m <- 100 #每个p重复m次实验
for (i in 1:length(p)) {
mst <- numeric(m) #储存每次实验的检验结果
  for (j in 1:m) {
    x <- matrix(rep(0,3*n),n,3)
    for (t in 1:n){
    sigma1 <- sample(c(1,10),size = 1,replace = TRUE,prob = c(p[i], 1-p[i]))
    sigma2 <- diag(sigma1,3,3)
    x[t,] <- mvrnorm(1,rep(0,3), sigma2)
    }
  mst[j] <- as.integer(n*b1d(x)/6 >= cv)#拒绝原假设为1，左侧拒绝域
  }
pow[i] <- mean(mst)#拒绝原假设的频率
}
plot(p, pow, type = "b",xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pow * (1-pow) / m) 
lines(p, pow+se, lty = 3)
lines(p, pow-se, lty = 3)

## ----echo=TRUE----------------------------------------------------------------
library(bootstrap)
set.seed(12121)
n <- nrow(scor)
m <- 1000
lambda.hat <- eigen(cov(scor))$values
theta.hat <- lambda.hat[1]/sum(lambda.hat)


theta.star <- numeric(m)
for (i in 1:m){
  index <- sample(1:n,n,replace = TRUE) 
  lambda.star <- eigen(cov(scor[index,]))$values
  theta.star[i] <- lambda.star[1]/sum(lambda.star)
}

c(theta.hat= theta.hat,boot.bias_theta=mean(theta.star) - theta.hat,boot.se_theta=sd(theta.star))

## ----echo=TRUE----------------------------------------------------------------
theta.j <- numeric(n)
for (j in 1:n){
  lambda <- eigen(cov(scor[-j,]))$values
  theta.j[j] <- lambda[1]/sum(lambda)
}
c(jack.bias_theta=(n-1)*(mean(theta.j) - theta.hat),jack.se_theta=sqrt((n-1)*mean((mean(theta.j)-theta.j)^2)))

## ----echo=TRUE----------------------------------------------------------------
library(boot)
set.seed(12121)
theta.boot <- function(data,i){
  lamdb <- eigen(cov(data[i,]))$values
  lamdb[1] / sum(lamdb)
}
boot.obj <- boot(scor,statistic = theta.boot,R = 2000)
print(boot.ci(boot.obj,type = c("perc","bca")))

## ----echo=TRUE----------------------------------------------------------------
set.seed(12345)
nor <- 999
x <- rchisq(10,10)
y <- rexp(10,0.1)
z <- c(x, y) 
N <- 1:20
cor <- numeric(nor)

cor0 <- cor(x, y, method = "spearman")
p0 <- cor.test(x, y,method = "spearman")$p.value
for (i in 1:nor) {
  n <- sample(N, size = 10, replace = FALSE)
  xs <- z[n]
  ys <- z[-n] #complement of x1
  cor[i] <- cor(xs, ys, method = "spearman")
}
p <- mean(abs(c(cor0, cor)) > abs(cor0))
c("cor.test" = p0, "permutation" = p)

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
NNF <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 <= n1); i2 <- sum(block2 > n1)
(i1 + i2) / (k * n)
}


p <- function(x,y){
  z <- rbind(x,y)
  N <- c(nrow(x),nrow(y))

  energyp.value <- eqdist.etest(z, sizes=N, R=999)$p.value
  
  ballp.value <- bd.test(x = x, y = y, num.permutations=999)$p.value
  
  boot.obj <- boot(z,statistic = NNF, R = 999, sim =   "permutation", sizes = N, k = 3)
  ann <- c(boot.obj$t0,boot.obj$t)
  NNp.value <- mean(ann>=ann[1])
  
  c(NNp.value, energyp.value, ballp.value)
}



## ----echo=TRUE----------------------------------------------------------------
library(MASS)
set.seed(12345)
e1 = c(0,0)
cov1 = matrix(c(1,0,0,1),nrow=2,ncol=2)
e2 = c(0,0)
cov2 = matrix(c(3,0,0,3),nrow=2,ncol=2)

p.values <- matrix(NA,100,3)
  for (i in 1:100){
    x <- mvrnorm(20,e1,cov1)
    y <- mvrnorm(20,e2,cov2)
    p.values[i,] <- p(x,y)
  }
  pwr1 <- colMeans(p.values<0.05)#alpha取0.05
  c("NN.power" = pwr1[1], "energy.power" = pwr1[2],"ball.power" = pwr1[3])


## ----echo=TRUE----------------------------------------------------------------
e1 = c(0,0)
cov1 = matrix(c(1,0,0,1),nrow=2,ncol=2)
e2 = c(1,1)
cov2 = matrix(c(2,0,0,2),nrow=2,ncol=2)

p.values <- matrix(NA,100,3)
  for (i in 1:100){
    x <- mvrnorm(20,e1,cov1)
    y <- mvrnorm(20,e2,cov2)
    p.values[i,] <- p(x,y)
  }
  pwr1 <- colMeans(p.values<0.05)#alpha取0.05
  c("NN.power" = pwr1[1], "energy.power" = pwr1[2],"ball.power" = pwr1[3])

## ----echo=TRUE----------------------------------------------------------------

p.values <- matrix(NA,100,3)
  for (i in 1:100){
    x <- as.matrix(rt(20,1,2),ncol=1)
    y <- as.matrix(rt(20,2,5),ncol=1)
    p.values[i,] <- p(x,y)
  }
  pwr1 <- colMeans(p.values<0.05)#alpha取0.05
  c("NN.power" = pwr1[1], "energy.power" = pwr1[2],"ball.power" = pwr1[3])

## ----echo=TRUE----------------------------------------------------------------

rbimodel <- function(n,e1,e2,sd1,sd2){
  index <- sample(1:2,n,replace=TRUE)
  x <- numeric(n)
  index1<- which(index==1)
  x[index1] <- rnorm(length(index1), e1, sd1)
  index2 <- which(index==2)
  x[index2] <- rnorm(length(index2), e2, sd2)
  return(x)
}

p.values <- matrix(NA,100,3)
  for (i in 1:100){
    x <- as.matrix(rbimodel(20,0,0,1,2),ncol=1)
    y <- as.matrix(rbimodel(20,1,1,1,2),ncol=1)
    p.values[i,] <- p(x,y)
  }
  pwr1 <- colMeans(p.values<0.05)#alpha取0.05
  c("NN.power" = pwr1[1], "energy.power" = pwr1[2],"ball.power" = pwr1[3])

## ----echo=TRUE----------------------------------------------------------------
e1 = c(0,0)
cov1 = matrix(c(1,0,0,1),nrow=2,ncol=2)
e2 = c(1,1)
cov2 = matrix(c(2,0,0,2),nrow=2,ncol=2)

p.values <- matrix(NA,100,3)
  for (i in 1:100){
    x <- mvrnorm(10,e1,cov1)
    y <- mvrnorm(100,e2,cov2)
    p.values[i,] <- p(x,y)
  }
  pwr1 <- colMeans(p.values<0.05)#alpha取0.05
  c("NN.power" = pwr1[1], "energy.power" = pwr1[2],"ball.power" = pwr1[3])

## -----------------------------------------------------------------------------
N <- 10000
X <- numeric(N)
b <- 1001      #discard the burn-in sample
X[1] <- rnorm(1,0,1)
for(i in 2:N){
  Xt <- X[i-1]
  Y <- rnorm(1,0,abs(Xt))
  r <- dt(Y,1)*dnorm(Xt,0,abs(Y))/dt(Xt,1)/dnorm(Y,0,abs(Xt))
  U <- runif(1)
  if(r > 1) r <- 1
  if(U <= r) X[i] <- Y
  else X[i] <- Xt
}
Y <- X[b:N]
deciles_sample <- quantile(Y, c(1:9)/10)   # 样本的十分位数
deciles_true <- qcauchy(c(1:9)/10,0,1)   # 真实分布的十分位数
deciles <- rbind(deciles_sample,deciles_true)
rownames(deciles) <- c("deciles_sample","deciles_true")
colnames(deciles) <- c(c(1:9)/10)
library(knitr)
kable(deciles)

## -----------------------------------------------------------------------------
a <- ppoints(100)
QR <- qcauchy(a,0,1)
Q <- quantile(Y, a)
qqplot(QR, Q,  main="",
xlab="Standard Cauchy Quantiles", ylab="Sample Quantiles")
abline(0,1,col = "red")

## -----------------------------------------------------------------------------
a <- 1
b <- 1
N <- 10000          #样本量
X <- matrix(0, N, 2)  #样本阵
X[1,] <- c(0,0.5)
for(i in 2:N){
  X2 <-  X[i-1, 2]
  X[i,1] <- rbinom(1,25,X2)
  X1 <- X[i,1]
  X[i,2] <- rbeta(1,X1+a,25-X1+b)
}
plot(X[,1],X[,2],xlab = "x",ylab = "y")

## -----------------------------------------------------------------------------
a <- 1
b <- 10
X <- matrix(0, N, 2)  #样本
X[1,] <- c(0,0.5)
for(i in 2:N){
  X2 <-  X[i-1, 2]
  X[i,1] <- rbinom(1,25,X2)
  X1 <- X[i,1]
  X[i,2] <- rbeta(1,X1+a,25-X1+b)
}
plot(X[,1],X[,2],xlab = "x",ylab = "y")

## -----------------------------------------------------------------------------
# 计算Gelman-Rubin statistic的函数
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)

  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
        }

## -----------------------------------------------------------------------------
# 生成标准柯西分布的Metropolis chain
# 提议函数仍取9.3中使用的对称正态分布 N(0,X[t]^2)
# X1为初始值
Standard_Cauchy_Chain <- function(N, X1){
  X <- numeric(N)
  X[1] <- X1    #初始值
  for(i in 2:N){
    Xt <- X[i-1]
    Y <- rnorm(1,0,abs(Xt))
    r <- dt(Y,1)*dnorm(Xt,0,abs(Y))/dt(Xt,1)/dnorm(Y,0,abs(Xt))
    U <- runif(1)
    if(r > 1) r <- 1
    if(U <= r) X[i] <- Y
    else X[i] <- Xt
  }
  return(X)
}

## -----------------------------------------------------------------------------
k <- 4      
N <- 8000
b <- 1000     #burn-in length
X1 <- c(0.1,0.2,0.1,0.2)    #初始值

# 生成4条样本
set.seed(12345)
X <- matrix(0, nrow = k, ncol = N)
for(i in 1:k){
  X[i,] <- Standard_Cauchy_Chain(N, X1[i])
}

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
# 四条样本的psi
for (i in 1:k)
  if(i==1){
    plot((b+1):N,psi[i, (b+1):N],ylim=c(-1,1), type="l",
         xlab='Index', ylab=bquote(phi))
  }else{
      lines(psi[i, (b+1):N], col=i)
  }
par(mfrow=c(1,1)) 

## -----------------------------------------------------------------------------
par(mfrow=c(1,1)) 
#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
# 生成二元随机变量的Gibbs sampler
# X1为初始值
Bivariate.Gibbs <- function(N, X1){
  a <- b <- 1
  X <- matrix(0, N, 2)
  X[1,] <- X1    #初始值
  for(i in 2:N){
    X2 <-  X[i-1, 2]
    X[i,1] <- rbinom(1,25,X2)
    X1 <- X[i,1]
    X[i,2] <- rbeta(1,X1+a,25-X1+b)
  }
  return(X)
}

## -----------------------------------------------------------------------------
k <- 4          
N <- 8000 
b <- 1000    #burn-in length
X1 <- cbind(c(2,7,10,15),runif(4)) #初始值

#生成4条样本，每个第一维的放在X中，第二维的放在Y中
set.seed(12345)
X <- matrix(0, nrow=k, ncol=N)
Y <- matrix(0, nrow=k, ncol=N)
for (i in 1:k){
  BG <- Bivariate.Gibbs(N, X1[i,])
  X[i, ] <- BG[,1]
  Y[i, ] <- BG[,2]
}

## -----------------------------------------------------------------------------
# 先考虑第一维样本X

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
# 再考虑第二维样本Y

#compute diagnostic statistics
psi <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
cal1 <- function(a,d,k,n){
  jc <- 1
  asquare <- sum(a^2)
    i <- 0:n
    j <- rep(c(1,-1), length=n+1)
    A <- ((j * asquare^(i+1)) / ((2*i+1)*(2*i+2)*factorial(i)*2^i)) * exp(lgamma((d+1)/2)+lgamma(i+1.5)-lgamma(i+d/2+1))
  return(A[k+1])
}

## -----------------------------------------------------------------------------
cal2 <- function(a,d,k,n){
  jc <- 1
  asquare <- sum(a^2)
    i <- 0:n
    j <- rep(c(1,-1), length=n+1)
    A <- ((j * asquare^(i+1)) / ((2*i+1)*(2*i+2)*factorial(i)*2^i)) * exp(lgamma((d+1)/2)+lgamma(i+1.5)-lgamma(i+d/2+1))
  return(c(A[k+1],sum(A)))
}


## -----------------------------------------------------------------------------
a <- c(1,2)
d <- 1
res1 <- cal2(a,d,10,10)
res1
res2 <- cal2(a,d,100,100)
res2
res3 <- cal2(a,d,200,200)
res3

## -----------------------------------------------------------------------------
integrand <- function(u,k0){
  (1+(u^2)/(k0-1))^(-k0/2)
}
g <- function(a,k){
  up <- sqrt(a^2 * (k-1)/(k-a^2))
  inte <- integrate(integrand,lower=0,upper=up,k0 = k)$value
  (2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2)))*inte
}
f <- function(a,k){
 g(a,k+1)-g(a,k)
}
f1 <- function(a,k){
  C <- numeric(length(a))
  for(i in 1:length(a)){
  C[i] <- f(a[i],k)
  }
  return(C)
}
k <- c(16:25,50,100,300)
a0 <- seq(0, 4, by=0.01)
plot(a0, f1(a0, k[1]), lty=1, col=1, type="l", xlim=c(0, 4),xlab="a", ylab="f(a|k)", main="f(a) with different k")
lines(a0, f1(a0, k[7]), xlim = c(0, 4), lty=2, col=2)
legend("topright", legend=c("k=16", "k=25"), col=1:2,lty=1:2)

## ----echo=TRUE----------------------------------------------------------------
sol <- function(k1){
 m <-uniroot(f,k=k1,lower = 1,upper = 2)$root
 m
}
B <- numeric(length(k))
for (i in 1:length(k)){
B[i] <- sol(k[i])
}
# 11.4
S = function(a,k){
 ck = sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

root = matrix(0,3,length(k))

for (i in 1:length(k)){
  root[2,i]=round(solve(k[i]),4)
}

root[1,] = k
root[3,] = B
rownames(root) = c('k','A(k)','B(k)')
root

## ----echo=TRUE----------------------------------------------------------------
LL <- function(lambda, y) {
f <- numeric(length(y))
for (i in 1:length(y)){
  if(y[i]==1){
    f[i] <- 1-pexp(y[i],rate=lambda)
  }else{
    f[i] <- dexp(y[i],rate=lambda)
  }
}
return( -sum(log(f)))
}
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
m <- 2000
lambda <- 1/mean(y)
opt <- optim(lambda, LL, y=y)
MLE <- opt$par
MLE


## ----echo=TRUE----------------------------------------------------------------
N <- 10000 #max. number of iterations
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
L <- 1/mean(y) #initial est. for lambdas
tol <- .Machine$double.eps^0.5
L.old <-L+1
for(i in 1:length(y)){
  if(y[i]==1){
    y[i] <- 0
  }
}
beta0 <- sum(y)
l1 <- length(y)
for (j in 1:N) {
L <-  l1/ (beta0+3*(L+1)/L)
if (sum(abs(L - L.old)/L.old) < tol) break
L.old <- L
}
print(list(lambda = L, iter = j, tol = tol))
EM <- L
c(MLE = MLE , EM = EM)

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
ex3la <- lapply(formulas, lm, data = mtcars)
ex3lf <- vector("list", length(formulas))
for (i in seq_along(formulas)){
  ex3lf[[i]] <- lm(formulas[[i]], data = mtcars)
}
rsq <- function(mod) summary(mod)$r.squared
sapply(ex3la,rsq)
sapply(ex3lf,rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
ex4la <- lapply(bootstraps, lm, formula = mpg ~ disp)
ex4lf <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  ex4lf[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
rsq <- function(mod) summary(mod)$r.squared
sapply(ex4la, rsq)
sapply(ex4lf, rsq)


## -----------------------------------------------------------------------------
library(microbenchmark)
library(Rcpp)
set.seed(12121)
cppFunction(
'NumericMatrix gibbsC(int N, int a, int b, int n) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  mat(0,0) = 0;
  mat(0,1) = 0.5;
  for(int i = 1; i < N; i++) {

    x = rbinom(1, n, mat(i-1,1))[0];
    mat(i, 0) = x;
    y = rbeta(1, x+a, n-x+b)[0];
    mat(i, 1) = y;
  }
  return(mat);
}
')

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

N <- 5000
C <- gibbsC(N,1,1,25)
R <- gibbsR(N,1,1,25)
Cg <- C[1001:N,]
Rg <- R[1001:N,]
a <- ppoints(100)
QC1 <- quantile(Cg[,1], a)
QR1 <- quantile(Rg[,1], a)

QC2 <- quantile(Cg[,2], a)
QR2 <- quantile(Rg[,2], a)

qqplot(QC1, QR1,  main="",
xlab="C.gibbs1 Quantiles", ylab="R.gibbs1 Quantiles")
abline(0,1,col = "red")

qqplot(QC2, QR2,  main="",
xlab="C.gibbs2 Quantiles", ylab="R.gibbs2 Quantiles")
abline(0,1,col = "red")

time <- microbenchmark(gibbsC(N,1,1,25),gibbsR(N,1,1,25))
summary(time)[,c(1,3,5,6)]


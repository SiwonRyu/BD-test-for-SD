#install.packages("foreach")
#install.packages("doParallel")
library("foreach")
library("doParallel")

rm(list=ls(all=TRUE)) 

# Set parameters
alpha <- 0.05
ngrid <- 100
B <- 100
R <- 1
n1 <- 100
n2 <- 100

# Design
design <- c(0,1,0.3,1)
sample1 <- matrix(rnorm(n1,design[1],design[2]), ncol=1)
sample2 <- matrix(rnorm(n2,design[3],design[4]), ncol=1)
grid <- matrix(seq(min(rbind(sample1,sample2)),max(rbind(sample1,sample2)),length=ngrid))
s <- 1

# BD-test function
# H0 : X1 s-th order SD F2
# Input parameter : sample1, sample2, grid, s(SD order)
# Output : 1x5 vector (BD-teststatistic, pvalues of m1,m2,b1,b2,b3)
bdtest <-function(sample1,sample2,grid,s){
  operator <- function(x,grid,s){
    return(t(apply(x,1,'<',grid)))
  }
  ecdf <- function(x,grid,s){
    return(matrix(colMeans(operator(x,grid,s)),nrow=1))
  }
  stat <-function(sample1,sample2,grid,s){
    n1 = dim(sample1)[1]
    n2 = dim(sample2)[1]
    stat = sqrt(n1*n2/(n1+n2)) * max( ecdf(sample1,grid,s) - ecdf(sample2,grid,s) )
    return(stat)
  }
  op_J <- function(x,grid,s){
    n = dim(x)[1]
    #u <- array(runif(n*B),dim = c(n,1,B))
    u <- array(rnorm(n*B,0,1),dim=c(n,1,B))
    y <- array(apply(operator(x,grid,s),1,'-',ecdf(x,grid,s)),dim=c(1,ngrid,B))
    op_J <- array(0, dim=c(n,ngrid,B))
    for(i in 1:B){
      for(j in 1:ngrid){
        op_J[,j,i]<- y[1,j,i]*u[,1,i]
      }
    }
    opJ <- array(sqrt(n)*colMeans(op_J),dim=c(1,ngrid,B))
    return(opJ)
  }
  multi1 <- function(sample1,sample2,grid,s){
    n1 = dim(sample1)[1]
    n2 = dim(sample2)[1]
    lbd = n2/(n1+n2)
    mstat <- sqrt(lbd)*op_J(sample1,grid,s)-sqrt(1-lbd)*op_J(sample2,grid,s)
    #return(array(apply(m,3,max),dim=c(1,1,B)))
    return(array(apply(mstat,3,max),dim=c(1,1,B)))
  }
  multi2 <- function(sample1,sample2,grid,s){
    mstat <- op_J(sample2,grid,s)
    return(array(apply(mstat,3,max),dim=c(1,1,B)))
  }
  boot1 <- function(sample1,sample2,grid,s){
    n2 = dim(sample2)[1]
    bstat <- array(0, dim=c(1,ngrid,B))
    for(i in 1:B){
      bsample2 <- array(sample(sample2, replace=T),dim=c(n2,1))
      bstat[,,i] <- sqrt(n2)*(ecdf(bsample2,grid,s)-ecdf(sample2,grid,s))
    }
    return(array(apply(bstat,3,max),dim=c(1,1,B)))
  }
  boot2 <- function(sample1,sample2,grid,s){
    n1 = dim(sample1)[1]
    n2 = dim(sample2)[1]
    bstat <- array(0, dim=c(1,ngrid,B))
    for(i in 1:B){
      bsample1 <- array(sample(rbind(sample1,sample2),n1, replace=T),dim=c(n1,1))
      bsample2 <- array(sample(rbind(sample1,sample2),n1, replace=T),dim=c(n2,1))
      bstat[,,i] <- sqrt(n1*n2/(n1+n2))*(ecdf(bsample1,grid,s)-ecdf(bsample2,grid,s))
    }
    return(array(apply(bstat,3,max),dim=c(1,1,B)))
  }
  boot3 <- function(sample1,sample2,grid,s){
    n1 = dim(sample1)[1]
    n2 = dim(sample2)[1]
    bstat <- array(0, dim=c(1,ngrid,B))
    for(i in 1:B){
      bsample1 <- array(sample(sample1, replace=T),dim=c(n1,1))
      bsample2 <- array(sample(sample2, replace=T),dim=c(n2,1))
      
      bstat[,,i] <- 
        sqrt(n1*n2/(n1+n2))*(
          (ecdf(bsample1,grid,s)-ecdf(sample1,grid,s))
          -(ecdf(bsample2,grid,s)-ecdf(sample2,grid,s))
        )
    }
    return(array(apply(bstat,3,max),dim=c(1,1,B)))
  }
  
  bd <- stat(sample1,sample2,grid,1)
  
  pval_m1 <- mean(multi1(sample1,sample2,grid,1)>bd)
  pval_m2 <- mean(multi2(sample1,sample2,grid,1)>bd)
  pval_b1 <- mean(boot1(sample1,sample2,grid,1)>bd)
  pval_b2 <- mean(boot2(sample1,sample2,grid,1)>bd)
  pval_b3 <- mean(boot3(sample1,sample2,grid,1)>bd)
  
  return(c(bd,pval_m1,pval_m2,pval_b1,pval_b2,pval_b3))
  
}
result = bdtest(sample1,sample2,grid,s)
result



# Simulation part
# init <- Sys.time()
# R = 10
# Cl <- makeCluster(10,type = "PSOCK")
# registerDoParallel(Cl)
# rej_pre <- foreach(i = 1:R) %dopar%{
#   sample1 <- matrix(rnorm(n1,design[1],design[2]), ncol=1)
#   sample2 <- matrix(rnorm(n2,design[3],design[4]), ncol=1)
#   grid <- matrix(seq(min(rbind(sample1,sample2)),max(rbind(sample1,sample2)),length=ngrid))
#   result <- bdtest(sample1,sample2,grid,s)
#   return(result[2:6])
#   #result[3]
# }
# stopCluster(Cl)
# rej = matrix(unlist(rej_pre), nrow = R, byrow = TRUE)
# sprintf("Rejection probability of m1 = %1.3f" , mean(rej[,1]))
# sprintf("Rejection probability of m2  = %1.3f", mean(rej[,2]))
# sprintf("Rejection probability of b1  = %1.3f", mean(rej[,3]))
# sprintf("Rejection probability of b2  = %1.3f", mean(rej[,4]))
# sprintf("Rejection probability of b3  = %1.3f", mean(rej[,5]))
# Sys.time() - init



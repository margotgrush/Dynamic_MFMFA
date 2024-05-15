################################################################################
#              Code to generate simulated synthetic data sets                  #
################################################################################
rm(list=ls())
# Simulaton study 1 
# Simulated data set with 3 clusters and 4 factors in each cluster
{
  library(MASS)
  library(mvtnorm)
  # simulated data set K = 3, H = 3
  # cluster weights
  w = c(1/3, 1/3, 1/3)
  #number of clusters, factors and variables
  G = 3; J=4; P=20
  # number of observations
  T=100
  # factors
  F <- t(mvrnorm(T, mu=rep(0,J), Sigma=diag(1,J)))
  # loading matrices
  Lambda=array(NA, dim=c(P,J,G))
  for (g in 1:G) {
    for (i in 1:P) { 
      Lambda[i,,g] <- mvrnorm(1, mu=rep(0,J), Sigma=diag(1,J))
    }
  }
  # generating idiosyncratic variances
  Sigma1 = array(NA, dim=c(P,P,G))
  for (g in 1:G) {Sigma1[,,g] <- diag(x=1/rgamma(P, 2, rate=1), P)}
  # generating cluster means
  mclg <- matrix(NA, nrow = G,ncol=P)
  for (g in 1:G) {mclg[g,] = mvrnorm(1, mu=rep(2*g-G-1,P), Sigma=diag(1,P))} # mixing the means so that clusters are not very distinguishable
  # filling in clusters with T observations
  Tg <- array(NA, dim=c(P,T,G))
  for ( g in 1:G) {
    for (t in 1:T) {
      Tg[,t,g] = mvrnorm(1, mu=mclg[g,]+Lambda[,,g]%*%F[,t], Sigma=Sigma1[,,g])
    }}
  # sample observation from different clusters to form the mixture according to weights
  {
    z <- sample(1:G, T, replace = TRUE, prob = w)
    Y = matrix(NA, nrow=P, ncol=T)
    for ( i in 1:T) {
      Y[,i] = Tg[,i,z[i]]
    }
  }
  #calculate initial cluster covariance matrices
  Omega0_orig <- array(NA, dim=c(P,P,G))
  for (g in 1:G) {
    Omega0_orig[,,g]  <- Lambda[,,g]%*%t(Lambda[,,g]) + Sigma1[,,g]
  }
  
  p=P
  Nz = tabulate(z,G)
  
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
  
  # calculate correlation matrix (adjust to rescaling)
  Omega0 <- array(NA, dim=c(P,P,G))
  for (g in 1:G) {
    Dm <- diag(Ys)
    Omega0[,,g]  <- solve(Dm)%*%Omega0_orig[,,g]%*%solve(Dm)
  }
  # rescale factor loading matrices
  Lambda0 <- array(NA, dim=c(P,J,G))
  for (g in 1:G) {
    Dm <- diag(Ys)
    Lambda0[,,g]  <- solve(Dm)%*%Lambda[,,g]
  }
  
  # rescale the matrix of uniquenesses
  Sigma_stand <- array(NA, dim=c(P,P,G))
  for (g in 1:G) {
    Dm <- diag(Ys)
    Sigma_stand[,,g]  <- solve(Dm)%*%Sigma1[,,g]%*%solve(Dm)
  }
}

# Simulaton study 2
# Simulated data set with unbalanced clusters - dense loadings
# 6 clusters of various sizes and with various number of factors
{
  library(MASS)
  library(mvtnorm)
  # simulated data set K = 6, H = drawn randomly from (0, 6)
  # cluster weights
  w = c(0.25, 0.25, 0.2, 0.15, 0.1, 0.05)
  #number of clusters, factors and variables
  G = 6; P=20
  # number of observations
  T=700
  # draw randomly J
  J = rep(0,G)
  J=floor(runif(6, min=1, max=7)) #max=(P-1)/2)
  print(c("number of factors",J))
  # factors
  F <- list()
  for (g in 1:G) {
    F[[g]] <- matrix(NA,nrow=J[g], ncol=T)
    F[[g]] <- t(mvrnorm(T, mu=rep(0,J[g]), Sigma=diag(1,J[g])))
  }
  # loading matrices
  Lambda=list()
  for (g in 1:G) {
    Lambda[[g]] <- matrix(NA,nrow=P, ncol=J[g])
    for (i in 1:P) { 
      Lambda[[g]][i,] <- mvrnorm(1, mu=rep(0,J[g]), Sigma=diag(2,J[g]))
    }
  }
  
  # generating idiosyncratic variances
  Sigma1 = array(NA, dim=c(P,P,G))
  for (g in 1:G) {Sigma1[,,g] <- diag(x=1/rgamma(P, 2, rate=1), P)}
  # generating cluster means
  mclg <- matrix(NA, nrow = G,ncol=P)
  for (g in 1:G) {mclg[g,] = mvrnorm(1, mu=rep(2*g-G-1,P), Sigma=diag(1,P))} # mixing the means so that clusters are not very distinguishable
  # filling in clusters with T observations
  Tg <- array(NA, dim=c(P,T,G))
  for ( g in 1:G) {
    for (t in 1:T) {
      Tg[,t,g] = mvrnorm(1, mu=mclg[g,]+Lambda[[g]]%*%F[[g]][,t], Sigma=Sigma1[,,g])
    }}
  # sample observation from different clusters to form the mixture according to weights
  {
    z <- sample(1:G, T, replace = TRUE, prob = w)
    Y = matrix(NA, nrow=P, ncol=T)
    for ( i in 1:T) {
      Y[,i] = Tg[,i,z[i]]
    }
  }
  
  p=P
  Nz = tabulate(z,G)
  print(Nz)
  
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
}


#Identification data set with block Lambda
{
  library(MASS)
  library(mvtnorm)
  set.seed(123)
  # simulated data set K = 3, H = 3
  # cluster weights
  w = c(1/3, 1/3, 1/3)
  #number of clusters, factors and variables
  G = 3; J=4; P=20
  # number of observations
  T=300
  # factors
  F <- t(mvrnorm(T, mu=rep(0,J), Sigma=diag(1,J)))
  # loading matrices
  Lambda=array(NA, dim=c(P,J,G))
  for (g in 1:G) {
    for (i in 1:P/2) { 
      for (h in 1:2) {
        # sample bi
        indic=c(1,0)
        probs=c(0.2, 1-0.2)
        bi=sample(indic, 1, prob=probs, replace=T)
        #sample factor loadings
        Lambda[i,h,g] <- (1 + 0.1*rnorm(1,0,1))*(-1)^bi
      }
      for (h in 3:4) {Lambda[i,h,g]=0}
    } # end loop for i until P/2
    
    for (i in (P/2+1):P) { 
      for (h in 3:4) {
        # sample bi
        indic=c(1,0)
        probs=c(0.2, 1-0.2)
        bi=sample(indic, 1, prob=probs, replace=T)
        #sample factor loadings
        Lambda[i,h,g] <- (1 + 0.1*rnorm(1,0,1))*(-1)^bi
      } 
      for (h in 1:2) {Lambda[i,h,g]=0}
    } # end loop for i from P/2
  } # end loop for g
  
  # generating idiosyncratic variances
  Sigma1 = array(NA, dim=c(P,P,G))
  for (g in 1:G) {Sigma1[,,g] <- diag(x=1/rgamma(P, 2, rate=1), P)}
  # generating cluster means
  mclg <- matrix(NA, nrow = G,ncol=P)
  for (g in 1:G) {mclg[g,] = mvrnorm(1, mu=rep(2*g-G-1,P), Sigma=diag(1,P))} # mixing the means so that clusters are not very distinguishable
  # filling in clusters with T observations
  Tg <- array(NA, dim=c(P,T,G))
  #Om <- array(NA, dim=c(P,P,G))
  for ( g in 1:G) {
    for (t in 1:T) {
      Tg[,t,g] = mvrnorm(1, mu=mclg[g,]+Lambda[,,g]%*%F[,t], Sigma=Sigma1[,,g])
    }}
  # sample observation from different clusters to form the mixture according to weights
  {
    z <- sample(1:G, T, replace = TRUE, prob = w)
    Y = matrix(NA, nrow=P, ncol=T)
    for ( i in 1:T) {
      Y[,i] = Tg[,i,z[i]]
    }
    #calculate initial cluster covariance matrices
    Omega0_orig <- array(NA, dim=c(P,P,G))
    for (g in 1:G) {
      Omega0_orig[,,g]  <- Lambda[,,g]%*%t(Lambda[,,g]) + Sigma1[,,g]
    }
  }
  
  p=P
  Nz = tabulate(z,G)
  
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
  
  # calculate correlation matrix (adjust to rescaling)
  Omega0 <- array(NA, dim=c(P,P,G))
  for (g in 1:G) {
    Dm <- diag(Ys)
    Omega0[,,g]  <- solve(Dm)%*%Omega0_orig[,,g]%*%solve(Dm)
  }
  # rescale factor loading matrices
  Lambda0 <- array(NA, dim=c(P,J,G))
  for (g in 1:G) {
    Dm <- diag(Ys)
    Lambda0[,,g]  <- solve(Dm)%*%Lambda[,,g]
  }
}
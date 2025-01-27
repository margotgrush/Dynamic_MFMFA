#################################################################################################################
#                  Sampler for Dynamic Mixture of Finite Mixtures of Factor Analysers Model                     #
#################################################################################################################
# libraries
{
  library(MASS)
  library(mvtnorm)
  library(bayesm)
  library(abind)
  library(MCMCpack)
  library(mvnfast)
  library(LaplacesDemon)
}
#################################################################################################################
rm(list=ls())
#set.seed(123) 
##################################################################################################################
#                                            Variables initialisatons                                            #
##################################################################################################################
# setting the number of iterations for the Gibbs sampler chain
Nsim = 50000
# setting the initial values and priors for cluster and factor parameters
{
###############################              mixture model parameters              ###############################                       
k0=12
k=rep(0, Nsim)
k[1] <- k0 # initiating k
k_plus = rep(0, Nsim)
k_plus[1] = k0
trunc = 100 # truncation point for updating k
### prior for component weights
# Dirichlet parameter
eta <- list()
eta[[1]] <- rep(1/k0,k0) # initiate eta from Uniform distribution
# generate weights from Dirichlet
#eta[[1]] = MCMCpack::rdirichlet(n=1, alpha=rep(1/k0, k0))
# prior for alpha for Dirichlet distributions for dynamic MFM
# we use F(6,3)
a=rep(0, Nsim)
adf1=6; adf2=3
a[1] <- adf2/(adf2-2) # initiating a1
s=4.1  # step size for MH step for sampling a
### initialising S_i latent allocations of observations Y_i 
clust <- kmeans(x=t(Y), centers = k0, nstart = 30, iter.max=300, algorithm="MacQueen")
S <- list()
S[[1]] <- clust$cluster # initiate the partition
### Setting the number of counts
N <- list()
N[[1]] = rep(0, k[1])
N[[1]] = tabulate(S[[1]], k[1])
# hyperparameters for BNB distribution
a_l=1; a_pi=4; b_pi=3
rho = rep(0, Nsim) # acceptance probabilities in the MH step for Dirichlet concentration parameter
val=rep(0, Nsim) 
##############################                 factor model parameters               ###############################   
# initial number of factors
H = floor((p-1)/2) # initiating H - maximum number of factors
#H = p # initiating from maximum number of factors p
a_stick = rep(0, Nsim) # letting hyperparameter for the stick breaking process be defined by data
# generalised Beta prime distribution for a_stick ~ p(a_a, b_a, 1, H)
a_a = 6; b_a = 2*H+1 # hyperparameters for generalised beta prime
a_stick[1] <- H*a_a/(b_a-1) # mean of the prior distribution
# for adaptation of H
H_plus=list()
H_plus[[1]] = rep(H, k_plus[1])
########## idiosyncratic variances for the factor model
ag =3 ; bxi = matrix(NA, nrow=p, ncol = Nsim)
axi <- 1
bg = rep(0,p)
R0 = rep(0,p)
for (i in 1:p) {R0[i] = range(Y[i,])[2]-range(Y[i,])[1]
bg[i] <- 100/R0[i]^2}
for (i in 1:p) {bxi[i,1] <- rgamma(1, ag, rate=bg[i])}

Xi <- list()
Xi[[1]] = matrix(NA, nrow=p, ncol=k[1])
  for (i in 1:p) {
    Xi[[1]][i,]=1/rgamma(k[1],axi,rate=bxi[i,1])
  }

########## factors
FC <- list()
FC_1 <- list(1:k_plus[1])
for (j in 1:k_plus[1]) {
  FC_1[[j]] = matrix(NA, nrow=H, ncol=N[[1]][j])
  FC_1[[j]]= t(mvnfast::rmvn(N[[1]][j], mu=rep(0,H), sigma=diag(H), isChol=FALSE))
}

########## factor loadings
### hyperparameters of the slab distribution
alphath = 3; a2 = 2; b2 = 1 # hyperparameters of the slab distribution
betath = rep(0, Nsim)
betath[1] <- a2/b2
### hyperparameter of the spike
alpha_inf = 21; as = 1; bs = 1 # hyperparameters of the spike so that mean=0.05
beta_inf = rep(0, Nsim)
beta_inf[1] <- as/bs

### tau probabilities of spike and slab
tau <- list()
tau[[1]] <-matrix(NA, nrow=k[1], ncol=H)
tau[[1]][,] <- rep(a_stick[1]/(a_stick[1]+H), k[1]*H)

# initiating theta
theta <- list()
theta[[1]]=matrix(NA, nrow=k[1], ncol=H)
for (j in 1:k[1]) {
# creating a vector of slab and spike probabilities
probs = rbind(1 - tau[[1]][j,], tau[[1]][j,])
# creating a vector of slab and spike draws
slab=rep(a2/(b2*(alphath-1)), H)
spike=rep(as/(bs*(alpha_inf-1)), H)
ssl=rbind(spike, slab)  
# sample initial values of thetas
for (h in 1:(H)) {
  theta[[1]][j,h]=sample(ssl[,h], 1, prob=probs[,h], replace=TRUE)
}}

# initialing factor loadings
L <- list()
L[[1]] = array(rep(0,p*H), dim=c(p, H, k[1]))
for (j in 1:k[1]) {
    for (h in 1:H) {
      for (i in 1:p) {
      L[[1]][i,h,j]=rnorm(1,0,sd=sqrt(theta[[1]][j,h]))
    }}}

# latent binary indicator I_h
I <- list()
I[[1]] <- matrix(1, nrow=k[1], ncol=H)
###############################           component distributions parameters           ##############################           
# setting initial cluster means and variances
### component means
b0 = rep(0,p)
for (i in 1:p) {b0[i] = median(Y[i,])}
# defining range of the data in the direction i
R0 = rep(0,p)
for (i in 1:p) {R0[i] = range(Y[i,])[2]-range(Y[i,])[1]}
B0 = as.matrix(diag(R0^2))
invB0=as.matrix(diag(1/R0^2))
############## generating means
means <- list()
means[[1]] = matrix(rep(0), nrow=k[1], ncol=p)
# setting cluster centers as initial means
means[[1]] <- clust$centers

### component variances
Omega <- list()
Omega[[1]] <- array(NA, dim=c(p,p,k[1]))

### estimate Omega from SFS and Lopes (2018)
v0=3
Omega_est=matrix(0,nrow=p, ncol=p)
# Estimating Omega
sum0=0
for (t in 1:T) {
  sum0 = sum0 + Y[,t]%*%t(Y[,t])
}

# for non-standardised data
#Cor_0 = (v0 + T/2)*(solve(v0*cor(t(Y)) + 0.5*sum0)) # when working with non-standardised data we estimate correlation matrix
# scaling Omega using the correlation matrix for unstandardised data
#Omega_est = diag(diag(cov(t(Y)))^(-1/2))%*%Cor_0%*%diag(diag(cov(t(Y)))^(-1/2))
#Omega_est = chol2inv(chol(Omega_est)) # for non-standardised data

# for standardised data
Omega_est = chol2inv(chol((v0 + T/2)*(solve(v0*diag(1,p) + 0.5*sum0)))) # for standardised data, we estimate covariance matrix
# assign the same Omega to each cluster
for (j in 1:k[1]) {Omega[[1]][,,j]=Omega_est}

} # end of setting initial values 

###################################################################################################################
#                                                 Gibbs sampler                                                   #                           
###################################################################################################################
for (g in 2:Nsim) { # open Gibbs sampler
  if (g %% 1000 == 0) print(paste(g, "MCMC iterations have finished"))
  ### BLOCK 1 - update the partition C
  # step 1 (a) classify observations and determine new partition
  S[[g]] = rep(0,T)
  
  mat = sapply(1:k[g-1], function(j) eta[[g-1]][j] * mvnfast::dmvn(t(Y), mu=means[[g-1]][j,], sigma=as.matrix(Omega[[g-1]][,,j])))
  S[[g]] = apply(mat, 1, function(x) sample(1:k[g-1], 1, prob = x,replace=TRUE))
  
  # step 1 (b) rearrange other parameters so that filled clusters come first
  N[[g]] = tabulate(S[[g]], k[g-1]) # number of elements in clusters
  k_plus[g] = sum(N[[g]]!= 0) # determine the number of filled clusters K+
   
  # relabel the components so that the first K+ clusters are non-empty
  # indices of non-empty clusters
  ikp <- which(N[[g]]>0)
  ke1 <- which(N[[g]]==0)
  indices <- c(ikp, ke1)
  # rearrange parameters so that nonempty clusters are ordered first
  N[[g]] <- N[[g]][indices, drop=FALSE]
  eta[[g-1]] <- eta[[g-1]][indices, drop=FALSE]
  
  # rearrange means and variances
  means[[g-1]] <- means[[g-1]][indices,, drop=FALSE]
  Omega[[g-1]] <- Omega[[g-1]][,,indices, drop=FALSE]
  #factor model parameters
  L[[g-1]] <- L[[g-1]][,,indices, drop=FALSE]
  theta[[g-1]] <- theta[[g-1]][indices,,drop=FALSE]
  Xi[[g-1]] <- Xi[[g-1]][,indices,drop=FALSE]
  tau[[g-1]] <- tau[[g-1]][indices,,drop=FALSE]

  S1=rep(0,T)
  for (i in 1:length(indices)) {
    S1[S[[g]]==i]=which(indices==i)
  }
  S[[g]]=S1

  
  ### BLOCK 2 - Run the CUSP sampler for each nonempty cluster
  FC_2 <- list()
  L[[g]] <- array(NA, dim=c(p,H,k_plus[g]))
  Xi[[g]] <- matrix(NA,p,k_plus[g])
  Omega[[g]] <- array(NA, dim=c(p,p,k_plus[g]))
  means[[g]] <- matrix(NA, nrow=k_plus[g], ncol=p)
  I[[g]] <- matrix(NA, nrow=k_plus[g], ncol=H)
  H_plus[[g]] <- rep(NA, k_plus[g])
  tau[[g]] <-matrix(NA, nrow=k_plus[g], ncol=H)
  theta[[g]] <- matrix(NA, nrow=k_plus[g], ncol=H)
  
  for (j in 1:k_plus[g]) { #start loop for all filled clusters
    
    #Step 1 - update factors
    FC_2[[j]]=matrix(NA, nrow=H, ncol=N[[g]][j])
    vecs = matrix(NA, nrow=p, ncol=N[[g]][j])
    vecs = Y[, S[[g]]==j,drop=FALSE]
    L[[g-1]][,,j] = as.matrix(L[[g-1]][,,j])
    sigma=diag(1/Xi[[g-1]][,j])
    W=diag(H)
    DT=chol2inv(chol(W+crossprod(L[[g-1]][,,j],sigma)%*%L[[g-1]][,,j]))
    for (t in 1:N[[g]][j]) {
      dT=as.matrix(DT)%*%crossprod(L[[g-1]][,,j],sigma)%*%(vecs[,t]-means[[g-1]][j,])
      FC_2[[j]][,t]=mvnfast::rmvn(1, mu=dT, sigma=as.matrix(DT), isChol=FALSE)
    }
    
    
  #Step 2 - update factor loadings
  vecy = matrix(NA, nrow=p, ncol=N[[g]][j])
  vecy =Y[, S[[g]]==j,drop=FALSE]-means[[g-1]][j,]  
  for (i in 1:p) {
      psi <- diag(1/theta[[g-1]][j,]) # inverse of the psi matrix
      LN=chol2inv(chol(psi+(1/Xi[[g-1]][i,j])*tcrossprod(FC_2[[j]],FC_2[[j]])))  
      lN=as.matrix(LN)%*%FC_2[[j]]%*%vecy[i,]*(1/Xi[[g-1]][i,j])  
      L[[g]][i,,j]=mvnfast::rmvn(1, mu=lN, sigma=as.matrix(LN), isChol=FALSE)
  }
  
    #Step 3 - update idiosyncratic variances
    sum2=rep(0,p)
    for (i in 1:p) {
      for (t in 1:N[[g]][j]) {
        sum2[i]=sum2[i]+(vecs[i,t]-means[[g-1]][j,i]-L[[g]][i,,j]%*%FC_2[[j]][,t])^2
      }
    }
    
    bN = rep(0,p)
    for (i in 1:p) {
      aN=axi+N[[g]][j]/2
      bN[i]=bxi[i,g-1]+sum2[i]/2
      Xi[[g]][i,j]=1/rgamma(1,aN,rate=bN[i])
    }
    
    # steps for adapting of the number of factors
    # Step 4 - sample the indicators Ih
    # spike t-distribution
    pt8 <- rep(0,H)
    df8 <- 2*alpha_inf
    mean8 <- rep(0,p)
    for (h in 1:H) {
     var8 <- (beta_inf[g-1]/alpha_inf)*diag(p)
     pt8[h] <- LaplacesDemon::dmvt(x=L[[g]][,h,j], mu=mean8, S=var8, df=df8, log=FALSE)*(H/(a_stick[g-1]+H))
     }
    
    # slab t-distribution
    ptth <- rep(0,H)
    dfth <- 2*alphath
    meanth <- rep(0,p)
    for (h in 1:H) {
      varth <- (betath[g-1]/alphath)*diag(p)
      ptth[h] <- LaplacesDemon::dmvt(x=L[[g]][,h,j], mu=meanth, S=varth, df=dfth, log=FALSE)*(a_stick[g-1]/(a_stick[g-1]+H))
    }
    
   # normalising and sampling binary indicators from a categorical distribution 
    sumpt <- pt8+ptth
    pt <- rbind(pt8/sumpt, ptth/sumpt) 
    for (h in 1:H) {
      I[[g]][j,h] <- sample(0:1, 1, prob=pt[,h], replace=TRUE)
    }
    
    # Step 5 - sampling the effective number of active columns in cluster j
    H_plus[[g]][j] <- sum(I[[g]][j,])
  
    # Step 6 - update unordered slab probabilities tau
    for (h in 1:H) {
    tau[[g]][j,h] <- rbeta(1, a_stick[g-1]/H + I[[g]][j,h], 2-I[[g]][j,h])
    }
    # Step 7 - update theta given I

      for (h in 1:H) {
        if (I[[g]][j,h]==0) {
          sumpl=0
          for (i in 1:p) {sumpl = sumpl+L[[g]][i,h,j]^2}
          theta[[g]][j,h]=1/rgamma(1, alpha_inf+p/2, rate=(beta_inf[g-1] + sumpl/2))} 
        
        else {
          sump=0
          for (i in 1:p) {sump = sump+L[[g]][i,h,j]^2}
          theta[[g]][j,h]=1/rgamma(1, alphath+p/2, rate=(betath[g-1] + sump/2))}
      }
    
    # Step 8. Update cluster means and covariances
      # Step 8 (a) - update cluster covariance
      Omega[[g]][,,j] = tcrossprod(L[[g]][,,j],L[[g]][,,j]) + diag(Xi[[g]][,j])

      # Step 8 (b) - update cluster means
      bk =rep(0,p); BK=matrix(NA,p,p)
      sum1=rep(0,p);
      for (t in 1:N[[g]][j]) {sum1=sum1+(vecs[,t]-L[[g]][,,j]%*%FC_2[[j]][,t])}
      BK = chol2inv(chol(invB0+N[[g]][j]*diag(1/Xi[[g]][,j])))
      bk = as.matrix(BK)%*%(invB0%*%b0 + (diag(1/Xi[[g]][,j]))%*%sum1)
      means[[g]][j,] = mvnfast::rmvn(1, mu=bk, sigma=as.matrix(BK), isChol=FALSE)
      
  } # end loop for all filled clusters  
  
    FC[[g]] = FC_2 # create a common list for all updated factor matrices of all filled clusters
    
    # Step 9 update rate hyperparameters for spike and slab
    sumth = 0
    sumth2 = 0
    for (j in 1:k_plus[g]) {
      for (h in 1:H) {
        if (I[[g]][j,h]==0) {
          sumth = sumth + 1/theta[[g]][j,h]
        } else {sumth2 = sumth2 + 1/theta[[g]][j,h]}
      }}
    H_inf = H*k_plus[g] - sum(H_plus[[g]])
    beta_inf[g] <- rgamma(1, as + H_inf*alpha_inf, rate = (bs + sumth))
    
    H_slab = sum(H_plus[[g]])
    betath[g] <- rgamma(1, a2 + H_slab*alphath, rate = (b2 + sumth2))
    
    
  # Step 10 - update hyperparameters of the idiosyncratic variances of the factor models
  # sum of variances
    sumvar = rep(0,p)
    for (i in 1:p) {
      for (j in 1:k_plus[g]) {
        sumvar[i] = sumvar[i] + 1/Xi[[g]][i,j]
      }
    }
    # sample bxi from Gamma distribution
    for (i in 1:p) {bxi[i,g] <- rgamma(1, ag+k_plus[g]*axi, rate = (bg[i] + sumvar[i]))}
    
    # Step 12 - update ESP parameter a_stick from generalised beta prime distribution
    X <- rbeta(1, a_a+H_slab, b_a+H_inf)
    a_stick[g] <- H*X/(1-X)
  
   ### Step 3 - draw new values of k and a

    # Step 3 (a) - draw new values of K from p(K|C,a)
    trunc = 100 # max number of new values for K
    trunc1 = trunc - k_plus[g]
    # function to calculate prior for K ~ BNB(a_j. a_pi, b_pi)
    lp_k <- function(x) lgamma(a_l+x-1) + lbeta(a_l+a_pi, x-1+b_pi) - lgamma(a_l) - lgamma(x) - lbeta(a_pi, b_pi)
    # functions to calculate the rest of the posterior for K
    lrat1 <- function(x) k_plus[g]*log(a[g-1]) + lgamma(x+1) - k_plus[g]*log(x) - lgamma(x-k_plus[g]+1)
    logsumK <- function(x) {
      sumK=0
      for (jj in 1:k_plus[g]) {sumK = sumK + (lgamma(N[[g]][jj]+a[g-1]/x)-lgamma(1+a[g-1]/x))}
      return(sumK)
    }  
    
    # calculate a vector pK of posterior probabilities for K
    lpK = rep(NA, trunc)
    lpK[k_plus[g]:trunc] <- lp_k(k_plus[g]:trunc) + lrat1(k_plus[g]:trunc) + logsumK(k_plus[g]:trunc)
    
    # normalising
    max_pK = max(lpK[k_plus[g]:trunc])
    lpK[k_plus[g]:trunc] = lpK[k_plus[g]:trunc] - max_pK
    pK = exp(lpK)
    pK[1:(k_plus[g]-1)]=0
    
    # sample new K
    k[g] = sample(1:trunc,1,prob=pK, replace=TRUE)
    
# Block 3 (b)
# MH step to sample a - Dirichlet parameter for the mixture weights
# generate proposal value of a
la_p <- log(a[g-1]) + rnorm(1, 0, s)
a_p <- exp(la_p)
# function to calculate prior distribution p(a)
logpa <- function(x) df(x, df1=adf1, df2=adf2, log=TRUE)
# function to calculate posterior distribution p(a|C,K)
#first ratio
logpa1 <- function(x) lgamma(x) - lgamma(T+x) + k_plus[g]*log(x)
#product of the second ratio
logsuma2 <- function(x) {
  sumK_a=0
  for (jj in 1:k_plus[g]) {sumK_a = sumK_a + (lgamma(N[[g]][jj]+x/k[g])-lgamma(1+x/k[g]))}
  return(sumK_a)
}  
# ratio of the proposal and the previous values
lrho <- logpa(a_p) - logpa(a[g-1]) + logpa1(a_p) - logpa1(a[g-1]) + logsuma2(a_p) - logsuma2(a[g-1])
rho[g] = exp(lrho)
val[g]=(runif(1) <= min(rho[g], 1))
a[g]=a[g-1] + (a_p-a[g-1])*val[g]
 
# Block 4
# Block 4 (a) add K - K+ empty components and sample their parameters from the prior
if (k[g]>k_plus[g]) {
ke=k[g]-k_plus[g]
N[[g]]<- N[[g]][N[[g]] != 0]
N[[g]] = c(N[[g]], rep(0,ke))

# generate locations and variances from the prior
# means
means_temp = matrix(0, nrow=ke, ncol=p)
for (j in 1:ke) {
means_temp[j,]=mvnfast::rmvn(1, mu=b0, sigma=B0, isChol=FALSE)
}
means[[g]] = rbind(means[[g]], means_temp)

# idiosyncratic variances
Xi_temp = matrix(0, nrow=p, ncol=ke)
for (j in 1:ke) {
  for (i in 1:p) {
  Xi_temp[i,j]=1/rgamma(1,axi,rate=bxi[i,g])
}}
Xi[[g]] = cbind(Xi[[g]], Xi_temp)

# tau 
tau_temp <-matrix(NA, nrow=ke, ncol=H)
tau_temp[,] <- rbeta(ke*H,a_stick[g]/H,1)
tau[[g]] = rbind(tau[[g]], tau_temp)

# thetas
theta_temp = matrix(NA, nrow=ke, ncol=(H))
for (j in 1:ke) {
  # creating a vector of slab and spike probabilities
  tau1 = 1-tau_temp[j,]
  probs = rbind(tau1,tau_temp[j,]) 
  # creating a vector of slab and spike draws
  slab=1/rgamma(H, alphath, rate=betath[g])
  spike=1/rgamma(H, alpha_inf, rate=beta_inf[g])
  ssl=rbind(spike, slab)
  for (h in 1:H) {
  theta_temp[j,h] = sample(ssl[,h], 1, prob=probs[,h], replace=TRUE)
}}
theta[[g]] = rbind(theta[[g]], theta_temp)

# factor loadings
L_temp = array(NA, dim=c(p, H, ke))

for (j in 1:ke) {
  for (i in 1:p) {
    L_temp[i,,j]=mvnfast::rmvn(1,mu=rep(0,H), sigma=diag(theta_temp[j,]), isChol=FALSE)
    }}
L[[g]] = abind(L[[g]],L_temp, along=3)

# cluster variances
Omega_temp <- array(NA, dim=c(p,p,ke))
for (j in 1:ke) {
Omega_temp[,,j] = tcrossprod(L_temp[,,j],L_temp[,,j]) + diag(Xi_temp[,j])
}
Omega[[g]] = abind(Omega[[g]],Omega_temp, along=3)
} 

# Block 4 (b) sample weights from multinomial Dirichlet
  et = rep(0,k[g])
  for (j in 1:k[g]) {et[j]=a[g]/k[g]+N[[g]][j]}
  eta[[g]] <- MCMCpack::rdirichlet(n=1, alpha=et)
} # ending Gibbs sampler
###################################################################################################################
#                                                    Results                                                      #
###################################################################################################################
par(mfrow = c(1,1))
Burnin = 0.2*Nsim
# checking Dirichlet parameter a
plot(a, type="l")
print(paste("mean value of a_nu", mean(a[Burnin:Nsim])))
print(paste("acceptance rate of a_nu", sum(val)/Nsim))
# checking IBP parameter a_stick
plot(a_stick, type="l")
mean(a_stick[Burnin:Nsim])
print(paste("mean value of a_tau", mean(a_stick[Burnin:Nsim])))
# checking k and k_plus
# k
plot(k[1:Nsim], type="l")
plot(k[Burnin:Nsim], type="l")
# k_plus
plot(k_plus[1:Nsim], type="l")
plot(k_plus[Burnin:Nsim], type="l")
# thinning k
k1 <- k_plus[Burnin:Nsim]
k.new = k1[seq(1, length(k1), 5)]
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
print(paste("mode of the number of clusters is", getmode(k.new)))

####Posterior of K:
par(mfrow = c(1,1))
p_K=tabulate(k[1:Nsim],nbins=max(k[1:Nsim]))/Nsim;round(p_K,digits=2);
barplot(p_K/sum(p_K), xlab = "K", names=1:length(p_K), ylab = "p(K|y)", xlim=c(0, 30), ylim=c(0, 0.6))
#K_ is estimator of K:
K_=which.max(tabulate(k,nbins=max(k)));K_
mean(k[1:Nsim])
median(k[1:Nsim])
quantile(k[1:Nsim],probs=c(0.025,0.975))
####################################################################################################################




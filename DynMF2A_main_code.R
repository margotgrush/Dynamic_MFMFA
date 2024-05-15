#########################################################################################################
#               Sampler for Dynamic Mixture of Finite Mixtures of Factor Analysers model                #
#########################################################################################################

# libraries
{
  library(MASS)
  library(mvtnorm)
  library(bayesm)
  library(abind)
  library(MCMCpack)
  library(mclust)
  library(openxlsx)
  library(FlexDir)
}

#########################################################################################################
#                                     Gibbs sampler - initialisatons                                    #
#########################################################################################################
# setting initial conditions and prior distributions
Nsim = 50000 # number of iterations for the Gibbs sampler to run
# setting the initial values for cluster and factor parameters
{
  ################################################################################################
  #                                  mixture model parameters                                    #
  ################################################################################################
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
  eta[[1]] = MCMCpack::rdirichlet(n=1, alpha=rep(1/k0, k0))
  # prior for alpha for Dirichlet distributions for dynamic MFM (a_M in the paper)
  # we use F(6,3)
  a=rep(0, Nsim)
  adf1=6; adf2=3
  a[1] <- adf2/(adf2-2) # initiating a1
  ### initialising Si latent allocations of observations Yi 
  clust <- kmeans(t(Y), centers = k0, nstart = 30, iter.max=300, algorithm="MacQueen")
  S <- list()
  #S[[1]] <- sample(1:k0, size=T, replace=TRUE, prob=eta[[1]])
  S[[1]] <- clust$cluster
  ### Setting the Number of counts
  N <- list()
  N[[1]] = rep(0, k[1])
  N[[1]] = tabulate(S[[1]], k[1])
  # hyperparameters for BNB distribution
  a_l=1; a_pi=4; b_pi=3
  rho = rep(0, Nsim) # acceptance probabilities in the MH step
  val=rep(0, Nsim) 
  # sampling a_stick in an MH step
  ratio = rep(0, Nsim) # acceptance probabilities in the MH step
  val1=rep(0, Nsim) 
  ################################################################################################
  #                                  factor model parameters                                     #
  ################################################################################################
  # initial number of factors
  H = floor((p-1)/2) # initiating H - maximum number of factors
  #H = p # initiating from maximum number of factors p
  a_stick = rep(0, Nsim) # letting hyperparameter for the stick breaking process be defined by data (a_B in the paper)
  alpha_1 = 2; alpha_2=0.11 # hyperparameters for the MH step for a_stick
  # gamma distribution as prior for a_stick
  a_a = 6; b_a = 2 # hyperparameters for Gamma
  a_stick[1] <- a_a/b_a # mean of the prior distribution
  # for adaptation of H
  H_plus=list()
  H_plus[[1]] = rep(H, k_plus[1])
  ########## idiosyncratic variances for the factor model
  ag =3 ; bxi = matrix(NA, nrow=p, ncol = Nsim)
  axi <- 1
  bg = rep(0,p)
  R0 = rep(0,p)
  for (i in 1:p) {R0[i] = range(Y[i,])[2]-range(Y[i,])[1];
  bg[i] <- 100/R0[i]^2}
  
  for (i in 1:p) {bxi[i,1] <- rgamma(1, ag, rate=bg[i])}
  Xi <- list()
  Xi[[1]] = matrix(NA, nrow=p, ncol=k[1])
  for (j in 1:k[1]) {
    for (i in 1:p) {
      Xi[[1]][i,j]=1/rgamma(1,axi,rate=bxi[i,1])
    }}
  Xi[[1]]
  
  ########## factors
  FC <- list()
  FC_1 <- list(1:k_plus[1])
  for (j in 1:k_plus[1]) {
    FC_1[[j]]= t(mvrnorm(N[[1]][j],rep(0,H), diag(H)))
  }
  FC[[1]] = FC_1
  
  ########## factor loadings
  # hyperparameters of the slab distribution
  alphath = 3; a2 = 2; b2 = 1 # hyperparameters of the slab distribution
  betath = rep(0, Nsim)
  betath[1] <- rgamma(1, a2, rate=b2)
  # hyperparameter of the spike
  alpha_inf = 21; as = 1; bs = 1 # hyperparameters of the spike so that mean=0.05
  beta_inf = rep(0, Nsim)
  beta_inf[1] <- rgamma(1, as, rate=bs)
  
  # tau probabilities of spike and slab
  tau <- list()
  tau[[1]] <-matrix(NA, nrow=k[1], ncol=H)
  tau[[1]][,] <- rbeta(k[1]*H,a_stick[1]/H,1)
  
  # initiating theta
  theta <- list()
  theta[[1]]=matrix(NA, nrow=k[1], ncol=(H))
  for (j in 1:k[1]) {
    # creating a vector of slab and spike probabilities
    tau1 = 1 - tau[[1]][j,]
    probs = rbind(tau[[1]][j,], tau1)
    # creating a vector of slab and spike draws
    slab=rep(a2/(b2*(alphath-1)), H)
    spike=rep(as/(bs*(alpha_inf-1)), H)
    ssl=rbind(slab, spike)  
    # sample initial values of thetas
    for (h in 1:(H)) {
      theta[[1]][j,h]=sample(ssl[,h], 1, prob=probs[,h], replace=T)
    }}
  
  # initialing factor loadings
  L <- list()
  L[[1]] = array(rep(0,p*H), dim=c(p, H, k[1]))
  for (j in 1:k[1]) {
    for (i in 1:p) {
      for (h in 1:H) {
        L[[1]][i,h,j]=rnorm(1,0,sqrt(theta[[1]][j,h]))
      }}}
  
  # latent binary indicator I_h
  I <- list()
  I[[1]] <- matrix(1, nrow=k_plus[1], ncol=H)
  
  ################################################################################################
  #                           component distributions parameters                                 #
  ################################################################################################
  # setting initial cluster means and variances
  ### component means
  b0 = rep(0,p)
  for (j in 1:p) {b0[j] = median(Y[j,])}
  # defining range of the data in the direction i
  R0 = rep(0,p)
  for (i in 1:p) {R0[i] = range(Y[i,])[2]-range(Y[i,])[1]}
  B0 = diag(R0^2)
  ############## generating means
  means <- list()
  means[[1]] = matrix(rep(0), nrow=k[1], ncol=p)
  for (j in 1:k[1]) {
    means[[1]][j,] <- mvrnorm(1, mu=b0, Sigma=B0)
  }  
  # setting cluster centers as initial means
  means[[1]] <- clust$centers # when means are set as cluster centers k_plus does not change
  
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
  #Cor_0 = (v0 + T/2)*(solve(v0*cor(t(Y)) + 0.5*sum0)) # when working with non-standardised data, Cor_0 is an inverse correlation matrix
  # scaling Omega using the correlation matrix for unstandardised data
  #Omega_est = diag(diag(cov(t(Y)))^(-1/2))%*%Cor_0%*%diag(diag(cov(t(Y)))^(-1/2))
  #Omega1 = chol2inv(chol(Omega_est)) # for non-standardised data
  
  # for standardised data
  Cor_0 = (v0 + T/2)*(solve(v0*diag(1,p) + 0.5*sum0)) # when working with standardised data, Cor_0 is a prevision matrix
  Omega1 = chol2inv(chol(Cor_0)) # for standardised data
  # assign the same Omega to each cluster
  for (j in 1:k[1]) {Omega[[1]][,,j]=Omega1}
  
} # end of setting initial values 

#########################################################################################################
#                                              Gibbs sampler                                            #
#########################################################################################################
for (g in 2:Nsim) { # open Gibbs sampler
  if (g %% 1000 == 0) print(paste(g, "MCMC iterations have finished"))
  #g=2
  ### BLOCK 1 - update the partition C
  # step 1 (a) classify observations and determine new partition
  S[[g]] = rep(0,T)
  
  mat = sapply(1:k[g-1], function(j) eta[[g-1]][j] * dmvnorm(t(Y), means[[g-1]][j,], as.matrix(Omega[[g-1]][,,j])))
  S[[g]] = apply(mat, 1, function(x) sample(1:k[g-1], 1, prob = x,replace=T))
  
  # step 1 (b)
  N[[g]] = tabulate(S[[g]], k[g-1])
  k_plus[g] = sum(N[[g]]!= 0)
  
  # relabeling the components such that the first K+ clusters are non-empty
  # indices of non-empty clusters
  ikp <- which(N[[g]]!=0)
  ke1 <- which(N[[g]]==0)
  indices <- c(ikp, ke1)
  # rearrange parameters so that nonempty clusters are ordered first
  N[[g]] <- N[[g]][indices]
  eta[[g-1]] <- eta[[g-1]][indices]
  
  # rearrange means and variances
  means[[g-1]] <- means[[g-1]][indices,]
  Omega[[g-1]] <- Omega[[g-1]][,,indices]
  #factor model parameters
  L[[g-1]] <- L[[g-1]][,,indices]
  theta[[g-1]] <- theta[[g-1]][indices,]
  Xi[[g-1]] <- Xi[[g-1]][,indices]
  tau[[g-1]] <- tau[[g-1]][indices,]
  #I[[g-1]] <- I[[g-1]][indices,]
  
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
    for (t in 1:N[[g]][j]) {
      sigma=diag(Xi[[g-1]][,j])
      W=diag(H)
      DT=chol2inv(chol(W+t(L[[g-1]][,,j])%*%chol2inv(chol(sigma))%*%L[[g-1]][,,j]))
      dT=DT%*%t(L[[g-1]][,,j])%*%chol2inv(chol(sigma))%*%(vecs[,t]-means[[g-1]][j,])
      FC_2[[j]][,t]=MASS::mvrnorm(1, mu=dT, Sigma=DT)
    }
    
    #Step 2 - update factor loadings
    vecy = matrix(NA, nrow=p, ncol=N[[g]][j])
    vecy =Y[, S[[g]]==j,drop=FALSE]-means[[g-1]][j,]  
    for (i in 1:p) {
      psi <- chol2inv(chol(diag(theta[[g-1]][j,])))
      LN=chol2inv(chol(psi+(1/Xi[[g-1]][i,j])*FC_2[[j]]%*%t(FC_2[[j]])))  
      lN=(LN%*%FC_2[[j]]*(1/Xi[[g-1]][i,j]))%*%vecy[i,]  
      L[[g]][i,,j]=MASS::mvrnorm(1, mu=lN, Sigma=LN)
    }
    
    #Step 3 - update idiosyncratic variances
    sum2=rep(0,p)
    for (i in 1:p) {
      for (t in 1:N[[g]][j]) {
        sum2[i]=sum2[i]+(vecs[i,t]-means[[g-1]][j,i]-t(L[[g]][i,,j])%*%FC_2[[j]][,t])^2
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
    for (h in 1:H) {
      mean8 <- rep(0,p)
      var8 <- diag(beta_inf[g-1]/alpha_inf,p,p)
      pt8[h] <- (mvtnorm::dmvt(L[[g]][,h,j], sigma=var8, df=df8, log=FALSE)*(H/(a_stick[g-1]+H)))
    }
    
    # slab t-distribution
    ptth <- rep(0,H)
    dfth <- 2*alphath
    for (h in 1:H) {
      meanth <- rep(0,p)
      varth <- diag((betath[g-1]/alphath),p,p)
      ptth[h] <- (mvtnorm::dmvt(L[[g]][,h,j], sigma=varth, df=dfth, log=FALSE)*(a_stick[g-1]/(a_stick[g-1]+H)))
    }
    
    # sampling binary indicators from a categorical distribution 
    pt <- rbind(pt8, ptth); xx=c(0,1)
    I[[g]][j,] <-apply(pt,2, function(x) sample(xx,1,prob=x, replace=T))
    
    
    # Step 7 - sampling the effective number of active columns
    H_plus[[g]][j] <- sum(I[[g]][j,])
    
    # Step 5 - update unordered slab probabilities tau
    for (h in 1:H) {
      tau[[g]][j,h] <- rbeta(1, a_stick[g-1]/H + I[[g]][j,h], 2-I[[g]][j,h])
    }
    # Step 6 - update theta given I
    
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
    
    
    # Step 8. Update cluster means and variances
    # Step 8 (a) - update cluster variance
    Omega[[g]][,,j] = L[[g]][,,j]%*%t(L[[g]][,,j]) + diag(Xi[[g]][,j])
    
    # Step 8 (b) - update cluster means
    
    bk =rep(0,p); BN=matrix(NA,p,p)
    sum1=rep(0,p);
    for (t in 1:N[[g]][j]) {sum1=sum1+(vecs[,t]-L[[g]][,,j]%*%FC_2[[j]][,t])}
    BK = chol2inv(chol(chol2inv(chol(B0))+N[[g]][j]*chol2inv(chol(diag(Xi[[g]][,j])))))
    bk = BK%*%(chol2inv(chol(B0))%*%b0 + (chol2inv(chol(diag(Xi[[g]][,j]))))%*%sum1)
    means[[g]][j,] = MASS::mvrnorm(1, mu=bk, Sigma=BK)
    
  } # end loop for all filled clusters  
  
  FC[[g]] = FC_2 # create a common list for all updated factor matrices of all filled clusters
  
  # Step 9 (a) - sample beta_inf
  
  sumth = 0
  for (j in 1:k_plus[g]) {
    for (h in 1:H) {
      if (I[[g]][j,h]==0) {
        sumth = sumth + 1/theta[[g]][j,h]
      } else {sumth=sumth}
    }}
  H_inf = H*k_plus[g] - sum(H_plus[[g]])
  beta_inf[g] <- rgamma(1, as + H_inf*alpha_inf, rate = (bs + sumth))
  
  # Step 9 (b) - sample betath
  sumth2 = 0
  for (j in 1:k_plus[g]) {
    for (h in 1:H) {
      if (I[[g]][j,h]==1) {
        sumth2 = sumth2 + 1/theta[[g]][j,h]
      } else {sumth2=sumth2}
    }}
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
  
  # Step 11 - MH step to sample a_stick
  s_a=1+alpha_1*(1-alpha_2)^H
  # generate proposal value for a_stick
  labp <- log(a_stick[g-1]) + rnorm(1,0,s_a)
  abp <- exp(labp)
  # function to calculate prior distribution of a_stick
  logpab <- function(x) dgamma(x, shape=a_a, rate=b_a, log=T)
  # calculation ratio of the proposal and previous values
  logratio <- H_slab*(log(abp/(abp+H)) - log(a_stick[g-1]/(a_stick[g-1]+H))) + H_inf*(log(H/(abp+H)) - log(H/(a_stick[g-1]+H))) + logpab(abp) - logpab(a_stick[g-1])
  ratio[g]=exp(logratio)
  val1[g] <- (runif(1) <= min(ratio[g],1))
  a_stick[g] <- a_stick[g-1] + (abp - a_stick[g-1])*val1[g]
  
  ### Step 3 - draw new values of k and a
  
  # Step 3 (a) - draw new values of K from p(K|C,a)
  trunc = 100 # max number of
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
  k[g] = sample(1:trunc,1,prob=pK, replace=T)
  
  # Block 3 (b)
  # MH step to sample A - Dirichlet parameter   
  s=1 # set s - standard deviation of the random walk proposal
  # generate proposal value of a
  la_p <- log(a[g-1]) + rnorm(1, 0, s)
  a_p <- exp(la_p)
  # function to calculate prior distribution p(a)
  logpa <- function(x) df(x, df1=adf1, df2=adf2, log=T)
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
      means_temp[j,]=MASS::mvrnorm(1, mu=b0, Sigma=B0)
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
      probs = rbind(tau_temp[j,], tau1) 
      # creating a vector of slab and spike draws
      slab=1/rgamma(H, alphath, rate=betath[g])
      spike=1/rgamma(H, alpha_inf, rate=beta_inf[g])
      ssl=rbind(slab, spike)
      for (h in 1:H) {
        theta_temp[j,h] = sample(ssl[,h], 1, prob=probs[,h], replace=T)
      }}
    theta[[g]] = rbind(theta[[g]], theta_temp)
    
    # factor loadings
    L_temp = array(NA, dim=c(p, H, ke))
    
    for (j in 1:ke) {
      for (i in 1:p) {
        L_temp[i,,j]=MASS::mvrnorm(1,mu=rep(0,H), Sigma=diag(theta_temp[j,]))
      }}
    L[[g]] = abind(L[[g]],L_temp)
    
    # cluster variances
    Omega_temp <- array(NA, dim=c(p,p,ke))
    for (j in 1:ke) {
      Omega_temp[,,j] = L_temp[,,j]%*%t(L_temp[,,j]) + diag(Xi_temp[,j])
    }
    Omega[[g]] = abind(Omega[[g]],Omega_temp)
    
  } 
  
  # Block 4 (b) sample weights from multinomial Dirichlet
  #gam = rep(0,k[g])
  et = rep(0,k[g])
  for (j in 1:k[g]) {et[j]=a[g]/k[g]+N[[g]][j]}
  
  eta[[g]] <- MCMCpack::rdirichlet(n=1, alpha=rep(et))
  
} # ending Gibbs sampler

#########################################################################################################
#                                                 Results                                               #
#########################################################################################################
#par(mfrow = c(1,1))
Burnin = 0.2*Nsim
k
k_plus 
# checking Dirichlet parameter a
plot(a, type="l")
mean(a[Burnin:Nsim])
# checking IBP parameter a_stick
plot(a_stick, type="l")
mean(a_stick[Burnin:Nsim])
# checking densities
plot(density(a))
plot(density(a_stick))
# checking k and k_plus
# k
plot(k[Burnin:Nsim], type="l")
mean(k[Burnin:Nsim])
# k_plus
plot(k_plus[Burnin:Nsim], type="l")
mean(k_plus[Burnin:Nsim])
k1 <- k_plus[Burnin:Nsim]
# thinning k (optional)
k.new = k1[seq(1, length(k1), 5)]
plot(k.new, type="l")
mean(k.new)
sd(k.new)
#mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
print(getmode(k_plus[Burnin:Nsim]))
print(getmode(k.new))

# checking convergence
plot(cumsum(k[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
plot(cumsum(k_plus[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
plot(cumsum(a[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
plot(cumsum(a_stick[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
plot(cumsum(beta_inf[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
plot(cumsum(betath[1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")
for (i in 1:p) {plot(cumsum(bxi[i,1:Nsim])/(1:Nsim), lwd=2, ty="l", xlim=c(0,Nsim), xlab="Iterations", ylab="")}
######################################################################################################################

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




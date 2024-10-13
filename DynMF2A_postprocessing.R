#########################################################################################################
#                                   Solving label switching problem                                     #
#########################################################################################################
# to run after DynMF2A_main_code.R
#libraries
{
library(mclust)
library(infinitefactor)
}
# remove matrices from the previous postprocessing run (optional)
{
remove(L_0, Xi_0, Omega_0, means_0, I_0, tau_0, H_plus_0, eta_0, N_0, S_0, theta_0, F_0,
       ALLMEANS, ALLMEANS1, ALLMEANSld, MLDET, F_reord, L_reord, Xi_reord,
       Omega_reord, means_reord, I_reord, tau_reord, H_plus_reord, eta_reord, N_reord, S_reord, theta_reord, F_r,
       MEANS, NN, ETA, SN, LLIST, FLIST, Xi_reord1, I_reord1, ITER, z_reord, ordering,
       tau_reord1, theta_reord1, N_reord1, LLIST_ord, FLIST_ord, LLIST_act, FLIST_act, OMEG, Omega_fin, MSE,
       I_reord1_ord, tau_reord1_ord, theta_reord1_ord)
}
# set burnin for the MCMC output
Burnin = 0.3*Nsim

########################################   Rearrange clusters   ##########################################
{ 
#mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# thinning - every 5th k
k1 <- k_plus[Burnin:Nsim]
k.new = k1[seq(1, length(k1), 5)]

# determine the mode of k_plus
#Mod <- getmode(k_plus[Burnin:Nsim])
Mod <- getmode(k.new)
print(Mod)

####################################       Step 1       ##############################
# choose only those draws where k_plus = Mod
L_0 <- list()
F_0 <- list()
Xi_0 <- list()
Omega_0 <- list()
means_0 <- list()
I_0 <- list()
H_plus_0 <- list()
eta_0 <- list()
N_0 <- list()
S_0 <- list()
theta_0 <- list()
tau_0 <- list()
gg=0
for (g in Burnin:Nsim) {
  if (k_plus[g]==Mod) {
    gg = gg+1
    # initialise vector/matrices/arrays
    F_1 <- list(1:Mod)
    for (j in 1:Mod) {F_1[[j]] = matrix(NA, nrow=H, ncol=N[[g]][j])}
    L_0[[gg]] = array(NA, dim = c(p,H,Mod))
    Xi_0[[gg]] = matrix(NA, nrow=p, ncol=Mod)
    Omega_0[[gg]] = array(NA, dim = c(p,p,Mod))
    means_0[[gg]] = matrix(NA, nrow=Mod, ncol=p)
    I_0[[gg]] = matrix(NA, nrow=Mod, ncol=H)
    tau_0[[gg]] = matrix(NA, nrow=Mod, ncol=H)
    theta_0[[gg]] = matrix(NA, nrow=Mod, ncol=H)
    H_plus_0[[gg]] = rep(NA, Mod)
    eta_0[[gg]] = rep(NA, Mod)
    N_0[[gg]] = rep(NA, Mod)
    S_0[[gg]] = rep(NA, T)
    
    # assign parameters where k_plus==mode(k_plus) 
    for (j in 1:Mod) {F_1[[j]] = FC[[g]][[j]]}
    F_0[[gg]] = F_1
    L_0[[gg]] = L[[g]][,,1:Mod]
    Xi_0[[gg]] = Xi[[g]][,1:Mod]
    Omega_0[[gg]] = Omega[[g]][,,1:Mod]
    means_0[[gg]] = means[[g]][1:Mod,]
    I_0[[gg]] = I[[g]][1:Mod,]
    tau_0[[gg]] = tau[[g]][1:Mod,]
    theta_0[[gg]] = theta[[g]][1:Mod,]
    H_plus_0[[gg]] = H_plus[[g]][1:Mod]
    eta_0[[gg]] = eta[[g]][1:Mod]
    N_0[[gg]] = N[[g]][1:Mod]
    S_0[[gg]] = S[[g]]}
}
# number of iterations for which k_plus = = mode(k_plus)
NM <- gg
print(c(NM, "iterations with mode of K+", Mod))

# check number of factors in each cluster 
KH <- matrix(unlist(H_plus_0) , nrow = NM, ncol=Mod, byrow=TRUE)
KHmode = apply(KH, 2, getmode)
print(c("number of factors in clusters", KHmode))

###################################        Step 2       ######################################
# create array of NM matrices with mode(K_plus) rows and p columns consisting of means
ALLMEANS <- array(unlist(means_0), dim=c(Mod,p,NM))
# create matrix
ALLMEANS1 <- apply(ALLMEANS, 2, c)

### as k-means clustering is not functionning well, add some information from Omega
# calculate determinant of Omega and log(det) of Omega
det_Omega <- matrix(NA, ncol=Mod, nrow=NM)
ALLMEANSld <- array(NA, dim=c(Mod, p+1, NM))

for (g in 1:NM) {
  for (k in 1:Mod) {
    det_Omega[g,k] <- base::det(Omega_0[[g]][,,k])
  }
ALLMEANSld[,,g] <- cbind(ALLMEANS[,,g], log(det_Omega[g,]))
}                      
# creating an KNM pool of all draws
MLDET <- apply(ALLMEANSld, 2, c)

# apply k-means clustering for mu + log(det)
{
mu_clustld = kmeans(MLDET, Mod, nstart = 30)
mu_clustld$cluster
# create classification matrix
kclassld <- matrix(mu_clustld$cluster, nrow=NM, ncol=Mod, byrow=TRUE)

# identifying non-permutations
mld <- c() # vector of row indices which are not permutations of 1:K+
for (i in 1:NM) {
  if (length(kclassld[i,]) != length(unique(kclassld[i,])))
  {mld <- c(mld,i)}
}
}

# calculate rate of non-permutations
non_perm_rate_mld = length(mld)/NM
print(c("non-perm rate", "mu+log(det)", non_perm_rate_mld))

### get unique labeling by reordering the draws
# initialise the reordered matrices
F_reord <- list()
L_reord <- list()
Xi_reord <- list()
Omega_reord <- list()
means_reord <- list()
I_reord <- list()
tau_reord <- list()
H_plus_reord <- list()
eta_reord <- list()
N_reord <- list()
S_reord <- list()
theta_reord <- list()
F_r <- list()

# for mu+log(det)
arr1=1:Mod
for (i in 1:NM) {
    # initialise vector/matrices/arrays
    F_r <- list(1:Mod)
    L_reord[[i]] = array(NA, dim = c(p,H,Mod))
    Xi_reord[[i]] = matrix(NA, nrow=p, ncol=Mod)
    Omega_reord[[i]] = array(NA, dim = c(p,p,Mod))
    means_reord[[i]] = matrix(NA, nrow=Mod, ncol=p)
    I_reord[[i]] = matrix(NA, nrow=Mod, ncol=H)
    tau_reord[[i]] = matrix(NA, nrow=Mod, ncol=H)
    theta_reord[[i]] = matrix(NA, nrow=Mod, ncol=H)
    H_plus_reord[[i]] = rep(NA, Mod)
    eta_reord[[i]] = rep(NA, Mod)
    N_reord[[i]] = rep(NA, Mod)
    S_reord[[i]] = rep(NA, T)
    
    # assign parameters where k_plus==mode(k_plus) 
    x1=arr1[order(kclassld[i,])]
    for (j in 1:Mod) {
      x=x1[j]
      F_r[[j]] = F_0[[i]][[x]]
    }
    F_reord[[i]] = F_r
    L_reord[[i]] = L_0[[i]][,,arr1[order(kclassld[i,])]]
    Xi_reord[[i]] = Xi_0[[i]][,arr1[order(kclassld[i,])]]
    Omega_reord[[i]] = Omega_0[[i]][,,arr1[order(kclassld[i,])]]
    means_reord[[i]] = means_0[[i]][arr1[order(kclassld[i,])],]
    I_reord[[i]] = I_0[[i]][arr1[order(kclassld[i,])],]
    tau_reord[[i]] = tau_0[[i]][arr1[order(kclassld[i,])],]
    theta_reord[[i]] = theta_0[[i]][arr1[order(kclassld[i,])],]
    H_plus_reord[[i]] = H_plus_0[[i]][arr1[order(kclassld[i,])]]
    eta_reord[[i]] = eta_0[[i]][arr1[order(kclassld[i,])]]
    N_reord[[i]] = N_0[[i]][arr1[order(kclassld[i,])]]
    # Rearrange latent allocations S[[i]]
    S_reord[[i]] <- kclassld[i,][S_0[[i]]]
} # end loop for reordering

### removing the draws which are non-permutations
# if using mu + log(det)
if (!is.null(mld)) {
F_reord = F_reord[-mld]
L_reord = L_reord[-mld]
Xi_reord = Xi_reord[-mld]
Omega_reord = Omega_reord[-mld]
means_reord = means_reord[-mld]
I_reord = I_reord[-mld]
tau_reord = tau_reord[-mld]
theta_reord = theta_reord[-mld]
H_plus_reord = H_plus_reord[-mld]
eta_reord = eta_reord[-mld]
N_reord = N_reord[-mld]
S_reord = S_reord[-mld]}
}
#######################################    Calculate results     #########################################
{
# calculate cluster means
NPERM <- length(means_reord)
MEANS <- array(unlist(means_reord), dim=c(Mod,p,NPERM))
MMEANS <- apply(MEANS,c(1:2),mean)

# number of observations in a cluster N
NN <- matrix(unlist(N_reord), ncol=Mod, byrow=TRUE)
NN1 = apply(NN, 2, getmode)

# calculate cluster weights ETA
ETA <- matrix(unlist(eta_reord), ncol=Mod, byrow=TRUE)
ETAM = apply(ETA, 2, mean)

# calculate partition
SN <- matrix(unlist(S_reord), ncol=T, byrow=TRUE)
SN1 = apply(SN, 2, getmode) # predicted partition
tabulate(SN1)

# if the data is simulated
if (exists("z")==TRUE) {
  print("confusion matrix")
  table(z,SN1)
  print(c("simulated weights",ETAM, "real weights", w))
  ARI = adjustedRandIndex(SN1,z)
  print(paste("ARI=", ARI, sep=""))
  # classification error
  Err = classError(SN1,z)
  print(paste("Error rate=", Err$errorRate, sep=""))

  ###  align the clustering order with the order of the simulated clusters
# relabelling partition and get the new permutation order

if (ARI==1) {
ordering <- SN1
if (sum(Y[1,z]==Y[1,ordering])!=T) {
  TAB <-table(z,ordering)
  T1 <- dim(TAB)[1]
  #T2 <- dim(TAB)[2]
  z_reord=rep(0,Mod)
  for ( i in 1:T1) {
      z_reord[i] <- which(TAB[i,]!=0)
  } 
  
# print new confusion matrix
ordering <- z_reord[SN1]
SN1 <- ordering
table(z, SN1)
} else {z_reord=1:Mod}}
else {z_reord=1:Mod}
} else {z_reord=1:Mod}

# check labelling for benchmark datasets
{
if (exists("part1")==TRUE) {
  print("confusion matrix")
  #confusion matrix
  table(part1, SN1)
# check ARI
ARI = adjustedRandIndex(SN1,part1)
print(paste("ARI=", ARI, sep=""))
# classification error
Err = classError(SN1,part1)
print(paste("Error rate=", Err$errorRate, sep=""))}

if (exists("part2")==TRUE) {
  print("confusion matrix")
  #confusion matrix
  table(part2, SN1)
  # check ARI
ARI = adjustedRandIndex(SN1,part2)
print(paste("ARI=", ARI, sep=""))
# classification error
Err = classError(SN1,part2)
print(paste("Error rate=", Err$errorRate, sep=""))}
}
}
#################################   Cluster-specific number of factors   ####################################
{
# calculate number of factors in each cluster
KH1 <- matrix(unlist(H_plus_reord) , ncol=Mod, byrow=TRUE)
KHmode1 = apply(KH1, 2, getmode)[z_reord]
print(c("Mode of the number of factors in clusters", KHmode1))

# check the distribution of the number of factors in each cluster
KH1 <- KH1[,z_reord]
par(mfrow=c(2,3))
for (j in 1:Mod) {
nam=tabulate(KH1[,j],nbins=max(KH1[,j]))/NPERM;round(nam,digits=2);
names(nam) <- paste("p_KH", j, sep = "")
barplot(nam/sum(nam), xlab = "H+", names=1:length(nam), ylab = "p(H+|y)", main = paste("Cluster", j, sep = " "))
}

# quantiles 95%
for (j in 1:Mod) {
  print(paste("Cluster", j, sep = " "))
print(quantile(KH1[,j],probs=c(0.025,0.5,0.975)))
}
remove(L_0, Xi_0, Omega_0, means_0, I_0, tau_0, H_plus_0, eta_0, N_0, S_0, theta_0, F_0)
}
################################# Postprocessing of cluster-specific parameters #################################
{
###################################################################################################################
###                       Choose only iterations with the correct number of factors                             ###
###################################################################################################################
{
LLIST <- list()
FLIST <- list()
for (j in 1:Mod) {
lname <- paste("L_list", j, sep = "")
lvalue <- list()
LLIST[[j]] <- assign(lname, lvalue)
fname <- paste("F_list", j, sep = "")
fvalue <- list()
FLIST[[j]] <- assign(fname, fvalue)
}

Xi_reord1 <- list()
I_reord1 <- list()
tau_reord1 <- list()
theta_reord1 <- list()
N_reord1 <- list()
gg=0
for (g in 1:NPERM) { 
  if (sum(KH1[g,1:Mod]==KHmode1)==Mod) { 
    gg=gg+1 
    
    for (j in 1:Mod) {
      LLIST[[j]][[gg]] <- L_reord[[g]][,,j]
      FLIST[[j]][[gg]] <- t(F_reord[[g]][[j]])
    }
    
    Xi_reord1[[gg]] <- Xi_reord[[g]]
    I_reord1[[gg]] <- I_reord[[g]]
    tau_reord1[[gg]] <- tau_reord[[g]]
    theta_reord1[[gg]] <- theta_reord[[g]]
    N_reord1[[gg]] <- N_reord[[g]]
  } 
} 

ITER=gg
}
###################################################################################################################
###                     determine the decreasing order statistics tau(1) > ... > tau(H)                         ###
###################################################################################################################
{
# initialise vector/matrices/arrays
# Lists of loadings and factors
LLIST_ord <- list()
FLIST_ord <- list()
for (j in 1:Mod) {
  lname <- paste("L_list_ord", j, sep = "")
  lvalue <- list()
  LLIST_ord[[j]] <- assign(lname, lvalue)
  fname <- paste("F_list_ord", j, sep = "")
  fvalue <- list()
  FLIST_ord[[j]] <- assign(fname, fvalue)
}

# other parameters
I_reord1_ord <- list()
tau_reord1_ord <- list()
theta_reord1_ord <- list()

for (g in 1:ITER) {
  
  for (j in 1:Mod) {
    LLIST_ord[[j]][[g]] <- matrix(NA, nrow=p, ncol=H)
    FLIST_ord[[j]][[g]] <- matrix(NA, nrow=N_reord1[[g]][j], ncol=H)
  }

  I_reord1_ord[[g]] = matrix(NA, nrow=Mod, ncol=H)
  tau_reord1_ord[[g]] = matrix(NA, nrow=Mod, ncol=H)
  theta_reord1_ord[[g]] = matrix(NA, nrow=Mod, ncol=H)

  perm = matrix(NA, nrow=Mod, ncol=H)
  for (j in 1:Mod) {
    perm[j,] <- order(tau_reord1[[g]][j,], decreasing=T)
    
    tau_reord1_ord[[g]][j,] <- tau_reord1[[g]][j,perm[j,]]
    I_reord1_ord[[g]][j,] <- I_reord1[[g]][j,perm[j,]]
    theta_reord1_ord[[g]][j,] <- theta_reord1[[g]][j,perm[j,]]
  }

  # reorder L and F
  for (j in 1:Mod) {
    LLIST_ord[[j]][[g]] <- LLIST[[j]][[g]][,perm[j,]]
    FLIST_ord[[j]][[g]] <- FLIST[[j]][[g]][,perm[j,]]
  }}
}
###################################################################################################################
###                Choose only active loadings and factors in factor loading and factor matrices                ###    
###################################################################################################################
{
LLIST_act <- list()
FLIST_act <- list()
for (j in 1:Mod) {
  lname <- paste("L_list_act", j, sep = "")
  lvalue <- list()
  LLIST_act[[j]] <- assign(lname, lvalue)
  fname <- paste("F_list_act", j, sep = "")
  fvalue <- list()
  FLIST_act[[j]] <- assign(fname, fvalue)
}

for (g in 1:ITER) {
  for (j in 1:Mod) {
    LLIST_act[[j]][[g]] <- LLIST_ord[[j]][[g]][,which(I_reord1_ord[[g]][j,]==1)]
    FLIST_act[[j]][[g]] <- FLIST_ord[[j]][[g]][,which(I_reord1_ord[[g]][j,]==1)]
  }
}

remove(ALLMEANS, ALLMEANS1, ALLMEANSld, MLDET, F_reord, L_reord, Xi_reord,
       Omega_reord, means_reord, I_reord, tau_reord, H_plus_reord, eta_reord, N_reord, S_reord, theta_reord, F_r,
       MEANS, NN, ETA, SN, LLIST, FLIST, LLIST_ord, FLIST_ord)
}
############################################################################################################
}
#################################################################################################################
###                            Calculate cluster-specific covariance matrices                                 ###
#################################################################################################################
# rearrange into NPERM arrays and Mod lists
# calculating Omega only with the active columns of Lambda
{
OMEG <- list()
for (j in 1:Mod) {OMEG[[j]] = array(NA, dim=c(p,p,ITER))}

Omega_fin <- array(NA, dim=c(p,p,Mod))

  for (j in 1:Mod) {
    for (g in 1:ITER) {
    OMEG[[j]][,,g] = tcrossprod(LLIST_act[[j]][[g]],LLIST_act[[j]][[g]]) + diag(Xi_reord1[[g]][,j])
    }
    Omega_fin[,,j] <- apply(OMEG[[j]],c(1:2),mean)
  }
Omega_fin <- Omega_fin[,,z_reord]
}
# plot heat matrices of Omegas
plotmat(Omega_fin[,,1])
plotmat(Omega_fin[,,2])
plotmat(Omega_fin[,,3])

# number of cluster-specific factors
print(c("Mode of the number of factors in clusters", KHmode1))

# partition
SN1

# original simulated covariance matrices for comparison
#plotmat(Omega0[,,1])
#plotmat(Omega0[,,2])
#plotmat(Omega0[,,3])

#########################    calculate MSE       ############################ 
MSE <- rep(0,Mod)
for (j in 1:Mod) {
for (i in 1:p) {
  for (l in i:p) {
    #MSE[j] = MSE[j] + mean((OMEG[[j]][i,l,] - Omega0[i,l,z_reord[j]])^2)/(p*(p+1)/2)
    MSE[j] = MSE[j] + mean((OMEG[[z_reord[j]]][i,l,] - Omega0[i,l,j])^2)/(p*(p+1)/2)
  }}}

# print MSE of the cluster covariance matrices
MSE <- MSE
print(MSE)
##############################################################################################################

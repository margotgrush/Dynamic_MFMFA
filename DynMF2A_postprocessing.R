#######################################################################################
#                      Solving label switching problem                                #
#######################################################################################
# to run after DynMF2A_main_code.R
# set burnin for the MCMC output
Burnin = 0.2*Nsim
#mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# thinning - every 5th k
k1 <- k_plus[Burnin:Nsim]
k.new = k1[seq(1, length(k1), 5)]

# determine the mode of k_plus
Mod <- getmode(k_plus[Burnin:Nsim])
Modn <- getmode(k.new)
print(Mod)

#####################################       Step 1       ############################################
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

ALLMEANSd <- array(NA, dim=c(Mod, p+1, NM))
ALLMEANSld <- array(NA, dim=c(Mod, p+1, NM))

for (g in 1:NM) {
  for (k in 1:Mod) {
    det_Omega[g,k] <- base::det(Omega_0[[g]][,,k])
  }
ALLMEANSd[,,g] <- cbind(ALLMEANS[,,g], det_Omega[g,])
ALLMEANSld[,,g] <- cbind(ALLMEANS[,,g], log(det_Omega[g,]))
}                      
# creating an KNM pool of all draws
MDET <- apply(ALLMEANSd, 2, c)
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
    #S_reord[[i]] = S_0[[i]] # i indicates iteration i=1, ..., NM
    S_reord[[i]] <- kclassld[i,][S_0[[i]]]
    #for (t in 1:T) {
    #for (j in 1:Mod) {
      #vecj=which(S_0[[i]]==j)
     # lenvecj=length(vecj)
      #S_reord[[i]]=base::replace(x=S_reord[[i]], list = vecj, values=rep(kclassld[i,j], lenvecj))
    
    #check that the inverse permutation is really the inverse
    #H_plus_reord[[i]] = H_plus_0[[i]][arr1[order(kclass[i,])]]
    #(H_plus_reord[[i]][kclassld[i,]]==H_plus_0[[i]])
      
      #if (S_0[[i]][t]==j) {replace(x=S_reord[[i]][t], list=1, values=kclass[i,j])} # end loop for if... replace
   # } # end loop for j
    #} # end loop for t
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


#####################################    Step 3. Results     #######################################
# calculate cluster means
NPERM <- length(means_reord)
MEANS <- array(unlist(means_reord), dim=c(Mod,p,NPERM))
MMEANS <- apply(MEANS,c(1:2),mean)
Means1 <- apply(MEANS, 2, c)

# number of observations in a cluster N
NN <- matrix(unlist(N_reord), ncol=Mod, byrow=TRUE)
NN1 = apply(NN, 2, getmode)
sum(NN1)
zG = tabulate(z, G) # compare with the simulated cluster sizes
print(c("discovered",NN1, "real", zG))

# calculate cluster weights ETA
ETA <- matrix(unlist(eta_reord), ncol=Mod, byrow=TRUE)
ETAM = apply(ETA, 2, mean)
sum(ETAM)
print(c("simulated",ETAM, "real", w))

# calculate partition
SN <- matrix(unlist(S_reord), ncol=T, byrow=TRUE)
SN1 = apply(SN, 2, getmode) # predicted partition
tabulate(SN1)
tabulate(z)
table(z,SN1)

###  align the clustering order with the order of the simulated clusters

# relabelling partition and get the new permutation order 
{
ordering <- SN1

if (sum(Y[1,z]==Y[1,ordering])!=T) {
  TAB <-table(z,ordering)
  T1 <- dim(TAB)[1]
  T2 <- dim(TAB)[2]
  z_reord=rep(0,Mod)
  for ( i in 1:T2) {
      z_reord[i] <- which(TAB[,i]!=0)
  }
}  

# print new confusion matrix
ordering <- z_reord[SN1]
SN1 <- ordering
table(z, SN1)
}

# check labelling for simulation datasets
#library(mclust)
z # true partition
ARI = adjustedRandIndex(SN1,z)
print(ARI)
# classification error
Err = classError(SN1,z)
print(Err)

# check labelling for benchmark datasets
part1 # true partition or part1
tabulate(part1)
ARI = adjustedRandIndex(SN1,part1)
print(ARI)
# classification error
Err = classError(SN1,part1)
print(Err)

#confusion matrix
table(part1, SN1)


# check labelling for benchmark datasets
part2 # true partition or part2
tabulate(part2)
ARI = adjustedRandIndex(SN1,part2)
print(ARI)
# classification error
Err = classError(SN1,part2)
print(Err)

#confusion matrix
table(part2, SN1)

# calculate number of factors in each cluster
KH1 <- matrix(unlist(H_plus_reord) , ncol=Mod, byrow=TRUE)
KHmode1 = apply(KH1, 2, getmode)[z_reord]
#KHmode2 = apply(KH1, 2, mean)
print(c("mode number of factors in clusters", KHmode1))
#print(c("mean number of factors in clusters", KHmode2))

# check the distribution of the number of factors in each cluster
par(mfrow=c(1,1))
# 1st cluster
p_KH1=tabulate(KH1[,1],nbins=max(KH1[,1]))/NPERM;round(p_KH1,digits=2);
barplot(p_KH1/sum(p_KH1), xlab = "H+", names=1:length(p_KH1), ylab = "p(H+|y)")
# 2nd cluster
p_KH2=tabulate(KH1[,2],nbins=max(KH1[,2]))/NPERM;round(p_KH2,digits=2);
barplot(p_KH2/sum(p_KH2), xlab = "H+", names=1:length(p_KH2), ylab = "p(H+|y)")
# 3d cluster
p_KH3=tabulate(KH1[,3],nbins=max(KH1[,3]))/NPERM;round(p_KH3,digits=2);
barplot(p_KH3/sum(p_KH3), xlab = "H+", names=1:length(p_KH3), ylab = "p(H+|y)")
# 4th cluster
p_KH4=tabulate(KH1[,4],nbins=max(KH1[,4]))/NPERM;round(p_KH4,digits=2);
barplot(p_KH4/sum(p_KH4), xlab = "H+", names=1:length(p_KH4), ylab = "p(H+|y)")
# 5th cluster
p_KH5=tabulate(KH1[,5],nbins=max(KH1[,5]))/NPERM;round(p_KH5,digits=2);
barplot(p_KH5/sum(p_KH5), xlab = "H+", names=1:length(p_KH5), ylab = "p(H+|y)")
# 6th cluster
p_KH6=tabulate(KH1[,6],nbins=max(KH1[,6]))/NPERM;round(p_KH6,digits=2);
barplot(p_KH6/sum(p_KH6), xlab = "H+", names=1:length(p_KH6), ylab = "p(H+|y)")
# 7th cluster
p_KH7=tabulate(KH1[,7],nbins=max(KH1[,7]))/NPERM;round(p_KH7,digits=2);
barplot(p_KH7/sum(p_KH7), xlab = "H+", names=1:length(p_KH7), ylab = "p(H+|y)")
# 8th cluster
p_KH8=tabulate(KH1[,8],nbins=max(KH1[,8]))/NPERM;round(p_KH8,digits=2);
barplot(p_KH8/sum(p_KH8), xlab = "H+", names=1:length(p_KH8), ylab = "p(H+|y)")

#all barplots in the same plot
par(mfrow=c(2,3))
barplot(p_KH1/sum(p_KH1), xlab = expression(paste(H[1])), names=1:length(p_KH1), ylab = expression(paste("p(", H[1], "|y)")), ylim= c(0,1))
barplot(p_KH2/sum(p_KH2), xlab = expression(paste(H[2])), names=1:length(p_KH2), ylab = expression(paste("p(", H[2], "|y)")), ylim= c(0,1))
barplot(p_KH3/sum(p_KH3), xlab = expression(paste(H[3])), names=1:length(p_KH3), ylab = expression(paste("p(", H[3], "|y)")), ylim= c(0,1))
barplot(p_KH4/sum(p_KH4), xlab = expression(paste(H[4])), names=1:length(p_KH4), ylab = expression(paste("p(", H[4], "|y)")), ylim= c(0,1))
barplot(p_KH5/sum(p_KH5), xlab = expression(paste(H[5])), names=1:length(p_KH5), ylab = expression(paste("p(", H[5], "|y)")), ylim= c(0,1))
barplot(p_KH6/sum(p_KH6), xlab = expression(paste(H[6])), names=1:length(p_KH6), ylab = expression(paste("p(", H[6], "|y)")), ylim= c(0,1))

# quantiles 95%
quantile(KH1[,1],probs=c(0.025,0.5,0.975))
quantile(KH1[,2],probs=c(0.025,0.5,0.975))
quantile(KH1[,3],probs=c(0.025,0.5,0.975))
quantile(KH1[,4],probs=c(0.025,0.5,0.975))
quantile(KH1[,5],probs=c(0.025,0.5,0.975))
quantile(KH1[,6],probs=c(0.025,0.5,0.975))
quantile(KH1[,7],probs=c(0.025,0.5,0.975))
quantile(KH1[,8],probs=c(0.025,0.5,0.975))

II0 <- array(unlist(I_reord), dim=c(Mod,H,NPERM))
II10 <- apply(II0,c(1:2),getmode)

remove(L_0, Xi_0, Omega_0, means_0, I_0, tau_0, H_plus_0, eta_0, N_0, S_0, theta_0, F_0)
###################################################################################################################
###                       Choose only iterations with the correct number of factors                             ###
###################################################################################################################
L_list1 <- list()
L_list2 <- list()
L_list3 <- list()
F_list1 <- list()
F_list2 <- list()
F_list3 <- list()
Xi_reord1 <- list()
I_reord1 <- list()
tau_reord1 <- list()
theta_reord1 <- list()
N_reord1 <- list()
gg=0
for (g in 1:NPERM) { # start loop for iterations
  if (sum(KH1[g,1:Mod]==KHmode1)==Mod) { # start loop for if
    gg=gg+1   
    L_list1[[gg]] <- L_reord[[g]][,,1]
    L_list2[[gg]] <- L_reord[[g]][,,2]
    L_list3[[gg]] <- L_reord[[g]][,,3]
    F_list1[[gg]] <- t(F_reord[[g]][[1]])
    F_list2[[gg]] <- t(F_reord[[g]][[2]])
    F_list3[[gg]] <- t(F_reord[[g]][[3]])
    
    Xi_reord1[[gg]] <- Xi_reord[[g]]
    I_reord1[[gg]] <- I_reord[[g]]
    tau_reord1[[gg]] <- tau_reord[[g]]
    theta_reord1[[gg]] <- theta_reord[[g]]
    N_reord1[[gg]] <- N_reord[[g]]
  } # end loop for if
} # end loop for iterations

ITER=gg

II <- array(unlist(I_reord1), dim=c(Mod,H,NPERM))
II1 <- apply(II,c(1:2),getmode)

###################################################################################################################
###                     determine the decreasing order statistics tau(1) > ... > tau(H)                         ###
###################################################################################################################
# initialise vector/matrices/arrays
# Lists of loadings and factors
L_list1_ord <- list()
L_list2_ord <- list()
L_list3_ord <- list()
F_list1_ord <- list()
F_list2_ord <- list()
F_list3_ord <- list()
# other matrices
I_reord1_ord <- list()
tau_reord1_ord <- list()
theta_reord1_ord <- list()

for (g in 1:ITER) {
  
  L_list1_ord[[g]] <- matrix(NA, nrow=p, ncol=H)
  L_list2_ord[[g]] <- matrix(NA, nrow=p, ncol=H)
  L_list3_ord[[g]] <- matrix(NA, nrow=p, ncol=H)
  F_list1_ord[[g]] <- matrix(NA, nrow=N_reord1[[g]][1], ncol=H)
  F_list2_ord[[g]] <- matrix(NA, nrow=N_reord1[[g]][2], ncol=H)
  F_list3_ord[[g]] <- matrix(NA, nrow=N_reord1[[g]][3], ncol=H)
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
  L_list1_ord[[g]] <- L_list1[[g]][,perm[1,]]
  L_list2_ord[[g]] <- L_list2[[g]][,perm[2,]]
  L_list3_ord[[g]] <- L_list3[[g]][,perm[3,]]
  
  F_list1_ord[[g]] <- F_list1[[g]][,perm[1,]]
  F_list2_ord[[g]] <- F_list2[[g]][,perm[2,]]
  F_list3_ord[[g]] <- F_list3[[g]][,perm[3,]]
  
}


# calculate binary matrix I
I_ord <- array(unlist(I_reord1_ord), dim=c(Mod,H,NPERM))
I_ord_mode <- apply(I_ord,c(1:2),getmode)

######################################################################################################
###                         Choose only active loadings and factors                                ###
######################################################################################################
# factor loadings
L_list1_act <- list()
L_list2_act <- list()
L_list3_act <- list()
for (g in 1:ITER) {
  L_list1_act[[g]] = L_list1_ord[[g]][,which(I_reord1_ord[[g]][1,]==1)]
  L_list2_act[[g]] = L_list2_ord[[g]][,which(I_reord1_ord[[g]][2,]==1)]
  L_list3_act[[g]] = L_list3_ord[[g]][,which(I_reord1_ord[[g]][3,]==1)]
}
# factors
F_list1_act <- list()
F_list2_act <- list()
F_list3_act <- list()
for (g in 1:ITER) {
  F_list1_act[[g]] = F_list1_ord[[g]][,which(I_reord1_ord[[g]][1,]==1)]
  F_list2_act[[g]] = F_list2_ord[[g]][,which(I_reord1_ord[[g]][2,]==1)]
  F_list3_act[[g]] = F_list3_ord[[g]][,which(I_reord1_ord[[g]][3,]==1)]
}

remove(ALLMEANS, ALLMEANS1, ALLMEANSd, ALLMEANSld, MDET, MLDET, F_reord, L_reord, Xi_reord,
       Omega_reord, means_reord, I_reord, tau_reord, H_plus_reord, eta_reord, N_reord, S_reord, theta_reord, F_r,
       MEANS, Means1, NN, ETA, SN, L_list1, L_list2, L_list3, F_list1, F_list2, F_list3)
########################################################################################################
###                      Calculate Omega only with active factor loadings                            ###
########################################################################################################
# rearrange into NPERM arrays and Mod lists

# calculating Omega only with the active columns of Lambda

OM1 = array(NA, dim=c(p,p,ITER))
OM2 = array(NA, dim=c(p,p,ITER))
OM3 = array(NA, dim=c(p,p,ITER))
for (g in 1:ITER) {
  OM1[,,g] = L_list1_act[[g]]%*%t(L_list1_act[[g]]) + diag(Xi_reord1[[g]][,1])
  OM2[,,g] = L_list2_act[[g]]%*%t(L_list2_act[[g]]) + diag(Xi_reord1[[g]][,2])
  OM3[,,g] = L_list3_act[[g]]%*%t(L_list3_act[[g]]) + diag(Xi_reord1[[g]][,3])
}

Omega_m1<- apply(OM1,c(1:2),mean)
Omega_m2<- apply(OM2,c(1:2),mean)
Omega_m3<- apply(OM3,c(1:2),mean)

# plot heat matrices of Omegas
library(infinitefactor)
plotmat(Omega_m1)
plotmat(Omega_m2)
plotmat(Omega_m3)

#plotmat(Omega0[,,1])
#plotmat(Omega0[,,2])
#plotmat(Omega0[,,3])

#########################    calculate MSE       ############################ 
MSE1=0
MSE2=0
MSE3=0

for (i in 1:p) {
  for (l in i:p) {
    MSE1 = MSE1 + mean((OM1[i,l,] - Omega0[i,l,z_reord[1]])^2)/(p*(p+1)/2)
    MSE2 = MSE2 + mean((OM2[i,l,] - Omega0[i,l,z_reord[2]])^2)/(p*(p+1)/2)
    MSE3 = MSE3 + mean((OM3[i,l,] - Omega0[i,l,z_reord[3]])^2)/(p*(p+1)/2)
  }
}

# print MSE of the cluster covariance matrices
MSE <- c(MSE1, MSE2, MSE3)[z_reord]
print(MSE)
##############################################################################################################

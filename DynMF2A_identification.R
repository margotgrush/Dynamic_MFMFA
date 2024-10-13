###################################################################################################################
###                                 Identification of factor loading matrices                                   ###
###################################################################################################################

# First run the Gibbs sample with DynMF2_Amain_code.R
# Second run the postprocessing procedure DynMF2A_postprocessing.R 
# to solve label switching, choose draws with only correct (equal to cluster mode)
# number of factors, and choose only those columns in matrices of factors and factor loadings, which correspond to active factors

# libraries
{
# MatchAlign identification required the package "infinitefactor"
library(infinitefactor)
}
  
# solve rotational invariance with MatchAlign algorithm  
aligned <- list(); MA <- list()
for (j in 1:Mod) {
aligned[[j]] <- jointRot(LLIST_act[[j]], FLIST_act[[j]])
MA[[j]] <- lmean(aligned[[j]]$lambda)
}

#reorder according to the initial cluster order - at the moment only works for the same number of factors in each cluster
MA <- MA[z_reord]

# plot heat maps of the identified factor loading matrices
plotmat(MA[[1]]) 
plotmat(MA[[2]])
plotmat(MA[[3]])
plotmat(MA[[4]])
plotmat(MA[[5]])
plotmat(MA[[6]])
plotmat(MA[[7]])

# ATTENTION: check that true and estimated factor loading matrices are in the same cluster order
plotmat(Lambda0[,,1])
plotmat(Lambda0[,,2])
plotmat(Lambda0[,,3])
plotmat(Lambda0[,,4])
plotmat(Lambda0[,,5])
plotmat(Lambda0[,,6])
plotmat(Lambda0[,,7])





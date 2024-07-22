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
require(infinitefactor)
# GLT rotation requires package "ecoometric.factor.identification"
require(devtools)
require(RcppArmadillo)
#devtools::install_github("hdarjus/econometric.factor.identification")
#install.packages("econometric.factor.identification_0.1.0.tar", type = "source", repos = NULL)
#library(econometric.factor.identification)
}
  
# solve rotational invariance with MatchAlign algorithm  
aligned1 = jointRot(L_list1_act, F_list1_act)
aligned2 = jointRot(L_list2_act, F_list2_act)
aligned3 = jointRot(L_list3_act, F_list3_act)

MA1 <- lmean(aligned1$lambda)
MA2 <- lmean(aligned2$lambda)
MA3 <- lmean(aligned3$lambda)

#reorder according to the initial cluster order - at the moment only works for the same number of factors in each cluster

MA <- array(c(MA1, MA2, MA3), dim=c(p,KHmode[1],Mod))
MA <-MA[,,z_reord]
MA1 <- MA[,,1]
MA2 <- MA[,,2]
MA3 <- MA[,,3]

# plot heat maps of the identified factor loading matrices
plotmat(MA1) 
plotmat(MA2)
plotmat(MA3)

#GLT transformations of the MatchAlign identified matrices
MAGLT1 <- t(to_glt(t(MA1), tolerance=0.000001))
MAGLT2 <- t(to_glt(t(MA2), tolerance=0.000001))
MAGLT3 <- t(to_glt(t(MA3), tolerance=0.000001))
plotmat(MAGLT1)
plotmat(MAGLT2)
plotmat(MAGLT3)

# ATTENTION: check that true and estimated factor loading matrices are in the same cluster order
plotmat(Lambda0[,,1])
plotmat(Lambda0[,,2])
plotmat(Lambda0[,,3])

# real simulated loading matrices after GLT rotation
L0GLT1 <- t(to_glt(t(Lambda0[,,1]), tolerance=0.000001))
L0GLT2 <- t(to_glt(t(Lambda0[,,2]), tolerance=0.000001))
L0GLT3 <- t(to_glt(t(Lambda0[,,3]), tolerance=0.000001))
plotmat(L0GLT1)
plotmat(L0GLT2)
plotmat(L0GLT3)




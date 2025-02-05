# Dynamic MFMFA
Code for the paper "Dynamic mixture of finite mixtures of factor analysers" by M. Grushanina and S. Frühwirth-Schnatter arXiv:2307.07045, 2023

## DynMFMFA_main_code.R
--> R code to run the Gibbs sampler algorithm for the dynamic mixture of finite mixtures of factor analysers model

## DynMFMFA_convergence.R
--> R code to check the algorithm convergence

## DynMFMFA_postprocessing.R
--> R code for the post processing of the MCMC draws which resulted from running DynMFMFA_main_code.R, i.e. solving lable switching. calculating partition, ARI, error rate, number of factors in each cluster, cluster covariance matrices and their MSE

## DynMFMFA_identification.R
--> R code for sloving rotation invariance of cluster-specific factor loaing matries (uses MatchAlign algorithm available in the R package "infinitefactor" and GLT rotation procedure available in the R package "econometric.factor.identification"

## Simulations_data.R
--> R code to reproduce simulated data sets

## real_data.R
--> R code to generate benchmark data sets

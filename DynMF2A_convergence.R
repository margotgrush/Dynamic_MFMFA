#########################################################################################################
#                                           Check results                                               #
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

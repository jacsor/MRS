
library(MASS)
set.seed(1)
p = 1
n_comp = 9
n_groups = 2
n_replicates = 10
n_obs = 30
mu0 = mvrnorm(n = 1, mu = rep(0,n_comp), Sigma = 0.8*diag(1,n_comp) + 0.2)
W = array(NA, dim=c(n_comp, n_groups, n_replicates ))


for(k in 1:n_replicates)
{
  for(j in 1:n_groups)
  {
    temp = mvrnorm(n = 1, mu = mu0, Sigma = 0.1*diag(1,n_comp) )
    W[,j,k] = exp(temp) / sum(exp(temp))       
  }
}

Z = array(NA, dim=c(n_groups, n_replicates, n_obs))

X = rep(NA, n_groups * n_replicates * n_obs)
G = rep(NA, n_groups * n_replicates * n_obs)
H = rep(NA, n_groups * n_replicates * n_obs)
it = 1
for(i in 1:n_obs)
{
  for(k in 1:n_replicates)
  {
    for(j in 1:n_groups)
    {
      Z[j,k,i] = sample(1:n_comp,1, prob=W[,j,k])
      X[it] = rnorm(1, Z[j,k,i], .2)
      if(j == 1 && Z[j,k,i]==1  )
        X[it] = X[it] + 1
      G[it] = j
      H[it] = k
      it = it + 1
    }
  }  
}



library(MRS)
ans_mrs = andova(X,G,H, K=6)
plot1D(ans_mrs, type = "eff", legend = T, group = 1)
plot1D(ans_mrs, type = "eff", legend = T, group = 2)
plot1D(ans_mrs, legend = T)

# Use Group 1 as baseline for computing effect sizes
ans_mrs2 = andova(X,G,H, K=6,baseline=1)
plot1D(ans_mrs2, type = "eff", legend = T, group = 1)
plot1D(ans_mrs2, type = "eff", legend = T, group = 2)

## Draw and plot posterior samples of effect sizes and the states
n_post_samples = 100
ans_mrs = andova(X,G,H, K=6,n_post_samples = n_post_samples)
for (sample_id in 1:n_post_samples) {

  ans_mrs$RepresentativeTree = ans_mrs$PostSamples[[sample_id]]
  ans_mrs$RepresentativeTree$EffectSizes[is.nan(ans_mrs$RepresentativeTree$EffectSizes)]=0
  plot1D(ans_mrs, type = "eff", legend = T, group = 1, main =paste("Sample",sample_id))
  # plot1D(ans_mrs, legend = T, group = 1, main =paste("Sample",sample_id))
  

}



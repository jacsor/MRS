
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
ans_mrs = mrs_nested(X,G,H, K=6)

plot1D(ans_mrs, type = "eff", legend = T, group = 1)

plot1D(ans_mrs, legend = T)

library(anovaDMMT)

Omega = matrix( range(X), nrow = 1)
tau2 = seq(0.001, 100, length=100)
init_state = c(0,1)
noise = matrix(rnorm(n_groups * n_replicates * n_obs))
for( i in 1:length(tau2))
{  
  ans = fitAnovaDMMTcpp(matrix(X),G, noise, isBinary = c(0,0), J=2, K = 6, init_state=init_state, Omega=Omega, rho_star = 0.5, tau2=tau2[i])  
  post_null[i] = ans$PosteriorNull
  like[i] = exp(ans$LogLikelihood)  
}


ans_mrs$PostGlobNull
plot1D(ans_mrs)


par(mfrow=c(1,2))
sum( like*post_null / sum(like) )
plot(like)
plot(post_null, ylim = c(0,1))


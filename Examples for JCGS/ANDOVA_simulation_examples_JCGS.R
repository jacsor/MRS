setwd("~/Dropbox/work/dANOVA/Code")

library(gtools)
library(MASS)
library(MRS)
library(cramer)
library(MTSKNN)
source("BF.PT.uni.R")

nu.vec = 10^seq(-1,4,length=10)
plot.roc = function(p, p0, fpr = seq(0,1,by=0.01), col="black", xlim=c(0,1), ylim=c(0,1), main="", lty=1,lwd=2) { ## helper function to plot ROC
  plot(fpr,ecdf(p)(quantile(p0,fpr,na.rm=TRUE)), type="l", xlab="False postive rate", ylab="True positive rate", col=col, xlim=xlim, ylim=ylim, main=main, lty=lty,lwd=lwd)
}
scenario.vec = c("Null","Local shift", "Local dispersion", "Global shift", "Global dispersion")


seed = 12345
set.seed(seed)

nsim = 500 ## Number of simulations
PostGlobNull.andova.mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec))
colnames(PostGlobNull.andova.mat) = scenario.vec
PostGlobNull.2s.mat = PostGlobNull.andova.mat
stat.cramer.mat = PostGlobNull.andova.mat
pval.cramer.mat = PostGlobNull.andova.mat
pval.knn.mat = PostGlobNull.andova.mat
BF.PT.mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec));colnames(BF.PT.mat)=scenario.vec
BF.CH.mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec));colnames(BF.CH.mat)=scenario.vec

## For each of the four different scenarios, fit six different methods. The paper reports results on seven methods. The method
## not fitted and reported here is ANOVA-DDP. The result that it is not included here is that it takes substantially more time
## to fit than the other six methods and we had to fit it on a computing cluster rather than on our desktop.
for (scenario in scenario.vec) {

  print(scenario)

  p = 1
  n_comp = 3
  n_groups = 2
  n_replicates = 4
  n_obs_tot = 500
  comp.mean.null = c(1,1.5,2.5)
  sd.vec = c(0.05,0.2,rep(0.1,n_comp-2))
  
  
  for (sim in 1:nsim) {
    print(sim)

    n_obs = matrix(NA,nrow=n_groups,ncol=n_replicates)
    for (j in 1:n_groups) {  
      n_obs[j,] = rmultinom(1,size=n_obs_tot,prob=rdirichlet(1,alpha=rep(1/n_replicates,n_replicates)*n_replicates))
    }
 
    mu0=rep(0,n_comp)
    W = array(NA, dim=c(n_comp, n_groups, n_replicates ))


    for(k in 1:n_replicates)
    {
      for(j in 1:n_groups)
      {
        temp = mvrnorm(n = 1, mu = mu0, Sigma = 1*diag(1,n_comp) )
        W[,j,k] = exp(temp) / sum(exp(temp))       
      }
    }

    X = rep(NA, sum(n_obs))
    G = rep(NA, sum(n_obs))
    H = rep(NA, sum(n_obs))
    
    it = 1

    for(k in 1:n_replicates) {
      for(j in 1:n_groups) {
        for (i in 1:n_obs[j,k]) {
          Z = sample(1:n_comp,1, prob=W[,j,k])
          X[it] = rnorm(1, comp.mean.null[Z], sd.vec[Z])
 
          ## Local shifts
          if (scenario == "Local shift") {
            if(j == 1 && Z==1  )
              X[it] = X[it] + 0.1 
          }
        
          if (scenario == "Global shift") {
            if (j==1) {
              X[it] = X[it] + 0.05
            }
          }
      
          if (scenario == "Local dispersion") {
            if (j == 1 && Z==1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*3)
            }
          }
      
          if (scenario == "Global dispersion") {
            if (j == 1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*2)
            }
          }
          G[it] = j
          H[it] = k
          it = it + 1
        }
      }
    }

    ans_andova = andova(X,G,H,K=11,nu=nu.vec)
    ans_mrs = mrs(X,G,K=11)

    PostGlobNull.andova.mat[sim,scenario] = log(ans_andova$PostGlobNull)
    print(ans_andova$PostGlobNull)
    PostGlobNull.2s.mat[sim,scenario] = log(ans_mrs$PostGlobNull)
    
    stat.cramer.mat[sim,scenario] = 1 - cramer.test(X[G==1],X[G==2],rep=1)$stat
    if (scenario == "Null") {
      pval.cramer.mat[sim,scenario] = cramer.test(X[G==1],X[G==2],rep=500)$p
    }
    pval.knn.mat[sim,scenario] = mtsknn(matrix(X[G==1]),matrix(X[G==2]),k=5)$P
    BF.PT.mat[sim,scenario] = BF.PT.uni(X[G==1],X[G==2])
    BF.CH.mat[sim,scenario] = BF.CH.uni(X[G==1],X[G==2])
  }
} 


## Figure 4 in Ma and Soriano (2017)
# Plot the tree plot for one simulation run under each scenario
# In the paper, the tree plot under the local shift and local dispersion scenarios are presented
set.seed(123456789)

xgrid = seq(0.5,4.5,by=0.001)
true.den = array(0,dim=c(length(xgrid),length(scenario.vec),n_groups,n_replicates),dimnames=list(xgrid,scenario.vec,1:n_groups,1:n_replicates))

n_obs = matrix(NA,nrow=n_groups,ncol=n_replicates)
for (j in 1:n_groups) {  
  n_obs[j,] = rmultinom(1,size=n_obs_tot,prob=rdirichlet(1,alpha=rep(1/n_replicates,n_replicates)*10))
}
mu0 = rep(0,n_comp)
W = array(NA, dim=c(n_comp, n_groups, n_replicates ))

for(k in 1:n_replicates)
{
  for(j in 1:n_groups)
  {
    temp = mvrnorm(n = 1, mu = mu0, Sigma = 1*diag(1,n_comp) )
    W[,j,k] = exp(temp) / sum(exp(temp))       
  }
}

for (scenario in scenario.vec) {
  
  comp.mean.null = comp.mean.null
  comp.sd.null = sd.vec
  
  X = rep(NA, sum(n_obs))
  G = rep(NA, sum(n_obs))
  H = rep(NA, sum(n_obs))
  
  it = 1
  
  for (j in 1:n_groups) {
    
    comp.mean = comp.mean.null
    comp.sd = comp.sd.null
    
    if (scenario == "Global shift" & j==1) {
      comp.mean = comp.mean.null + 0.05
    }
    
    if (scenario == "Global dispersion" & j==1) {
      comp.sd = comp.sd.null*2
    }
    
    if (scenario == "Local shift" & j==1) {
      comp.mean[1] = comp.mean.null[1] + 0.1
    }
    
    if (scenario == "Local dispersion" & j==1) {
      comp.sd[1] = comp.sd.null[1]*3
    }
    
    for (k in 1:n_replicates) {
      for (comp in 1:n_comp) {
        true.den[,scenario,j,k] = true.den[,scenario,j,k]+W[comp,j,k] * 
          dnorm(xgrid,mean=comp.mean[comp],sd=comp.sd[comp])
      }
      
      for (i in 1:n_obs[j,k]) {
        Z = sample(1:n_comp,1, prob=W[,j,k])
        X[it] = rnorm(1, comp.mean[Z], comp.sd[Z])
        G[it] = j
        H[it] = k
        it = it + 1
      }
      
    }
    
  }
  
  ans_andova = andova(X,G,H,K=11,nu=nu.vec)
  ans_mrs = mrs(X,G,K=11)
  
  ## The following commands plot the tree plots
  plot1D(ans_andova,main=paste(scenario,": ANDOVA"),legend=TRUE)
  plot1D(ans_mrs,main=paste(scenario,": ms-BB w/ inf nu"),legend=TRUE)
  
}


## Figure 3 in Ma and Soriano (2017)
## Plot the true density along with ROC curves for each scenario
# pdf(file="Figures/ROC_curves.pdf",width=8,height = 11)
par(mfrow=c(length(scenario.vec)-1,2))
for (scenario in scenario.vec[-1]) {
  xlim=c(0.5,3)
  ylim=c(0,max(true.den[,scenario,,]))
  for (j in 1:n_groups) {
    for (k in 1:n_replicates) {
      if (k>1 | (j==2 & k==1)) { par(new=TRUE)}
      plot(x=xgrid,y=true.den[,scenario,j,k],xlab="x",ylab="Density",type="l",col=3-j,lty=j,xlim=xlim,ylim=ylim,main=scenario)
    }  
  }
  legend("topright",legend=c("Group 1","Group 2"),col=c("black","red"), lty=c(2,1),lwd=2,cex=0.8) 
  
  plot.roc(PostGlobNull.andova.mat[,scenario],PostGlobNull.andova.mat[,"Null"],main=scenario)
  par(new=TRUE)
  plot.roc(PostGlobNull.2s.mat[,scenario],PostGlobNull.2s.mat[,"Null"],col="red",lty=2)
  par(new=TRUE)
  plot.roc(stat.cramer.mat[,scenario],stat.cramer.mat[,"Null"],col="purple",lty=3)
  par(new=TRUE)
  plot.roc(pval.knn.mat[,scenario],pval.knn.mat[,"Null"],col="blue",lty=4)
  par(new=TRUE)
  plot.roc(BF.PT.mat[,scenario],BF.PT.mat[,"Null"],col="green",lty=5)
  par(new=TRUE)
  plot.roc(BF.CH.mat[,scenario],BF.CH.mat[,"Null"],col="orange",lty=6)
  abline(a=0,b=1,lty="dashed",col="gray")
  legend("bottomright",legend=c("ANDOVA","ms-BB w/ inf nu","Cramer","KNN","PT","CH"),col=c("black","red","purple","blue","green","orange"), lty=c(1,2,3,4,5,6),lwd=2,cex=0.7)
}
# dev.off()



## Figure 2 in Ma and Soriano (2017)
## Plot the true density under the null scenario as well as the histograms of the different statistic under the null

# pdf(file="Figures/null_stat_hist.pdf",width=8,height = 7)
par(mfrow=c(4,2))

scenario="Null"
xlim=c(0.5,3)
ylim=c(0,max(true.den[,scenario,,]))
for (j in 1:n_groups) {
  for (k in 1:n_replicates) {
    if (k>1 | (j==2 & k==1)) { par(new=TRUE) }
    if (k==1 & j==1) {xlab = "x";ylab="Density";main=scenario}
    else {xlab="";ylab="";main=""}
    plot(x=xgrid,y=true.den[,scenario,j,k],xlab=xlab,ylab=ylab,type="l",col=3-j,lty=j,xlim=xlim,ylim=ylim,main=main)
  }  
}
legend("topright",legend=c("Group 1","Group 2"),col=c("black","red"), lty=c(2,1),lwd=2) 

xlim=c(0,1)
hist(exp(PostGlobNull.andova.mat[,"Null"]),breaks=20,xlim=xlim,freq=TRUE,xlab="Posterior null probability",main="ANDOVA")
hist(exp(PostGlobNull.2s.mat[,"Null"]),breaks=20,xlim=c(0,1),xlab="Posterior null probability",main="ms-BB w/ inf nu")
hist(pval.cramer.mat[,"Null"],breaks=20,xlim=c(0,1),xlab="p-value",main="Cramer")
hist(pval.knn.mat[,"Null"],breaks=20,xlim=c(0,1),xlab="p-value",main="KNN")
hist(1/(1+exp(-BF.PT.mat[,"Null"])),breaks=20,xlim=xlim,xlab="Posterior null probability",main="PT",freq=TRUE)
hist(1/(1+exp(-BF.CH.mat[,"Null"])),breaks=20,xlim=xlim,xlab="Posterior null probability",main="CH",freq=TRUE)
# dev.off()


## A sensitivity analysis on the maximum depth of scanning 
K.vec = 6:12
PostGlobNull.andova.array = array(NA,dim=c(nsim,length(scenario.vec),length(K.vec)))
dimnames(PostGlobNull.andova.array)=list(1:nsim,scenario.vec,as.character(K.vec))                                  
dimnames(PostGlobNull.andova.array)                                  

for (scenario in scenario.vec) {
  
  print(scenario)
  
  p = 1
  n_comp = 3
  n_groups = 2
  n_replicates = 4
  n_obs_tot = 500
  comp.mean.null = c(1,1.5,2.5)
  sd.vec = c(0.05,0.2,rep(0.1,n_comp-2))
  
  
  for (sim in 1:nsim) {
    print(sim)
    
    n_obs = matrix(NA,nrow=n_groups,ncol=n_replicates)
    for (j in 1:n_groups) {  
      n_obs[j,] = rmultinom(1,size=n_obs_tot,prob=rdirichlet(1,alpha=rep(1/n_replicates,n_replicates)*n_replicates))
    }
    
    mu0=rep(0,n_comp)
    W = array(NA, dim=c(n_comp, n_groups, n_replicates ))
    
    
    for(k in 1:n_replicates)
    {
      for(j in 1:n_groups)
      {
        temp = mvrnorm(n = 1, mu = mu0, Sigma = 1*diag(1,n_comp) )
        W[,j,k] = exp(temp) / sum(exp(temp))       
      }
    }
    
    X = rep(NA, sum(n_obs))
    G = rep(NA, sum(n_obs))
    H = rep(NA, sum(n_obs))
    
    it = 1
    
    for(k in 1:n_replicates) {
      for(j in 1:n_groups) {
        for (i in 1:n_obs[j,k]) {
          Z = sample(1:n_comp,1, prob=W[,j,k])
          X[it] = rnorm(1, comp.mean.null[Z], sd.vec[Z])
          
          ## Local shifts
          if (scenario == "Local shift") {
            if(j == 1 && Z==1  )
              X[it] = X[it] + 0.1 
          }
          
          if (scenario == "Global shift") {
            if (j==1) {
              X[it] = X[it] + 0.05
            }
          }
          
          if (scenario == "Local dispersion") {
            if (j == 1 && Z==1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*3)
            }
          }
          
          if (scenario == "Global dispersion") {
            if (j == 1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*2)
            }
          }
          G[it] = j
          H[it] = k
          it = it + 1
        }
      }
    }
    
    for (K in K.vec) {
      ans_andova = andova(X,G,H,K=K,nu=nu.vec)
      PostGlobNull.andova.array[sim,scenario,as.character(K)] = log(ans_andova$PostGlobNull)
    }
  }
} 

# save(PostGlobNull.andova.array,file="results_sensitivity_on_K.Rdata")

## Plot the true density along with ROC curves for each scenario
# pdf(file="Figures/ROC_curves_sensitivity_on_K.pdf",width=8,height = 6)
par(mfrow=c(length(scenario.vec)/2,2))
for (scenario in scenario.vec[-1]) {

  for (i in c(1,3,5,7)) { ## The ROC curves are all essentially overlapping. We plot four of them to increase readability.
    if (i > 1) par(new=TRUE)
    plot.roc(PostGlobNull.andova.array[,scenario,as.character(K.vec[i])],PostGlobNull.andova.array[,"Null",as.character(K.vec[i])],main=scenario,col=(i+1)/2,lty=(i+1)/2)
  }
  
  abline(a=0,b=1,lty="dashed",col="gray")
  legend("bottomright",legend=paste("K=",K.vec[c(1,3,5,7)]),col=1:4, lty=1:4,lwd=2,cex=0.7)
}
# dev.off()

# pdf(file="Figures/PJAPs_sensitivity_on_K.pdf",width=8,height = 8)
par(mfrow=c(2,2))
for (scenario in scenario.vec[-1]) {
  plot(PostGlobNull.andova.array[,scenario,"6"],PostGlobNull.andova.array[,scenario,"12"],main=scenario,xlab="Posterior null probability (log), K=6",ylab="Posterior null probability (log), K=12")
  abline(a=0,b=1,lty="solid",col="gray")
}
# dev.off()

## Three group comparison
PostGlobNull.andova.3mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec))
colnames(PostGlobNull.andova.3mat) = scenario.vec
PostGlobNull.3s.3mat = PostGlobNull.andova.3mat
BF.PT.3mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec));colnames(BF.PT.3mat)=scenario.vec
BF.CH.3mat = matrix(NA,nrow=nsim,ncol=length(scenario.vec));colnames(BF.CH.3mat)=scenario.vec

## For each of the four different scenarios, fit six different methods. The paper reports results on seven methods. The method
## not fitted and reported here is ANOVA-DDP. The result that it is not included here is that it takes substantially more time
## to fit than the other six methods and we had to fit it on a computing cluster rather than on our desktop.
for (scenario in scenario.vec) {
  
  print(scenario)
  
  p = 1
  n_comp = 3
  n_groups = 3
  n_replicates = 4
  n_obs_tot = 500
  comp.mean.null = c(1,1.5,2.5)
  sd.vec = c(0.05,0.2,rep(0.1,n_comp-2))
  
  
  for (sim in 1:nsim) {
    print(sim)
    
    n_obs = matrix(NA,nrow=n_groups,ncol=n_replicates)
    for (j in 1:n_groups) {  
      n_obs[j,] = rmultinom(1,size=n_obs_tot,prob=rdirichlet(1,alpha=rep(1/n_replicates,n_replicates)*n_replicates))
    }
    
    mu0=rep(0,n_comp)
    W = array(NA, dim=c(n_comp, n_groups, n_replicates ))
    
    
    for(k in 1:n_replicates)
    {
      for(j in 1:n_groups)
      {
        temp = mvrnorm(n = 1, mu = mu0, Sigma = 1*diag(1,n_comp) )
        W[,j,k] = exp(temp) / sum(exp(temp))       
      }
    }
    
    X = rep(NA, sum(n_obs))
    G = rep(NA, sum(n_obs))
    H = rep(NA, sum(n_obs))
    
    it = 1
    
    for(k in 1:n_replicates) {
      for(j in 1:n_groups) {
        for (i in 1:n_obs[j,k]) {
          Z = sample(1:n_comp,1, prob=W[,j,k])
          X[it] = rnorm(1, comp.mean.null[Z], sd.vec[Z])
          
          ## Local shifts
          if (scenario == "Local shift") {
            if(j == 1 && Z==1  )
              X[it] = X[it] + 0.1 
          }
          
          if (scenario == "Global shift") {
            if (j==1) {
              X[it] = X[it] + 0.05
            }
          }
          
          if (scenario == "Local dispersion") {
            if (j == 1 && Z==1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*3)
            }
          }
          
          if (scenario == "Global dispersion") {
            if (j == 1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*2)
            }
          }
          G[it] = j
          H[it] = k
          it = it + 1
        }
      }
    }
    
    ans_andova = andova(X,G,H,K=11,nu=nu.vec)
    ans_mrs = mrs(X,G,K=11)
    
    PostGlobNull.andova.3mat[sim,scenario] = log(ans_andova$PostGlobNull)
    print(ans_andova$PostGlobNull)
    PostGlobNull.3s.3mat[sim,scenario] = log(ans_mrs$PostGlobNull)
    
    # stat.cramer.3mat[sim,scenario] = 1 - cramer.test(X[G==1],X[G==2],rep=1)$stat
    # if (scenario == "Null") {
    #   pval.cramer.3mat[sim,scenario] = cramer.test(X[G==1],X[G==2],rep=500)$p
    # }
    # pval.knn.3mat[sim,scenario] = mtsknn(matrix(X[G==1]),matrix(X[G==2]),k=5)$P
    BF.PT.3mat[sim,scenario] = BF.PT.uni.3sample(X[G==1],X[G==2],X[G==3])
    BF.CH.3mat[sim,scenario] = BF.CH.uni.3sample(X[G==1],X[G==2],X[G==3])
  }
} 
# save(PostGlobNull.3s.3mat,PostGlobNull.andova.3mat,BF.CH.3mat,BF.CH.3mat,true.den,file="results_3group.Rdata")

## Repeat the simulation for three group comparisons (Figure 3)
## Plot the true density along with ROC curves for each scenario
# pdf(file="Figures/ROC_curves_3group.pdf",width=8,height = 11)
par(mfrow=c(length(scenario.vec)-1,2))
for (scenario in scenario.vec[-1]) {
  xlim=c(0.5,3)
  ylim=c(0,max(true.den[,scenario,,]))
  for (j in 1:n_groups) {
    for (k in 1:n_replicates) {
      if (k>1 | (j>1 & k==1)) { par(new=TRUE)}
      plot(x=xgrid,y=true.den[,scenario,j,k],xlab="x",ylab="Density",type="l",col=4-j,lty=j,xlim=xlim,ylim=ylim,main=scenario)
    }  
  }
  legend("topright",legend=c("Group 1","Group 2","Group 3"),col=c(1,2,3), lty=c(3,2,1),lwd=2,cex=0.8) 
  
  plot.roc(PostGlobNull.andova.3mat[,scenario],PostGlobNull.andova.3mat[,"Null"],main=scenario)
  par(new=TRUE)
  plot.roc(PostGlobNull.3s.3mat[,scenario],PostGlobNull.3s.3mat[,"Null"],col="red",lty=2)
  # par(new=TRUE)
  # plot.roc(stat.cramer.3mat[,scenario],stat.cramer.3mat[,"Null"],col="purple",lty=3)
  # par(new=TRUE)
  # plot.roc(pval.knn.3mat[,scenario],pval.knn.3mat[,"Null"],col="blue",lty=4)
  par(new=TRUE)
  plot.roc(BF.PT.3mat[,scenario],BF.PT.3mat[,"Null"],col="green",lty=5)
  par(new=TRUE)
  plot.roc(BF.CH.3mat[,scenario],BF.CH.3mat[,"Null"],col="orange",lty=6)
  abline(a=0,b=1,lty="dashed",col="gray")
  legend("bottomright",legend=c("ANDOVA","ms-BB w/ inf nu","PT","CH"),col=c("black","red","green","orange"), lty=c(1,2,5,6),lwd=2,cex=0.7)
}
# dev.off()



## Repeat the simulation for three group comparison (Figure 2)
## Plot the true density under the null scenario as well as the histograms of the different statistic under the null

# pdf(file="Figures/null_stat_hist_3group.pdf",width=8,height = 7)
par(mfrow=c(3,2))

scenario="Null"
xlim=c(0.5,3)
ylim=c(0,max(true.den[,scenario,,]))
for (j in 1:n_groups) {
  for (k in 1:n_replicates) {
    if (k>1 | (j>1 & k==1)) { par(new=TRUE) }
    if (k==1 & j==1) {xlab = "x";ylab="Density";main=scenario}
    else {xlab="";ylab="";main=""}
    plot(x=xgrid,y=true.den[,scenario,j,k],xlab=xlab,ylab=ylab,type="l",col=4-j,lty=j,xlim=xlim,ylim=ylim,main=main)
  }  
}
legend("topright",legend=c("Group 1","Group 2","Group 3"),col=c(3,2,1), lty=c(1,2,3),lwd=2) 

xlim=c(0,1)
hist(exp(PostGlobNull.andova.3mat[,"Null"]),breaks=20,xlim=xlim,freq=TRUE,xlab="Posterior null probability",main="ANDOVA")
hist(exp(PostGlobNull.3s.3mat[,"Null"]),breaks=20,xlim=c(0,1),xlab="Posterior null probability",main="ms-BB w/ inf nu")
# hist(pval.cramer.mat[,"Null"],breaks=20,xlim=c(0,1),xlab="p-value",main="Cramer")
# hist(pval.knn.mat[,"Null"],breaks=20,xlim=c(0,1),xlab="p-value",main="KNN")
hist(1/(1+exp(-BF.PT.3mat[,"Null"])),breaks=20,xlim=xlim,xlab="Posterior null probability",main="PT",freq=TRUE)
hist(1/(1+exp(-BF.CH.3mat[,"Null"])),breaks=20,xlim=xlim,xlab="Posterior null probability",main="CH",freq=TRUE)
# dev.off()


## Assessing approximation accuracy for the five scenarios
n_grid_nu_prec = 100
nu_vec_prec = 10^seq(-1,4,length=n_grid_nu_prec)
n_grid_theta_prec = 100


plot.roc = function(p, p0, fpr = seq(0,1,by=0.01), col="black", xlim=c(0,1), ylim=c(0,1), main="", lty=1,lwd=2) { ## helper function to plot ROC
  plot(fpr,ecdf(p)(quantile(p0,fpr,na.rm=TRUE)), type="l", xlab="False postive rate", ylab="True positive rate", col=col, xlim=xlim, ylim=ylim, main=main, lty=lty,lwd=lwd)
}
scenario.vec = c("Null","Local shift", "Local dispersion", "Global shift", "Global dispersion")


seed = 12345
set.seed(seed)
par(mfrow=c(length(scenario.vec),2))

for (scenario in scenario.vec) {
  
  print(scenario)
  
  p = 1
  n_comp = 3
  n_groups = 2
  n_replicates = 4
  n_obs_tot = 500
  comp.mean.null = c(1,1.5,2.5)
  sd.vec = c(0.05,0.2,rep(0.1,n_comp-2))
  
    n_obs = matrix(NA,nrow=n_groups,ncol=n_replicates)
    for (j in 1:n_groups) {  
      n_obs[j,] = rmultinom(1,size=n_obs_tot,prob=rdirichlet(1,alpha=rep(1/n_replicates,n_replicates)*n_replicates))
    }
    
    mu0=rep(0,n_comp)
    W = array(NA, dim=c(n_comp, n_groups, n_replicates ))
    
    
    for(k in 1:n_replicates)
    {
      for(j in 1:n_groups)
      {
        temp = mvrnorm(n = 1, mu = mu0, Sigma = 1*diag(1,n_comp) )
        W[,j,k] = exp(temp) / sum(exp(temp))       
      }
    }
    
    X = rep(NA, sum(n_obs))
    G = rep(NA, sum(n_obs))
    H = rep(NA, sum(n_obs))
    
    it = 1
    
    for(k in 1:n_replicates) {
      for(j in 1:n_groups) {
        for (i in 1:n_obs[j,k]) {
          Z = sample(1:n_comp,1, prob=W[,j,k])
          X[it] = rnorm(1, comp.mean.null[Z], sd.vec[Z])
          
          ## Local shifts
          if (scenario == "Local shift") {
            if(j == 1 && Z==1  )
              X[it] = X[it] + 0.1 
          }
          
          if (scenario == "Global shift") {
            if (j==1) {
              X[it] = X[it] + 0.05
            }
          }
          
          if (scenario == "Local dispersion") {
            if (j == 1 && Z==1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*3)
            }
          }
          
          if (scenario == "Global dispersion") {
            if (j == 1) {
              X[it] = rnorm(1,comp.mean.null[Z],sd.vec[Z]*2)
            }
          }
          G[it] = j
          H[it] = k
          it = it + 1
        }
      }
    }
    
    ans_andova_precise = andova(X,G,H,K=K,nu=nu_vec_prec,n_grid_theta = n_grid_theta_prec,method="riemann")
    ans_andova_laplace = andova(X,G,H,K=K,nu=nu_vec_prec,n_grid_theta = n_grid_nu_prec)
    ans_andova_riemann = andova(X,G,H,K=K,nu=10^seq(-1,4,length=10),n_grid_theta = n_grid_theta_prec,method="riemann")
    
    plot(ans_andova_precise$RepresentativeTree$AltProbs,ans_andova_laplace$RepresentativeTree$AltProbs,main=scenario,xlab="Precise PMAPs",ylab="INLA PMAPs")
    abline(a=0,b=1,col="gray")
    
    plot(ans_andova_precise$RepresentativeTree$AltProbs,ans_andova_riemann$RepresentativeTree$AltProbs,main=scenario,xlab="Precise PMAPs",ylab="Riemann PMAPs with 10 grid points")
    abline(a=0,b=1,col="gray")
    
}




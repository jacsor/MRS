
post.pars <- function(M, alpha=1, y)
{ ## posterior pars for a PT prior with
  ## PT(B, A) with (right) boundaries of the partitioning sequence B
  ## defined in B[[1]]..B[[M]] for levels 1..M
  ## alpha[ep1...epsm] = alpha*m^2
  
  n <-length(y)
  neps <- list(1:M)
  a <- list(1:M)
  B <- as.list(1:M)
  for(m in 1:M){
    q <- 1/(2**m)*(1:(2**m-1))
    ## part sets (right boundaries of sets)
    B[[m]] <- 2^(-m)*seq(1,2^m)
  }
  
  for(m in 1:M)
  {
    neps[[m]] <- 0*(1:2**m)
    a[[m]] <- alpha*m**2*rep(1,2**m)
    #a[[m]] <- 0
    for(i in 1:n)
    {
      j <- which(B[[m]]>y[i])[1]
      neps[[m]][j] <- neps[[m]][j]+1
    }
    a[[m]] <- a[[m]]+neps[[m]]
  }
  return(a)
  ## initialize counts
  ## prior beta parameters
  ## update the counts
  ## updated beta pars
}



marg.log.like <- function(C, alpha=1)
{
  # C are posterior pseudo-counts
  log_like <- 0
  M <- length(C)
  for(m in 1:M)
  {
    for(j in seq(1,2^m, by=2))
    {
      log_like <- log_like + lgamma(C[[m]][j]) + lgamma(C[[m]][j+1]) + lgamma(2*alpha*m^2) - lgamma(C[[m]][j] + C[[m]][j+1]) - 2*lgamma(alpha*m^2)
      ## log_like <- log_like + lgamma(C[[m]][j]) + lgamma(C[[m]][j+1]) + lgamma(2*alpha) - lgamma(C[[m]][j] + C[[m]][j+1]) - 2*lgamma(alpha)
    }
  }
  return(log_like)
}


BF.PT.uni <- function(y1, y2, M=10, alpha=1)
{
  y.range = range(c(y1,y2))
  y1_norm = (y1 - y.range[1])/(y.range[2]-y.range[1])*0.9999
  y2_norm = (y2 - y.range[1])/(y.range[2]-y.range[1])*0.9999
  
  C1 = post.pars(M=M, y=y1_norm)
  C2 = post.pars(M=M, y=y2_norm)
  C0 = post.pars(M=M, y=c(y1_norm,y2_norm))
  
  return (marg.log.like(C0) - marg.log.like(C1) - marg.log.like(C2) )
}


BF.CH.uni <- function(y1, y2, M=10)
{
  ### null hypothesis
  y = c(y1, y2)
  sigma = sd(y)
  mu = mean(y)
  y_norm = pnorm(y, mean = mu, sd = sigma)
  
  alpha = exp( 14/19*( 1:M - 1) -7 )
  null.like = 0
  for( i in 1:length(alpha))
  {
    coeff = post.pars(M=M, y=y_norm, alpha = alpha[i])
    like = marg.log.like(coeff, alpha = alpha[i])
    if( i == 1)
      null.like = like
    else
    {
      if( like > null.like)
        null.like = like
    }
  }
  null.like = null.like + sum( dnorm(y, mean = mu, sd = sigma, log = TRUE)  )
  
  ### alternative hypothesis 
  sigma1 = sd(y1); sigma2 = sd(y2)
  mu1 = mean(y1); mu2 = mean(y2) 
  y1_norm = pnorm(y1, mean = mu1, sd = sigma1); y2_norm = pnorm(y2, mean = mu2, sd = sigma2) 
  
  alpha = exp( 14/19*( 1:M - 1) -7 )
  alt.like.1 = 0; alt.like.2 = 0
  for( i in 1:length(alpha))
  {
    coeff.1 = post.pars(M=M, y=y1_norm, alpha = alpha[i])
    like.1 = marg.log.like(coeff.1, alpha = alpha[i])
    coeff.2 = post.pars(M=M, y=y2_norm, alpha = alpha[i])
    like.2 = marg.log.like(coeff.2, alpha = alpha[i])
    if( i == 1)
    {
      alt.like.1 = like.1
      alt.like.2 = like.2
    }
    else
    {
      if( like.1 > alt.like.1)
        alt.like.1 = like.1
      if( like.2 > alt.like.2)
        alt.like.2 = like.2
    }
  }
  alt.like.1 = alt.like.1 + sum( dnorm(y1, mean = mu1, sd = sigma1, log = TRUE)  )
  alt.like.2 = alt.like.2 + sum( dnorm(y2, mean = mu2, sd = sigma2, log = TRUE)  )
  
  return (null.like - alt.like.1 - alt.like.2)  
  
}

## PT for 3 samples
BF.PT.uni.3sample <- function(y1, y2, y3, M=10, alpha=1)
{
  y.range = range(c(y1,y2,y3))
  y1_norm = (y1 - y.range[1])/(y.range[2]-y.range[1])*0.9999
  y2_norm = (y2 - y.range[1])/(y.range[2]-y.range[1])*0.9999
  y3_norm = (y3 - y.range[1])/(y.range[2]-y.range[1])*0.9999
  
  C1 = post.pars(M=M, y=y1_norm)
  C2 = post.pars(M=M, y=y2_norm)
  C3 = post.pars(M=M, y=y3_norm)
  C0 = post.pars(M=M, y=c(y1_norm,y2_norm,y3_norm))
  
  return (marg.log.like(C0) - marg.log.like(C1) - marg.log.like(C2) -marg.log.like(C3))
}



## Chen and Hanson for 3 samples
BF.CH.uni.3sample <- function(y1, y2, y3, M=10)
{
  ### null hypothesis
  y = c(y1, y2, y3)
  sigma = sd(y)
  mu = mean(y)
  y_norm = pnorm(y, mean = mu, sd = sigma)
  
  alpha = exp( 14/19*( 1:M - 1) -7 )
  null.like = 0
  for( i in 1:length(alpha))
  {
    coeff = post.pars(M=M, y=y_norm, alpha = alpha[i])
    like = marg.log.like(coeff, alpha = alpha[i])
    if( i == 1)
      null.like = like
    else
    {
      if( like > null.like)
        null.like = like
    }
  }
  null.like = null.like + sum( dnorm(y, mean = mu, sd = sigma, log = TRUE)  )
  
  ### alternative hypothesis 
  sigma1 = sd(y1); sigma2 = sd(y2); sigma3= sd(y3)
  mu1 = mean(y1); mu2 = mean(y2); mu3 = mean(y3)
  y1_norm = pnorm(y1, mean = mu1, sd = sigma1); 
  y2_norm = pnorm(y2, mean = mu2, sd = sigma2); 
  y3_norm = pnorm(y3, mean = mu3, sd = sigma3);
  
  alpha = exp( 14/19*( 1:M - 1) -7 )
  alt.like.1 = 0; alt.like.2 = 0; alt.like.3 = 0
  
  for( i in 1:length(alpha))
  {
    coeff.1 = post.pars(M=M, y=y1_norm, alpha = alpha[i])
    like.1 = marg.log.like(coeff.1, alpha = alpha[i])
    coeff.2 = post.pars(M=M, y=y2_norm, alpha = alpha[i])
    like.2 = marg.log.like(coeff.2, alpha = alpha[i])
    coeff.3 = post.pars(M=M, y=y3_norm, alpha = alpha[i])
    like.3 = marg.log.like(coeff.3, alpha = alpha[i])
    
    if( i == 1)
    {
      alt.like.1 = like.1
      alt.like.2 = like.2
      alt.like.3 = like.3
    }
    else
    {
      if( like.1 > alt.like.1)
        alt.like.1 = like.1
      if( like.2 > alt.like.2)
        alt.like.2 = like.2
      if (like.3 > alt.like.3)
        alt.like.3 = like.3
    }
  }
  alt.like.1 = alt.like.1 + sum( dnorm(y1, mean = mu1, sd = sigma1, log = TRUE)  )
  alt.like.2 = alt.like.2 + sum( dnorm(y2, mean = mu2, sd = sigma2, log = TRUE)  )
  alt.like.3 = alt.like.3 + sum( dnorm(y3, mean = mu3, sd = sigma3, log = TRUE)  )
  
  return (null.like - alt.like.1 - alt.like.2 - alt.like.3)  
  
}




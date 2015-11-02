\name{MRS-package}
\alias{MRS-package}
\alias{MRS}
\docType{package}
\title{
Multi-resolution multi-sample comparison using MRS}
\description{
This package fits the Multi Resolution Scanning algorithm (MRS) for comparison across distributions.
The model is based on a nonparametric process taking the form of a Markov model that transitions between a "null" and an "alternative" state on a multi-resolution partition tree of the sample space. The MRS effectively detects and characterizes a variety of underlying differences. These differences can be visualized using several plotting functions. }
\details{
\tabular{ll}{
Package: \tab MRS \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-23\cr
License: \tab GPL (>= 3) \cr
}

}
\author{
\tabular{ll}{
Authors: \tab Jacopo Soriano and Li Ma \cr
Maintainer: \tab Jacopo Soriano <jacopo.soriano@gmail.com> \cr
}
}
\references{
Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
 through acceptance-rejection Markov trees. \emph{Preprint}. \url{http://arxiv.org/abs/1404.3753}
}

\examples{
## EXAMPLE 1D ##

set.seed(1)
p = 1
n1 = 200
n2 = 200
mu1 = matrix( c(0,10), nrow = 2, byrow = TRUE)
mu2 = mu1; mu2[2] = mu1[2] + .01
sigma = c(1,.1)

Z1 = sample(2, n1, replace=TRUE, prob=c(0.9, 0.1))
Z2 = sample(2, n2, replace=TRUE, prob=c(0.9, 0.1))
X1 = mu1[Z1] + matrix(rnorm(n1*p), ncol=p)*sigma[Z1]
X2 = mu2[Z2] + matrix(rnorm(n2*p), ncol=p)*sigma[Z1]
X = rbind(X1, X2)
G = c(rep(1, n1), rep(2,n2))

ans = mrs(X, G, K=10)
out = summary(ans)
plot1D(ans, type = "prob")

## EXAMPLE 2D ##

set.seed(1)
p = 2
n1 = 200
n2 = 200
mu1 = matrix( c(9,9,0,4,-2,-10,3,6,6,-10), nrow = 5, byrow=TRUE)
mu2 = mu1; mu2[2,] = mu1[2,] + 1

Z1 = sample(5, n1, replace=TRUE)
Z2 = sample(5, n2, replace=TRUE)
X1 = mu1[Z1,] + matrix(rnorm(n1*p), ncol=p)
X2 = mu2[Z2,] + matrix(rnorm(n2*p), ncol=p)
X = rbind(X1, X2)
colnames(X) = c(1,2)
G = c(rep(1, n1), rep(2,n2))

ans = mrs(X, G, K=10)
plot2D(ans, type = "prob", legend = TRUE)


}
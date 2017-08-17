Probabilistic multi-resolution scanning for cross-sample variation
================================

This package fits the MRS algorithm for comparison across probability distributions. 
The model is based on a nonparametric process taking the form of a Markov model that transitions 
between a "null" and a "alternative" state on a multi-resolution partition tree of the sample space.
MRS effectively detects and characterizes a variety of underlying differences. 
These differences can be visualized using several plotting functions.

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('MaStatLab/MRS')
```

### Use
There are six functions in this package, and their descriptions are provided in the help files

```S
ans = mrs(X, G)
ans = andova(X, G, H)
summary(ans)
plot1D(ans)
plot2D(ans)
plotTree(ans)
```

### Reference
Soriano J. and Ma L. (2017). Probabilistic multi-resolution scanning for two-sample differences. Journal of the Royal Statistical Society. Series B. Vol. 79, No. 2, 547-572.

Ma L. and Soriano J. (2017). Analysis of distributional variation through multi-scale Beta-Binomial modeling.

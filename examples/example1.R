library(MRS)
set.seed(1) 
n = 20
p = 2
X = matrix(c(runif(p*n/2),rbeta(p*n/2, 1, 4)), nrow=n, byrow=TRUE)
G = c(rep(1,n/2), rep(2,n/2))
G0 = sample(2, n, replace = TRUE)


ans = mrs(X=X, G=G, K = 4)
ans$PostGlobNull
ans$LogLikelihood


ans = mrs(X=X, G=G0, K = 4)
ans$PostGlobNull
ans$LogLikelihood


ans = mrs(X=X, G=G, K = 4, return_global_null= FALSE)
ans = mrs(X=X, G=G, K = 4, return_tree= FALSE)
ans = mrs(X=X, G=G, K = 4, return_global_null= FALSE, return_tree = FALSE)


n = 1000
M = 5
class_1 = sample(M, n, prob= 1:5, replace=TRUE  )
class_2 = sample(M, n, prob = 5:1, replace=TRUE )

Y_1 = rnorm(n, mean=class_1, sd = .2)
Y_2 = rnorm(n, mean=class_2, sd = .2)

Y = matrix( c(Y_1, Y_2), ncol = 1)
G = c(rep(1,n),rep(2,n))
H = sample(3,2*n, replace = TRUE  )

ans = mrs_nested(Y, G, H)
ans$PostGlobNull
plot1D(ans)



####################


n = 1000
M = 5
H = sample(M,2*n, replace = TRUE  )

Y_1 = rnorm(n, mean=H, sd = .2)
Y_2 = rnorm(n, mean=H, sd = .2)

Y = matrix( c(Y_1, Y_2), ncol = 1)
G = c(rep(1,n),rep(2,n))

ans = mrs_nested(Y, G, H)
ans$PostGlobNull
plot1D(ans)

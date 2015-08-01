


theta0 = seq(0.01, 0.99, by = 0.01)
N = 100
data_0 = sample( seq(10, 100, by = 10), N, replace = TRUE)
data_1 = 10*data_0
nu = 1
alpha = 0.5
h = rep(NA, length(theta0))
for( i in 1:length(theta0))
  h[i] = eval_h(theta0[i], data_0, data_1, nu, alpha)
plot(theta0, h)
plot(theta0, exp(h) )


data_0[1] 
data_1[1] 
nu = 1E10
y = lbeta( nu*theta0 + data_0[1], nu*(1-theta0) + data_1[1]) - lbeta(nu*theta0, nu*(1-theta0))
plot(theta0, y, ylab = "" )
theta0[which.max(y)]
max(y)
plot(exp(theta0-max(theta0)), y, ylab = "", xlab = "" )




theta0 = seq(0.01, 0.99, by = 0.01)
N = 100
data_0 = sample( seq(10, 100, by = 10), N, replace = TRUE)
data_1 = 10*data_0 + sample( seq(-5, 5), N, replace = TRUE)
nu = exp(seq(-1, 5))
alpha = 0.5
plot(1, type="n", axes=T, xlab="theta", ylab="", ylim = c(-23000,-18500), xlim = c(0,1))
for(j in 1:length(nu))
{
  for( i in 1:length(theta0))
    points( theta0[i], eval_h(theta0[i], data_0, data_1, nu[j], alpha), col = j)
}

theta0 = seq(0.01, 0.99, by = 0.01)
N = 100
data_0 = sample( seq(10, 100, by = 10), N, replace = TRUE)
data_1 = data_0 + sample( seq(-5, 5), N, replace = TRUE)
nu = exp(seq(-1, 5))
alpha = 0.5
plot(1, type="n", axes=T, xlab="theta", ylab="", ylim = c(-10000,-7000), xlim = c(0,1))
for(j in 1:length(nu))
{
  for( i in 1:length(theta0))
    points( theta0[i], eval_h(theta0[i], data_0, data_1, nu[j], alpha), col = j)
}




theta0 = seq(0.01, 0.99, by = 0.01)
N = 2

data_0 = sample( seq(5, 10), N, replace = TRUE)
data_1 = sample( seq(5, 30), N, replace = TRUE) # data_0 + 
nu = exp(seq(-1, 5))
alpha = 0.5
plot(1, type="n", axes=T, xlab="theta", ylab="", ylim = c(-10,10), xlim = c(0,1))
output = matrix(NA, ncol = length(nu), nrow = length(theta0))
for(j in 1:length(nu))
{
  for( i in 1:length(theta0))
    output[i,j] = eval_h(theta0[i], data_0, data_1, nu[j], alpha) 
}


out = output - max(output)

plot(theta0, exp(out[,1]), col = 1, ylim = range(exp(out)), type = 'l')
for(j in 2:length(nu))
  lines(theta0, exp(out[,j]), col = j)


newtonMethod(data_0, data_1, nu = 20, alpha)
newtonMethod(data_0, data_1, nu = 148.41, alpha)

sum(data_0 + data_1)

range(exp( output[,7]/N ))
range(exp( output[,2]/N ))

par(mfrow = c(1,1))
newtonMethod(data_0, data_1, nu = nu[7], alpha)

a = -93.9148 / 2
b = - 2*a*0.336623
c = -15.6891 - a*0.336623*0.336623 - b * 0.336623
plot(theta0, output[,7]/N, col = "red")
lines(theta0, a*theta0^2 + b*theta0+c, col = "red")


newtonMethod(data_0, data_1, nu = nu[2], alpha)

a = -9.6148 / 2
b = - 2*a*0.438886
c = -16.8217 - a*0.438886*0.438886 - b * 0.438886
points(theta0, output[,2]/N)
lines(theta0, a*theta0^2 + b*theta0+c )



plot(theta0, exp(output[,7]), col = "red")
plot(theta0, dnorm(theta0, 0.336623, sqrt(-1/a)  ), col = "red")

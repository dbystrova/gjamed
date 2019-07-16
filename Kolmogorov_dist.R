
N=10000

x=rnorm(N)

max(abs(pnorm(sort(x)) - (1:N)/N ))

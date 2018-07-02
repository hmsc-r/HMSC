library(ggplot2)

N = 10000
m = 2
r = 1e9
sdEps = 1e-3

psi = m + rnorm(N,0,sdEps)
y1 = rpois(N,exp(psi))
p = exp(psi)/(r+exp(psi))


psiTrue = log(p/(1-p))
mTrue = m - log(r)
# sdPsiTrue = sd(psiTrue)
y2 = rnbinom(N, r, prob=1-p )
Y = rbind(data.frame(y=y1,m="poisson"),data.frame(y=y2,m="NB"))

ggplot(Y, aes(x=y, fill=m)) + geom_histogram(alpha=0.2, position="identity")


y = y2
w = rpg(N, y+r, psiTrue)
prec = sdEps^-2
sigmaZ = (prec + w)^-1
muZ = sigmaZ*((y-r)/2 + prec*mTrue)

# muZ
z = rnorm(length(y), muZ, sqrt(sigmaZ))

plot(psiTrue, muZ)
abline(0,1,col="red")

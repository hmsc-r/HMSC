


# create 2 random levels and specify both data and priors for them
rL1 = HmscRandomLevel$new(N=ny, priors="default")
rL2 = HmscRandomLevel$new(data=S)
rL2$setPriors(nfMax=10, mu=3, a1=5, b1=1, a2=3, b2=1, alphapw=cbind(p=c(0,1,2,3),w=c(0.4,0.2,0.2,0.2)))


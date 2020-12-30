makedata = function(ns, ny, hierarchical=FALSE, spatial=FALSE){
  rho = 0.5
  sigma = 0.3
  study.design = data.frame("sampling.unit" = factor(sprintf('su_%.3d',1:ny)),
                            "plot" = factor(sprintf('plot_%.3d',rep(1:10,ny/10))))
  xy = cbind(runif(ny), runif(ny))

  X.categorical = factor(sample(c("A","B","C"), ny, replace=TRUE))
  X.covariate = rnorm(n=ny)
  X.data = data.frame(X.categorical, X.covariate)
  X.formula = ~ X.categorical + X.covariate
  X = model.matrix(X.formula, data=X.data)
  nc = dim(X)[2]

  Tr.categorical = factor(sample(c("A","B","C"), ns, replace=TRUE))
  Tr.covariate = rnorm(n=ns)
  Tr.data = data.frame(Tr.categorical, Tr.covariate)
  Tr.formula = ~ Tr.categorical + Tr.covariate
  Tr = model.matrix(Tr.formula, data=Tr.data)
  nt = dim(Tr)[2]

  C = matrix(0,nrow=ns, ncol=ns)
  for (i in 1:ns){
    for (j in 1:ns){
      if(floor((i-1)/5)==floor((j-1)/5)){
        C[i,j] = 0.9
      }
      if(i==j){C[i,j] = 1}
    }
  }

  V = (sigma/2)^2*diag(nc)
  gamma = matrix(rnorm(n = nc*nt, mean=0, sd=sigma), ncol=nt, nrow=nc)
  mu = tcrossprod(gamma,Tr)
  Si = kronecker(V,rho*C + (1-rho)*diag(ns))
  beta = matrix(mvrnorm(n=1, mu=as.vector(mu), Sigma=Si), ncol=ns, nrow=nc)
  LF = X%*%beta

  eta = matrix(0, ncol=2, nrow=ny)
  if (spatial){
    di = as.matrix(dist(xy))
    Si.alpha = exp(-di/0.5)
    eta[,1] = mvrnorm(mu=rep(0,ny), Sigma=Si.alpha)
    Si.alpha = exp(-di/0.1)
    eta[,2] = mvrnorm(mu=rep(0,ny), Sigma=Si.alpha)
  } else {
    eta[,1] = rnorm(n=ny)
    if (hierarchical){
      tmp = rnorm(n = 10)
      eta[,2] = rep(tmp,ny/10)
    } else {
      eta[,2] = rnorm(n=ny)
    }
  }
  lambda = matrix(rnorm(n = 2*ns, mean=0, sd=sigma), ncol=ns, nrow=2)
  LR = eta%*%lambda

  L = LF + LR
  eps = matrix(rnorm(n = ny*ns, mean=0, sd=2*sigma), ncol=ns, nrow=ny)
  Y = L + eps
  all.data = list(study.design = study.design, X.data = X.data, X.formula = X.formula, Tr.data = Tr.data, Tr.formula = Tr.formula, C = C, Y = Y, xy = xy)
  all.parameters = list(gamma = gamma, beta = beta, rho = rho, V = V, eta = eta, lambda = lambda, sigma = sigma, L = L, mu = mu, LF = LF, LR = LR, eps = eps)
  list(all.data, all.parameters)
}

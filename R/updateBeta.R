updateBeta = function(Z,Gamma,iV,iSigma,Eta,Lambda, X,Tr,Pi){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   LFix = 0
   LRan = vector("list", nr)
   for(r in 1:nr){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }
   S = Z - (LFix + Reduce("+", LRan))

   Q = crossprod(X)
   Mu = tcrossprod(Gamma,Tr)
   isXTS = crossprod(X,S) * matrix(iSigma,nc,ns,byrow=TRUE)

   # iVTen = rep.tensor(to.tensor(iV), ns, pos=3)
   # QTen = rep.tensor(to.tensor(Q), ns, pos=3)

   iVTen = to.tensor(iV)
   QTen = to.tensor(Q)

   iSigmaTen = to.tensor(iSigma)
   names(iSigmaTen) = "ns"
   iUBetaTen = iVTen + QTen*iSigmaTen
   RiUBetaTen = chol.tensor(iUBetaTen, "I1", "I2")

   RiUBetaArray = to.matrix.tensor(RiUBetaTen,i="lambda",j="I1",by="ns")
   Beta = matrix(NA, nc, ns)
   for(j in 1:ns){
      # RUBeta = chol(Matrix(iV + Q/sigma[j]))
      RiUBeta = RiUBetaArray[,,j]
      UBeta = chol2inv(RiUBeta);
      mBeta = UBeta %*% (iV%*%Mu[,j] + isXTS[,j]);
      Beta[,j] = mBeta + backsolve(RiUBeta, rnorm(nc))
   }
   return(Beta)

   # # mul.tensor(RiUBetaTen, i="lambda", RiUBetaTen, by="ns")
   # mIntTen = to.tensor(iV%*%Mu + isXTS)
   # names(mIntTen)[2] = "ns"
   # randTen = to.tensor(matrix(rnorm(nc*ns),nc,ns))
   # names(randTen) = c("lambda","ns")
   # tmpTen = solve.tensor(RiUBetaTen, mIntTen, i=c("I1"), by="ns");
   # Beta = solve.tensor(RiUBetaTen, tmpTen + randTen, i="lambda", by="ns")
   # return(to.matrix.tensor(Beta,i="I1",j="ns"))
}


combineParameters = function(Beta,Gamma,iV,rho,iSigma,Eta,Lambda,Alpha,Psi,Delta, rhopw){
   RiV = chol(iV)
   V = chol2inv(RiV)
   sigma = 1/iSigma
   par = list(Beta=Beta, Gamma=Gamma, V=V, rho=rhopw[rho,1], sigma=sigma, Eta=Eta, Lambda=Lambda, Alpha=Alpha, Psi=Psi, Delta=Delta)
}

Hmsc$set("private", "combineParameters", combineParameters, overwrite=TRUE)



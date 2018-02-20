combineParameters = function(Beta,Gamma,iV,iSigma,Eta,Lambda,Psi,Delta){
   V = solve(iV)
   sigma = 1/iSigma
   par = list(Beta=Beta, Gamma=Gamma, V=V, sigma=sigma, Eta=Eta, Lambda=Lambda, Psi=Psi, Delta=Delta)
}

Hmsc$set("private", "combineParameters", combineParameters, overwrite=TRUE)



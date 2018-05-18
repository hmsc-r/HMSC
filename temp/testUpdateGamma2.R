V0 = UGamma0
V

inv = function(A){
   return(chol2inv(chol(A)))
}
cholL = function(A){
   return(t(chol(A)))
}

Vz1 = kronecker(diag(ns),X)%*%(kronecker(tcrossprod(Tr),V0)+kronecker(diag(ns),V))%*%kronecker(diag(ns),t(X)) + kronecker(diag(ns),diag(ny))
Vz2 = kronecker(diag(ns),diag(ny)) + tcrossprod(kronecker(diag(ns),X)%*%t(chol(kronecker(tcrossprod(Tr),V0)+kronecker(diag(ns),V))))
Vz = Vz2
iVz1 = chol2inv(chol(Vz))

iV = chol2inv(chol(V))
iV0 = chol2inv(chol(V0))
XX = crossprod(X)
TT = crossprod(Tr)

iVz2 = kronecker(diag(d^-1), diag(ny)) - tcrossprod( kronecker(diag(d^-1),X) %*%
      cholL(inv( inv( kronecker(tcrossprod(Tr),V0)+kronecker(diag(ns),V) ) + kronecker(diag(d^-1),XX)) ) )



iVz3 = kronecker(diag(d^-1), diag(ny)) - tcrossprod( kronecker(diag(d^-1),X) %*%
      cholL(inv( kronecker(diag(ns),iV)+kronecker(diag(d^-1),XX) -
            tcrossprod(kronecker(Tr,iV) %*% cholL(inv( kronecker(diag(nt),iV0)+kronecker(TT,iV) ) ) ))))


iP = chol2inv(chol(iV + XX))
LiP = cholL(iP)
R = inv( kronecker(diag(nt),iV0) + kronecker(TT,iV-tcrossprod(iV%*%LiP)) )
LR = cholL(R)
iVz4 = kronecker(diag(ns),diag(ny)-tcrossprod(X%*%LiP)) - tcrossprod(kronecker(Tr,X%*%iP%*%iV) %*% LR)


hist(iVz4%*%Vz - diag(ns*ny))

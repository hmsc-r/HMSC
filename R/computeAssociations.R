#' @title Hmsc$computeAssociations
#'
#' @description Computes the species association matrices
#' @param start
#'
#'@examples
#'

computeAssociations = function(start=1){
   m = self
   OmegaCor = vector("list", m$nr)
   postList=poolMcmcChains(m$postList, start = start)
   getOmegaCor = function(a,r=r)
      return(cov2cor(crossprod(a$Lambda[[r]])))
   for (r in 1:m$nr){
      OmegaCor1 = lapply(postList, getOmegaCor, r=r)
      mOmegaCor1 = apply(abind(OmegaCor1,along=3),c(1,2),mean)
      OmegaCor2 = lapply(OmegaCor1, function(a) return(a>0))
      support1 = apply(abind(OmegaCor2,along=3),c(1,2),mean)
      colnames(mOmegaCor1) = m$spNames
      rownames(mOmegaCor1) = m$spNames
      colnames(support1) = m$spNames
      rownames(support1) = m$spNames
      tmp = list()
      tmp$mean = mOmegaCor1
      tmp$support = support1
      OmegaCor[[r]] = tmp
   }
   return(OmegaCor)
}

Hmsc$set("public", "computeAssociations", computeAssociations, overwrite=TRUE)

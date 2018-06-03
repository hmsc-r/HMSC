#' @title Hmsc$computeVariancePartitioning
#'
#' @description Computes variance partitionint w.r.t. given grouping of fixed efffects and levels of random effects
#' @param group
#' @param groupnames
#' @param start
#'
#' @examples
#'

computeVariancePartitioning = function(group, groupnames, start=1){
   m = self

   ny = m$ny
   nc = m$nc
   nt = m$nt
   ns = m$ns
   np = m$np
   nr = m$nr

   ngroups = max(group);
   fixed = matrix(0,nrow=ns,ncol=1);
   fixedsplit = matrix(0,nrow=ns,ncol=ngroups);
   random = matrix(0,nrow=ns,ncol=nr);
   traitR2 = 0
   cM = cov(m$X)

   postList=poolMcmcChains(m$postList, start=start)

   for (i in 1:m$samples){
      fixed1 = matrix(0,nrow=ns,ncol=1);
      fixedsplit1 = matrix(0,nrow=ns,ncol=ngroups);
      random1 = matrix(0,nrow=ns,ncol=nr);

      Beta = postList[[i]]$Beta
      Gamma = t(postList[[i]]$Gamma)
      Lambdas = postList[[i]]$Lambda

      mu = t(m$Tr%*%Gamma)
      f = m$X%*%Beta
      a = m$X%*%mu
      res1 = 0
      res2 = 0
      for (j in 1:ny){
         res1 = res1 + cov(a[j,],f[j,])^2
         res2 = res2 + (cov(a[j,],a[j,])*cov(f[j,],f[j,]))
      }
      traitR2 = traitR2 + res1/res2

      for (j in 1:ns){
         ftotal = t(Beta[,j])%*%cM%*%Beta[,j]
         fixed1[j] = fixed1[j] + ftotal
         for (k in 1:ngroups){
            sel = (group==k)
            fpart = t(Beta[sel,j])%*%cM[sel,sel]%*%Beta[sel,j]
            fixedsplit1[j,k] = fixedsplit1[j,k] + fpart
         }
      }

      for (level in 1:nr){
         Lambda = Lambdas[[level]]
         nf = dim(Lambda)[[1]]
         for (factor in 1:nf){
            random1[,level] = random1[,level] + t(Lambda[factor,])*Lambda[factor,]
         }
      }
      if (nr>0){
         tot = fixed1+rowSums(random1)
         fixed = fixed + fixed1/tot
         for (level in 1:nr){
            random[,level] = random[,level] + random1[,level]/tot
         }
      } else
      {
         fixed = fixed + matrix(1,nrow=ns,ncol=1)
      }
      for (k in 1:ngroups){
         fixedsplit[,k] =  fixedsplit[,k] + fixedsplit1[,k]/rowSums(fixedsplit1)
      }
   }
   fixed = fixed/m$samples
   random = random/m$samples
   fixedsplit = fixedsplit/m$samples
   traitR2 = traitR2/m$samples


   vals = matrix(0,nrow=ngroups+nr,ncol=ns)
   for (i in 1:ngroups){
     vals[i,] = fixed*fixedsplit[,i]

   }
   for (i in 1:nr){
     vals[ngroups+i,] = random[,i]
   }

   VP = list()
   VP$vals = vals
   VP$traitR2 = traitR2
   VP$group = group
   VP$groupnames = groupnames

   return(VP)
}

Hmsc$set("public", "computeVariancePartitioning", computeVariancePartitioning, overwrite=TRUE)

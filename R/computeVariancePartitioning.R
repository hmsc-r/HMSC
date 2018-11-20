#' @title Hmsc$computeVariancePartitioning
#'
#' @description Computes variance partitions with respect to given grouping of fixed efffects and levels of random effects
#' @param group vector of numeric values corresponding to group identifiers in groupnames
#' @param groupnames vector of names for each random and fixed effect
#' @param start index of first MCMC sample included
#'
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

computeVariancePartitioning = function(hM, group, groupnames, start=1){
   ny = hM$ny
   nc = hM$nc
   nt = hM$nt
   ns = hM$ns
   np = hM$np
   nr = hM$nr

   ngroups = max(group);
   fixed = matrix(0,nrow=ns,ncol=1);
   fixedsplit = matrix(0,nrow=ns,ncol=ngroups);
   random = matrix(0,nrow=ns,ncol=nr);
   traitR2 = 0
   cM = cov(hM$X)

   postList=poolMcmcChains(hM$postList, start=start)

   geta = function(a)
      return(hM$X%*%(t(hM$Tr%*%t(a$Gamma))))
   la=lapply(postList, geta)
   getf = function(a)
      return(hM$X%*%(a$Beta))
   lf=lapply(postList, getf)

   for (i in 1:hM$samples){
      fixed1 = matrix(0,nrow=ns,ncol=1);
      fixedsplit1 = matrix(0,nrow=ns,ncol=ngroups);
      random1 = matrix(0,nrow=ns,ncol=nr);
      Beta = postList[[i]]$Beta
      Lambdas = postList[[i]]$Lambda

      f = lf[[i]]
      a = la[[i]]

      a=a-matrix(rep(rowMeans(a),hM$ns),ncol=hM$ns)
      f=f-matrix(rep(rowMeans(f),hM$ns),ncol=hM$ns)
      res1 = sum((rowSums((a*f))/(hM$ns-1))^2)
      res2 = sum((rowSums((a*a))/(hM$ns-1))*(rowSums((f*f))/(hM$ns-1)))
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

      for (level in seq_len(nr)){
         Lambda = Lambdas[[level]]
         nf = dim(Lambda)[[1]]
         for (factor in 1:nf){
            random1[,level] = random1[,level] + t(Lambda[factor,])*Lambda[factor,]
         }
      }
      if (nr>0){
         tot = fixed1+rowSums(random1)
         fixed = fixed + fixed1/tot
         for (level in seq_len(nr)){
            random[,level] = random[,level] + random1[,level]/tot
         }
      } else{
         fixed = fixed + matrix(1,nrow=ns,ncol=1)
      }
      for (k in 1:ngroups){
         fixedsplit[,k] =  fixedsplit[,k] + fixedsplit1[,k]/rowSums(fixedsplit1)
      }
   }
   fixed = fixed/hM$samples
   random = random/hM$samples
   fixedsplit = fixedsplit/hM$samples
   traitR2 = traitR2/hM$samples


   vals = matrix(0,nrow=ngroups+nr,ncol=ns)
   for (i in 1:ngroups){
      vals[i,] = fixed*fixedsplit[,i]

   }
   for (i in seq_len(nr)){
      vals[ngroups+i,] = random[,i]
   }

   VP = list()
   VP$vals = vals
   VP$traitR2 = traitR2
   VP$group = group
   VP$groupnames = groupnames

   return(VP)
}

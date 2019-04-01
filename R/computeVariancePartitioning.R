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
#' returns an object VP with components VP$vals, VP$R2T, VP$group and VP$groupnames.
#'
#' @details
#'
#' The vector \code{group} has one value for each column of the matrix \code{hM$X}, describing the index of the
#' group to which this column is to be included. The names of the group are given by \code{groupnames}. The output object
#' \code{VP$vals} gives the variance proportion for each group and species. The output object \code{VP$R2T} has shows the
#' variance among species explained by traits, measured for species niches (\code{VP$R2T$Beta}) and species occurrences
#' (\code{VP$R2T$Y})
#'
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' VP = computeVariancePartitioning(m, group=c(1,1,1,2,2), groupnames = c("habitat quality","climate"))
#' VP$R2T
#' plotVariancePartitioning(m,VP)
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
   R2T.Y = 0
   R2T.Beta = rep(0,nc)

   switch(class(hM$X),
          matrix = {
            cMA = cov(hM$X)
          },
          list = {
            cMA = lapply(hM$X, cov)
          }
   )

   postList=poolMcmcChains(hM$postList, start=start)

   geta = function(a){
     switch(class(hM$X),
            matrix = {res = hM$X%*%(t(hM$Tr%*%t(a$Gamma)))},
            list = {
              res = matrix(NA,hM$ny,hM$ns)
              for(j in 1:hM$ns) res[,j] =  hM$X[[j]]%*%(t(hM$Tr[j,]%*%t(postList[[1]]$Gamma)))
            }
     )
     return(res)
   }
   la=lapply(postList, geta)

   getf = function(a){
     switch(class(hM$X),
            matrix = {res = hM$X%*%(a$Beta)},
            list = {
             res = matrix(NA,hM$ny,hM$ns)
             for(j in 1:hM$ns) res[,j] = hM$X[[j]]%*%a$Beta[,j]
            }
     )
     return(res)
   }
   lf=lapply(postList, getf)

   lmu=lapply(postList, gemu)
   gemu = function(a){
      res = t(hM$Tr%*%t(a$Gamma))
      return(res)
   }

      gebeta = function(a){
      res = a$Beta
      return(res)
   }
   lbeta=lapply(postList, gebeta)


   for (i in 1:hM$samples){
      for (k in 1:nc){
         R2T.Beta[k] = R2T.Beta[k] + cor(lbeta[[i]][k,],lmu[[i]][k,])^2
      }

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
      R2T.Y = R2T.Y + res1/res2
      for (j in 1:ns){
         switch(class(hM$X), matrix = {cM = cMA}, list = {cM = cMA[[j]]})
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
   R2T.Y = R2T.Y/hM$samples
   R2T.Beta = R2T.Beta/hM$samples


   vals = matrix(0,nrow=ngroups+nr,ncol=ns)
   for (i in 1:ngroups){
      vals[i,] = fixed*fixedsplit[,i]

   }
   for (i in seq_len(nr)){
      vals[ngroups+i,] = random[,i]
   }

   names(R2T.Beta)=hM$covNames
   colnames(vals)=hM$spNames
   leg = groupnames
   for (r in 1:hM$nr) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
   }
   rownames(vals)=leg

   VP = list()
   VP$vals = vals
   VP$R2T = list(Beta=R2T.Beta,Y=R2T.Y)
   VP$group = group
   VP$groupnames = groupnames

   return(VP)
}

#' @title computeVariancePartitioning
#'
#' @description Computes variance components with respect to given grouping of fixed effects and levels
#' of random effects
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param group vector of numeric values corresponding to group
#'     identifiers in groupnames. If the model was defined with
#'     \code{XData} and \code{XFormula}, the default is to use model
#'     terms.
#' @param groupnames vector of names for each group of fixed
#'     effect. Should match \code{group}. If the model was defined
#'     with \code{XData} and \code{XFormula}, the default is to use
#'     the labels of model terms.
#' @param start index of first MCMC sample included
#' @param na.ignore logical. If TRUE, covariates are ignored for sites
#'     where the focal species is NA when computing
#'     variance-covariance matrices for each species
#'
#'
#' @return
#' returns an object VP with components VP$vals, VP$R2T, VP$group and VP$groupnames.
#'
#' @details
#' The vector \code{group} has one value for each column of the matrix \code{hM$X}, describing the index of the
#' group in which this column is to be included. The names of the group are given by \code{groupnames}. The output object
#' \code{VP$vals} gives the variance proportion for each group and species. The output object \code{VP$R2T} gives the
#' variance among species explained by traits, measured for species' responses to covariates (\code{VP$R2T$Beta}) and species occurrences
#' (\code{VP$R2T$Y})
#'
#' @seealso
#' Use \code{\link{plotVariancePartitioning}} to display the result object.
#'
#' @examples
#' # Partition the explained variance for a previously fitted model
#' # without grouping environmental covariates
#' VP = computeVariancePartitioning(TD$m)
#'
#' # Partition the explained variance for a previously fitted model
#' # while grouping the two environmental variables together
#' VP = computeVariancePartitioning(TD$m, group=c(1,1), groupnames = c("Habitat"))
#'
#' @importFrom stats cov cor terms
#' @export

computeVariancePartitioning =
    function(hM, group=NULL, groupnames=NULL, start=1, na.ignore=FALSE)
{
   ny = hM$ny
   nc = hM$nc
   nt = hM$nt
   ns = hM$ns
   np = hM$np
   nr = hM$nr

   if(is.null(group)){
      ## default: use terms
      if(nc > 1){
         if (is.null(hM$XFormula))
             stop("no XFormula: you must give 'group' and 'groupnames'")
         group = attr(hM$X, "assign")
         if (group[1] == 0) # assign (Intercept) to group
            group[1] <- 1
         groupnames = attr(terms(hM$XFormula, data = hM$XData), "term.labels")
      } else {
         group = c(1)
         groupnames = hM$covNames[1]
      }

   }
   ngroups = max(group);
   fixed = matrix(0,nrow=ns,ncol=1);
   fixedsplit = matrix(0,nrow=ns,ncol=ngroups);
   random = matrix(0,nrow=ns,ncol=nr);

   R2T.Y = 0
   R2T.Beta = rep(0,nc)

   #If na.ignore=T, convert XData to a list
   if(na.ignore){
      xl=list()
      for(s in 1:ns){
         xl[[s]]=hM$X
      }
      hM$X=xl
   }

   switch(class(hM$X)[1L], matrix = {
      cMA = cov(hM$X)
   }, list = {

      if(na.ignore){
         cMA = list()
         for(s in 1:ns){cMA[[s]]=cov(hM$X[[s]][which(hM$Y[,s]>-Inf),])}
      }
      else{cMA = lapply(hM$X, cov)}
   })

   postList=poolMcmcChains(hM$postList, start=start)

   geta = function(a){
     switch(class(hM$X)[1L],
            matrix = {res = hM$X%*%(t(hM$Tr%*%t(a$Gamma)))},
            list = {
              res = matrix(NA,hM$ny,hM$ns)
              for(j in 1:hM$ns) res[,j] =  hM$X[[j]]%*%(t(hM$Tr[j,]%*%t(a$Gamma)))
            }
     )
     return(res)
   }
   la=lapply(postList, geta)

   getf = function(a){
     switch(class(hM$X)[1L],
            matrix = {res = hM$X%*%(a$Beta)},
            list = {
             res = matrix(NA,hM$ny,hM$ns)
             for(j in 1:hM$ns) res[,j] = hM$X[[j]]%*%a$Beta[,j]
            }
     )
     return(res)
   }
   lf=lapply(postList, getf)

   gemu = function(a){
      res = t(hM$Tr%*%t(a$Gamma))
      return(res)
   }
   lmu=lapply(postList, gemu)


      gebeta = function(a){
      res = a$Beta
      return(res)
   }
   lbeta=lapply(postList, gebeta)

   poolN <- length(postList) # pooled chains
   for (i in seq_len(poolN)){
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
         switch(class(hM$X)[1L], matrix = {cM = cMA}, list = {cM = cMA[[j]]})
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
   fixed = fixed/poolN
   random = random/poolN
   fixedsplit = fixedsplit/poolN
   R2T.Y = R2T.Y/poolN
   R2T.Beta = R2T.Beta/poolN

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
   for (r in seq_len(nr)) {
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

#' @title computeInitialParameters
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'



computeDataParameters = function(){
   parList = list()

   if(!is.null(self$C)){
      iQg = array(NA, c(self$ns,self$ns,nrow(self$rhopw)))
      RiQg = array(NA, c(self$ns,self$ns,nrow(self$rhopw)))
      detQg = rep(NA, nrow(self$rhopw))
      if(any(self$rhopw[,1] < 0))
         iC = chol2inv(chol(self$C))
      for(rg in 1:nrow(self$rhopw)){
         rho = self$rhopw[rg,1]
         if(rho >= 0){
            rhoC = rho*self$C;
         } else{
            rhoC = (-rho)*iC;
         }
         Q = rhoC + (1-abs(rho))*diag(self$ns)
         RQ = chol(Q);
         iQg[,,rg] = chol2inv(RQ)
         RiQg[,,rg] = t(backsolve(RQ, diag(self$ns)))
         detQg[rg] = 2*sum(log(diag(RQ)))
      }
   } else{
      iQg = NULL
      detQg = NULL
      RiQg = NULL
   }

   rLPar = vector("list", self$nr)
   for(r in 1:self$nr){
      if(self$rL[[r]]$sDim > 0){
         alphapw = self$rL[[r]]$alphapw
         np = self$np[r]
         s = self$rL[[r]]$s[levels(self$dfPi[,r]),]
         alphaN = nrow(alphapw)
         distance = as.matrix(dist(s))

         iWg = array(NA, c(np,np,alphaN))
         RiWg = array(NA, c(np,np,alphaN))
         detWg = rep(NA, alphaN)
         for(ag in 1:alphaN){
            alpha = alphapw[ag,1]
            if(alpha==0){
               W = diag(np)
            } else{
               W = exp(-distance/alpha)
            }
            RW = chol(W)
            iW = chol2inv(RW)

            iWg[,,ag] = iW
            RiWg[,,ag] = chol(iW)
            detWg[ag] = 2*sum(log(diag(RW)))
         }
         rLPar[[r]] = list(iWg=iWg, RiWg=RiWg, detWg=detWg)
      }
   }
   parList$iQg = iQg
   parList$RiQg = RiQg
   parList$detQg = detQg
   parList$rLPar = rLPar

   return(parList)
}

Hmsc$set("private", "computeDataParameters", computeDataParameters, overwrite=TRUE)


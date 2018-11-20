#' @title computeInitialParameters
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
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


computeDataParameters = function(hM){
   parList = list()

   if(!is.null(hM$C)){
      iQg = array(NA, c(hM$ns,hM$ns,nrow(hM$rhopw)))
      RiQg = array(NA, c(hM$ns,hM$ns,nrow(hM$rhopw)))
      detQg = rep(NA, nrow(hM$rhopw))
      if(any(hM$rhopw[,1] < 0))
         iC = chol2inv(chol(hM$C))
      for(rg in 1:nrow(hM$rhopw)){
         rho = hM$rhopw[rg,1]
         if(rho >= 0){
            rhoC = rho*hM$C;
         } else{
            rhoC = (-rho)*iC;
         }
         Q = rhoC + (1-abs(rho))*diag(hM$ns)
         RQ = chol(Q);
         iQg[,,rg] = chol2inv(RQ)
         RiQg[,,rg] = t(backsolve(RQ, diag(hM$ns)))
         detQg[rg] = 2*sum(log(diag(RQ)))
      }
   } else{
      iQg = NULL
      detQg = NULL
      RiQg = NULL
   }

   rLPar = vector("list", hM$nr)
   for(r in seq_len(hM$nr)){
      if(hM$rL[[r]]$sDim > 0){
         alphapw = hM$rL[[r]]$alphapw
         np = hM$np[r]
         s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
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


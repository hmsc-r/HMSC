#' @title Hmsc$computeMarginalEffects
#'
#' @description Computes ...
#' @param ngrid
#' @param prob
#' @param clist
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

computeMarginalEffects = function(ngrid=20, prob=c(0.025,0.5,0.975), clist=1:self$nc){
   m = self
   nc = length(clist)

   pred = vector("list", nc)
   for(ii in 1:nc){
      i = clist[[ii]]
      X = matrix(NA, ncol=m$nc, nrow=ngrid)
      x = m$X[,i]
      mi = min(x)
      ma = max(x)
      xx = seq(mi, ma, length.out=ngrid)
      X[,i] = xx

      for(j in 1:m$nc){
         if(j!=i){
            xv = m$X[,2]
            yv = m$X[,j]
            mylm = lm(yv~xv)
            yy = predict(mylm, newdata=data.frame(xv=xx))
            X[,j] = yy
         }
      }

      dfPi = matrix(NA, ngrid, m$nr)
      for(r in seq_len(m$nr)){
         dfPi[,r] = sprintf('new_unit', 1:(ngrid))
      }
      dfPi = as.data.frame(dfPi)
      postList = poolMcmcChains(m$postList)

      rL = vector("list", m$nr)
      for (r in seq_len(m$nr)){
         tmp = m$rL[[r]]
         rL1 = tmp$clone()
         units = rL1$pi
         units1 = c(units,"new_unit")
         xydata = rL1$s
         if(!is.null(xydata)){
            nxy = dim(xydata)[1]
            xydata1 = matrix(NA, nrow=nxy+1, ncol=rL1$sDim)
            xydata1[1:nxy,] = xydata
            xydata1[nxy+1,] = colMeans(xydata)
            rownames(xydata1) = units1
            colnames(xydata1) = colnames(xydata)
            rL1$s = xydata1
         }
         rL1$pi = units1
         rL1$N = rL$N+1
         rL[[r]] = rL1
      }

      predY1 = m$predict(post=postList, X=X, dfPiNew=dfPi, rL=rL, expected=TRUE)
      for (j in 1:length(predY1)){
         colnames(predY1[[j]]) = m$spNames
      }
      predS1 = lapply(predY1, rowSums)
      predT1 = lapply(predY1, function(a) (a%*%m$Tr)/matrix(rep(rowSums(a),m$nt),ncol=m$nt))

      qpredY1 = apply(abind(predY1,along=3),c(1,2),quantile,prob=c(0,prob))
      qpredS1 = apply(abind(predS1,along=2),c(1),quantile,prob=c(0,prob))
      qpredT1 = apply(abind(predT1,along=3),c(1,2),quantile,prob=c(0,prob))

      qpredS1[1,] = xx
      rownames(qpredS1)[[1]] = m$covNames[[i]]

      for (j in 1:m$nt){
         qpredT1[1,,j] = xx
      }
      rownames(qpredT1)[1] = m$covNames[[i]]

      for (j in 1:m$ns){
         qpredY1[1,,j] = xx
      }
      rownames(qpredY1)[1] = m$covNames[[i]]

      pred1=list()
      pred1$Y = qpredY1
      pred1$S = qpredS1
      pred1$Tr = qpredT1

      pred[[ii]] = pred1
   }
   return(pred)
}

Hmsc$set("public", "computeMarginalEffects", computeMarginalEffects, overwrite=TRUE)

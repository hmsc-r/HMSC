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
#' @export

computeMarginalEffects = function(hM, ngrid=20, prob=c(0.025,0.5,0.975), clist=1:self$nc){
   nc = length(clist)

   pred = vector("list", nc)
   for(ii in 1:nc){
      i = clist[[ii]]
      X = matrix(NA, ncol=hM$nc, nrow=ngrid)
      x = hM$X[,i]
      mi = min(x)
      ma = max(x)
      xx = seq(mi, ma, length.out=ngrid)
      X[,i] = xx

      for(j in 1:hM$nc){
         if(j!=i){
            xv = hM$X[,2]
            yv = hM$X[,j]
            mylm = lm(yv~xv)
            yy = predict(mylm, newdata=data.frame(xv=xx))
            X[,j] = yy
         }
      }

      dfPi = matrix(NA, ngrid, hM$nr, dimnames=c(NULL, hM$rLNames))
      for(r in seq_len(hM$nr)){
         dfPi[,r] = sprintf('new_unit', 1:(ngrid))
      }
      dfPi = as.data.frame(dfPi)
      postList = poolMcmcChains(hM$postList)

      rL = vector("list", hM$nr)
      for (r in seq_len(hM$nr)){
         rL1 = hM$rL[[r]]
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

      predY1 = predict(hM, post=postList, X=X, dfPiNew=dfPi, rL=rL, expected=TRUE)
      for (j in 1:length(predY1)){
         colnames(predY1[[j]]) = hM$spNames
      }
      predS1 = lapply(predY1, rowSums)
      predT1 = lapply(predY1, function(a) (a%*%hM$Tr)/matrix(rep(rowSums(a),hM$nt),ncol=hM$nt))

      qpredY1 = apply(abind(predY1,along=3),c(1,2),quantile,prob=c(0,prob))
      qpredS1 = apply(abind(predS1,along=2),c(1),quantile,prob=c(0,prob))
      qpredT1 = apply(abind(predT1,along=3),c(1,2),quantile,prob=c(0,prob))

      qpredS1[1,] = xx
      rownames(qpredS1)[[1]] = hM$covNames[[i]]

      for (j in 1:hM$nt){
         qpredT1[1,,j] = xx
      }
      rownames(qpredT1)[1] = hM$covNames[[i]]

      for (j in 1:hM$ns){
         qpredY1[1,,j] = xx
      }
      rownames(qpredY1)[1] = hM$covNames[[i]]

      pred1=list()
      pred1$Y = qpredY1
      pred1$S = qpredS1
      pred1$Tr = qpredT1

      pred[[ii]] = pred1
   }
   return(pred)
}

#' @title combineHmscChains
#'
#' @description Combines different chains of Hmsc model fitted with different arguments in sampleMcmc(..., indChains, ...) call
#' @export

combineHmscChains = function(hMList, alignPost=TRUE){
   lenVecPostList = sapply(hMList, function(m) length(m$postList))
   if(length(unique(lenVecPostList)) != 1){
      stop("lengths of postList fields are not identical")
   }
   lenPostList = unique(lenVecPostList)
   indChainCount = rep(0, lenPostList)
   hM = hMList[[1]]
   for(i in 1:length(hMList)){
      hasChainVec = as.logical(1-sapply(hMList[[i]]$postList, is.null))
      indChainCount = indChainCount + hasChainVec
      hM$postList[which(hasChainVec)] = hMList[[i]]$postList[which(hasChainVec)]
   }
   print(indChainCount)
   if(any(indChainCount != 1)){
      tmpStr = sprintf("[%s]",paste(sprintf("%d", indChainCount),collapse=","))
      stop(sprintf("combined model chains do not form a partition of 1:nChains \n individual chains were fitted %s times", tmpStr))
   }
   if(alignPost==TRUE){
      for(i in 1:5){
         hM = alignPosterior(hM)
      }
   }
   return(hM)
}

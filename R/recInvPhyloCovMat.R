# @title recInvPhyloCovMat
#
#' @importFrom Matrix Matrix Diagonal bdiag
#'
#'
recInvPhyloCovMat = function(s,rho){
   # print(s)
   openPar = 0
   startInd = 1
   compNum = 0
   compList = list()
   for(i in 1:nchar(s)){
      if(substr(s,i,i)=="("){
         openPar = openPar + 1
      } else if(substr(s,i,i)==")"){
         openPar = openPar - 1
      } else if(substr(s,i,i)=="," && openPar==0){
         compNum = compNum + 1
         compList[[compNum]] = substr(s,startInd,i-1)
         startInd = i+1
      }
   }
   compNum = compNum + 1
   compList[[compNum]] = substr(s,startInd,nchar(s))
   # print(compList)
   matList = selfMatList = leafFlagVecList = vector("list", length(compList))
   for(i in 1:length(compList)){
      sTmp = compList[[i]]
      if(length(grep("\\(",sTmp)) == 0){
         sParts = unlist(strsplit(gsub("^(.*)[:](.*)$", "\\1€\\2", sTmp),"€"))
         # print(sParts)
         nodeName = sParts[1]
         bL = as.numeric(sParts[2])
         mat11 = Matrix((bL*rho + (1-rho))^-1)
         matList[[i]] = mat11
         rownames(matList[[i]]) = colnames(matList[[i]]) = nodeName
         leafFlagVecList[[i]] = 1
      } else{
         sParts = unlist(strsplit(gsub("^\\((.*)\\)(.*)[:](.*)$", "\\1€\\2€\\3", sTmp),"€"))
         # print(sParts)
         nodeName = sParts[2]
         bL = as.numeric(sParts[3])
         mat11 = (bL*rho)^-1
         recRes = recInvPhyloCovMat(sParts[1],rho)
         matRec = recRes[[1]]
         childInd = recRes[[2]]
         childMatList = recRes[[3]]
         childFlagVec = recRes[[4]]
         tmp1 = Matrix(0,nrow(matRec),1)
         for(j in 1:length(childMatList)){
            tmp1[childInd[j],1] = childMatList[[j]]
         }
         tmp2 = Reduce("+",childMatList)
         mat1 = rbind(matRec, -t(tmp1))
         mat2 = rbind(-tmp1, mat11+tmp2)
         matList[[i]] = cbind(mat1,mat2)
         rownames(matList[[i]]) = colnames(matList[[i]]) = c(rownames(matRec), nodeName)
         leafFlagVecList[[i]] = c(childFlagVec,0)
      }
      selfMatList[[i]] = mat11
   }
   mat = bdiag(matList)
   rownames(mat) = colnames(mat) = unlist(lapply(matList, rownames))
   childInd = cumsum(unlist(lapply(matList,nrow)))
   leafFlagVec = Reduce(c, leafFlagVecList)
   # print(mat)
   # print(rownames(mat))
   return(list(mat=mat,childInd=childInd,childMatList=selfMatList,leafFlagVec=leafFlagVec))
}

#' @title versionUpgrade
#'
#' @description Attempts to upgrade the model object version to be in line with the latest Hmsc version
#' @export

versionUpgrade = function(hM){
   if(compareVersion(as.character(hM$HmscVersion), "3.0")>=0 && compareVersion(as.character(hM$HmscVersion), "3.1")==-1){
      hM = versionUpgrade_3.1(hM)
   }
   # if(compareVersion(as.character(hM$HmscVersion), "3.1")>=0 && compareVersion(as.character(hM$HmscVersion), "3.2")==-1){
   #    hM = versionUpgrade_3.2(hM)
   # }
   hM$HmscVersion = packageVersion("Hmsc")
   return(hM)
}

versionUpgrade_3.1 = function(hM){
   hM$rhoAlphaDP = 0
   # if(hM$postList){
   #    for(chain in 1:length(postList)){
   #       for(sam in 1:length(postList[[chain]])){
   #          postList[[chain]][[sam]]$rho = rep(postList[[chain]][[sam]]$rho, hM$nc)
   #       }
   #    }
   # }
   #TODO that to do with parameters computed in computeDataParList?
   hM$HmscVersion = package_version("3.1.1")
   return(hM)
}


versionUpgrade_3.2 = function(hM){
   "............"
   return(hM)
}

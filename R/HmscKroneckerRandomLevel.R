#' @title Create an \code{Hmsc} kronecker random level structure
#'
#' @description Specifies the a kronecker random factor level
#' @param rLList
#'
#' @return a \code{HmscKroneckerRandomLevel}-class object that can be used for \code{Hmsc}-class object construction
#'
#' @details Nothing yet.
#'
#'
#' @seealso [setPriors.Hmsc()]
#'
#' @examples
#' # Setting a random level with 50 units
#' rL = HmscRandomLevel(units=TD$studyDesign$sample)
#'
#' # Setting a spatial random level
#' rL = HmscRandomLevel(sData=TD$xycoords)
#'
#' # Setting a covariate-dependent random level.
#' rL = HmscRandomLevel(xData=data.frame(x1=rep(1,length(TD$X$x1)),x2=TD$X$x2))
#'
#' @importFrom plyr mdply
#'
#' @export

HmscKroneckerRandomLevel = function(rLList, sepStr="=", etaMethod="R", alphaMethod="R", alphaMarginalMethod=FALSE,
                                    logDetKstBatchSize=0, cgIterN=NULL) {
   kRL = structure(list(rLList=NULL, kN=NULL, pi=NULL, sDim=NULL, xDim=NULL, N=NULL, sepStr=sepStr, #
                        nfMax=NULL, nfMin=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphaPrior=NULL,
                        etaMethod=etaMethod, alphaMethod=alphaMethod, alphaMarginalMethod=alphaMarginalMethod,
                        logDetKstBatchSize=logDetKstBatchSize, cgIterN=cgIterN),
                   class=c("HmscKroneckerRandomLevel","HmscRandomLevel"))
   kRL$rLList = rLList #TODO introduce some checks
   if(length(rLList) != 2){
      stop("HmscKroneckerRandomLevel: currently the number of elements in Kronecker random level must be exactly 2.")
   }
   kRL$kN = length(kRL$rLList)
   dfTemp = expand.grid(lapply(rev(rLList), function(rL) rL$pi))[rev(1:length(rLList))]
   kRL$pi = as.factor(mdply(dfTemp, paste, sep=sepStr)[,length(kRL$rLList)+1])
   kRL$N = Reduce(prod, lapply(rLList, function(rL) rL$N))

   if(Reduce(all, lapply(rLList, function(rL) rL$sDim==0))){
      stop("HmscKroneckerRandomLevel: all kronecker random level elements are non-spatial - use a simple random level instead")
      kRL$sDim = 0
   } else{
      kRL$sDim = Inf
   }
   kRL$xDim = 0

   alphaPrior = lapply(kRL$rLList, function(a) a$alphapw)
   kRL = setPriors(kRL, nu=kRL$rLList[[1]]$nu, a1=kRL$rLList[[1]]$a1, b1=kRL$rLList[[1]]$b1, a2=kRL$rLList[[1]]$a2, b2=kRL$rLList[[1]]$b2,
                   nfMax=kRL$rLList[[1]]$nfMax, nfMin=kRL$rLList[[1]]$nfMin,
                   alphaPrior=alphaPrior)
   kRL$call <- match.call()
   return(kRL)
}

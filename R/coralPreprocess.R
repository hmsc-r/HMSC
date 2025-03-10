#' @title coralPreprocess
#'
#' @description Splits typical Hmsc data into common and rare partitions for CORAL analysis
#'
#' @param Y community matrix of both common and rare species
#' @param spNames.common vector of species that are considered common in CORAL analysis
#' @param TrData trait matrix for both common and rare species
#' @param phyloTree phylogeny tree covering both common and rare species
#' @param Taxa dataframe with
#' @param TaxaFormula formula for calculating phylogeny from taxonomy, as in \code{as.phylo.formula(...)}
#'
#' @return
#' A named list containing \code{Y.xyz} and \code{TrData.xyz} parts for common and rare partitions,
#' phylogeny tree \code{phyloTree.common} for common part and
#' phylogeny similarity matrix \code{C.common.rare} between common and rare species
#'
#' @details
#' This functions implies that all column names of \code{Y} that are not listed in \code{spNames.common} argument
#' are considered rare species in the context of CORAL analysis.
#'
#' Exactly one argument of phyloTree and Taxa must be specified.
#'
#' @importFrom ape as.phylo.formula rtree vcv
#' @importFrom stats as.formula
#'
#' @export


coralPreprocess = function(Y, spNames.common, TrData=NULL, phyloTree=NULL, Taxa=NULL, TaxaFormula=NULL){
   if(!is.null(phyloTree) && !is.null(Taxa)){
      stop("only one of phyloTree and Taxa arguments can be specified")
   }
   if(is.null(phyloTree) && is.null(Taxa)){
      stop("exactly one of phyloTree and Taxa arguments shall be specified for CORAL approach")
   }
   if(!is.null(Taxa) && is.null(TaxaFormula)){
      warning("no TaxaFormula specified, using all columns of Taxa")
      TaxaFormula = as.formula(paste0("~", paste(colnames(Taxa), collapse="/")))
   }
   ind.common = match(spNames.common, colnames(Y))
   ns.common = length(ind.common)
   ind.rare = setdiff(1:ncol(Y), ind.common)
   spNames.rare = colnames(Y)[ind.rare]
   ns.rare = length(ind.rare)
   Y.common = Y[,ind.common]
   Y.rare = Y[,match(spNames.rare, colnames(Y))]
   if(!is.null(TrData)){
      TrData.common = TrData[match(spNames.common, rownames(TrData)),]
      TrData.rare = TrData[match(spNames.rare, rownames(TrData)), ]
   } else{
      TrData.common = NULL
      TrData.rare = NULL
   }
   if(!is.null(phyloTree)){
      phyloTree.common = keepTipRoot(phyloTree, spNames.common)
      phyloTreePart = keepTipRoot(phyloTree, c(spNames.common, spNames.rare))
      C = vcv(phyloTreePart, model="Brownian", corr=TRUE)
      C.common.rare = C[spNames.common, spNames.rare]
   } else{
      Taxa.common = Taxa[match(spNames.common, rownames(TrData)),]
      phyloTree.common = as.phylo.formula(TaxaFormula, data=Taxa.common)
      phyloTree.common$edge.length = rep(1, length(phyloTree.common$edge.length))
      Taxa.rare = Taxa[match(spNames.rare, rownames(TrData)),]
      C12 = matrix(0, ns.common, ns.rare)
      for(k in 1:ncol(Taxa)){
         C12 = C12 + (matrix(Taxa.common[,k],ns.common,ns.rare) == matrix(Taxa.rare[,k],ns.common,ns.rare,byrow=TRUE))
      }
      C12 = C12 / ncol(Taxa)
      colnames(C12) = spNames.rare
      rownames(C12) = spNames.common
      C.common.rare = C12[, match(spNames.rare, colnames(C12))]
   }
   return(list(Y.common=Y.common, TrData.common=TrData.common, phyloTree.common=phyloTree.common,
               Y.rare=Y.rare, TrData.rare=TrData.rare, C.common.rare=C.common.rare))
}

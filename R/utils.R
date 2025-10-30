#' @importFrom ape bind.tree keep.tip drop.tip

keepTipRoot = function(phyloTree, tipNames, tmpTipName="_extra_tip_at_root_"){
   tmpTree1 = bind.tree(phyloTree, rtree(1, tip.label=tmpTipName))
   tmpTree2 = keep.tip(tmpTree1, c(tmpTipName,tipNames), trim.internal=TRUE)
   tmpTree3 = drop.tip(tmpTree2, tmpTipName, trim.internal=FALSE, collapse.singles=FALSE)
   return(tmpTree3)
}

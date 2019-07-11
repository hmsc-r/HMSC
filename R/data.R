#' Simulated data and a fitted Hmsc model for a small species community.
#'
#' This dataset contains simulated occurence data for 4 species in 50 sampling units.
#' The data is based on a hierarchical study design consisting of 50 sampling units in 10 georeferenced plots.
#' Occurences of 4 species were simulated using one continuous environmental variable (x1) and spatial autocorrelation between the plots.
#' Response of species to the environment are related to one species trait which is fully phylogenetically structured.
#' This dataset is used for the examples and package testing.
#' The variables are as follows:
#'
#' @format A list of 12 objects
#' \describe{
#'   \item{ns}{Number of species in the dataset}
#'   \item{units}{Number of sampling units}
#'   \item{plots}{Number of plots}
#'   \item{X}{A 3 by 50 environmental matrix consisting of one continuous and one categorical variable. Also includes intercept column}
#'   \item{phy}{A list containing the simulated phylogenetic tree for 4 species}
#'   \item{C}{A 4 by 4 phylogenetic variance covariance matrix}
#'   \item{Tr}{A 4 by 3 trait matrix with one phylogenetically phylogenetically structured continuous trait, one categroical trait and an intercept}
#'   \item{xycoords}{simulated 2 dimensional coordinates}
#'   \item{studyDesign}{Sampling unit and plot IDs}
#'   \item{Y}{Simulated species occurences}
#'   \item{m}{A fitted Hmsc object with 100 posterior samples}
#' }
"TD"

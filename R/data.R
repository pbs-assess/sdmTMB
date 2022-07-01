#' Example data for Pacific Cod
#'
#' @format A data frame.
"pcod"

#' Example data for Pacific Cod (years 2011 and after)
#'
#' @format A data frame.
"pcod_2011"

#' Example SPDE mesh for the Pacific Cod data (years 2011 and after)
#'
#' @format A list object of class `sdmTMBmesh`.
"pcod_mesh_2011"

#' Example 2x2km prediction grid for Queen Charlotte Sound
#'
#' @format A data frame.
"qcs_grid"

#' BC coastline data from ropensci/rnaturalearthhires
#'
#' @format An sf data frame.
"bc_coast"

#' Tree locations from Barro Colorado Island
#'
#' @format A dataframe, from the spatstat.data package.
#'
#' This is an is an object of class \code{"ppp"} representing the
#' point pattern of tree locations. See \code{\link[spatstat.geom]{ppp.object}}
#' for details of the format.
#'
#' The dataset \code{bei} gives the positions of 3605 trees
#' of the species \emph{Beilschmiedia pendula} (Lauraceae)
#' in a 1000 by 500 metre rectangular sampling region
#' in the tropical rainforest of Barro Colorado Island.
#' All spatial coordinates are given in metres.
#
#' These data are part of a much larger dataset containing the positions of
#' hundreds of thousands of trees belong to thousands of species;
#' see Hubbell and Foster (1983), Condit, Hubbell and Foster (1996)
#' and Condit (1998).
#'
#' The present data were analysed by Moller and Waagepetersen (2007)
#'
#' Condit, R. (1998) Tropical Forest Census Plots.
#' Springer-Verlag, Berlin and R.G. Landes Company, Georgetown, Texas.
#'
#' Condit, R., Hubbell, S.P and Foster, R.B. (1996)
#' Changes in tree species abundance in a neotropical forest: impact of
#' climate change. Journal of Tropical Ecology, 12,
#' 231--256.
#'
#' Hubbell, S.P and Foster, R.B. (1983)
#' Diversity of canopy trees in a neotropical forest and implications for
#' conservation. In: Tropical Rain Forest: Ecology and Management
#' (eds. S.L. Sutton, T.C. Whitmore and A.C. Chadwick),
#' Blackwell Scientific Publications, Oxford, 25--41.
#'
#' Moller, J. and Waagepetersen, R.P. (2007)
#' Modern spatial point process modelling and inference (with discussion).
#' Scandinavian Journal of Statistics 34, 643--711.
#'
#' @source \url{https://github.com/spatstat/spatstat.data}
"bei"

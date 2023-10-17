#' Example fish survey data
#'
#' @description
#' Various fish survey datasets.
#'
#' @format `pcod`: Trawl survey data for Pacific Cod in Queen Charlotte Sound. A
#'   data frame.
#' @rdname surveydata
"pcod"

#' @format `pcod_2011`: A version of `pcod` for years 2011 and after (smaller
#'   for speed). A data frame.
#' @rdname surveydata
"pcod_2011"

#' @format `pcod_mesh_2011`: A mesh pre-built for `pcod_2011` for examples. A
#'   list of class `sdmTMBmesh`.
#' @rdname surveydata
"pcod_mesh_2011"

#' @format `qcs_grid` A 2x2km prediction grid for Queen Charlotte Sound. A data
#'   frame.
#' @rdname surveydata
"qcs_grid"

#' @format `dogfish`: Trawl survey data for Pacific Spiny Dogfish on West Coast
#'   Vancouver Island. A data frame.
#' @rdname surveydata
"dogfish"

#' @format `yelloweye`: Survey data for Yelloweye Rockfish from the Hard Bottom
#'   Longline Survey (South) off West Coast Vancouver Island.
#' @rdname surveydata
"yelloweye"

#' @format `hbll_s_grid`: A survey domain grid to go with `yelloweye`. A data frame.
#' @rdname surveydata
"hbll_s_grid"

#' @format `wcvi_grid`: A survey domain grid to go with `dogfish`. A data frame.
#' @rdname surveydata
"wcvi_grid"

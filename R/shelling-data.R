#' Arabidopsis QTL data on gravitropism
#'
#' Output from 5,000 replications of the 4-populations shelling model from Flache NETLOGO model.
#' We vary three parameters \code{density}, \code{X..similar.wanted}
#' (\code{X..} is the R friendly way of writing percent) and \code{radiusNeighborhood}.
#' The other 77 columns summary statistics are summary statistics (outputs) of the model
#' (plus one column containing random-seed for replication)
#'
#' @docType data
#'
#' @usage data(shelling)
#'
#' @format An object of class \code{"tibble"};
#'
#' @keywords datasets
#'
#' @references Flache, A., & de Matos Fernandes, C. A. (2021). Agent-based computational models. In G. Manzo (Ed.), Research Handbook on Analytical Sociology. Cheltenham: Edward Elgar.
#'
#' Ernesto Carrella (under revision). No free lunch when estimating simulation parameters
#'
#'
#' @source \href{https://www.comses.net/codebases/7c562b23-4964-4d28-862c-1c8b254fd6ad/releases/1.0.0/}{COMSES netlogo code}
#'
#' @examples
#' data(shelling)
#' rf.shelling<-
#' cross_validate_random_forest(shelling,
#'                              ngroup=5,
#'                              parameter_colnames = c("density","radiusNeighborhood","X..similar.wanted"),
#'                              summary_statistics_colnames =
#'                                c("trend.percent.unhappy","trend.percent.unhappy.Red","trend.percent.unhappy.Blue",
#'                                  "trend.percent.similar"),
#'                              fast=TRUE)
#' rf.shelling$performance
"shelling"



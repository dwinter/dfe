#' dfe.
#'
#' @name dfe
#' @docType package
NULL
#' Fitness data for the control lines of a mutation accumulation experiment
#' 
#' These data represent the mean fitness of the ancestral lines used in a
#' mutation accumulation. The most common use of this data will be to estimate
#' the 'experimental variance' (Ve) of the the fitness assay
#' \code{var(ma_control)}
#'
#' @format a vector of floating point fitness values
"ma_control"

#' Fitness data for a set of mutation accumulation lines
#' 
#' These data represent the mean fitness of a set of mutation accumulation 
#' lines. These data can be used to calculate the change in fitness of these i
#' lines relative to  mean fitness of the control data: \code{ma_experimental -
#' mean(ma_control)}
#' @format a vector of floating point fitness values
"ma_experimental"

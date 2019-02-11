#' bioacoustics: detect and extract automatically acoustic features in Zero-Crossing files and audio recordings
#'
#' @description
#'
#' bioacoustics contains all the necessary functions to read Zero-Crossing files and audio recordings of various formats,
#' filter noisy files, display audio signals, detect and extract automatically acoustic features
#' for further analysis such as species identification based on classification of animal vocalizations.
#'
#' @details
#'
#' bioacoustics is subdivided into three main components:
#'
#' \itemize{
#' \item Read, write and manipulate acoustic recordings.
#' \item Display what's inside acoustic recordings, whether to plot or just extract metadata.
#' \item Analyse audio recordings in batch in search of specific vocalizations and extract acoustic features.
#' }
#'
#' To learn more about bioacoustics, start with the introduction vignette:
#' `vignette("introduction", package = "bioacoustics")`
#'
#' @useDynLib bioacoustics, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

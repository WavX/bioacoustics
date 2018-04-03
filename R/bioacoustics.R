#' bioacoustics: detect and extract automatically acoustic features in Zero-Crossing files and audio recordings
#'
#'@description
#'
#' bioacoustics contains all the necessary functions to read Zero-Crossing files and audio recordings of various formats,
#' filter noisy files, display audio signals, detect and extract automatically acoustic features
#' for further analysis such as species identification based on classification of animal vocalizations.
#'
#'@details
#'
#' bioacoustics is subdivided into three main components:
#'
#' \itemize{
#' \item Read, extract data (not yet implemented), display, and write Zero-Crossing files.
#' \item Stand-alone tools to display, convert or resample MP3, WAV, and WAC files.
#' \item Read MP3, WAV or WAC files, filter, and extract automatically acoustic features.
#' }
#'
#' To learn more about bioacoustics, start with the vignette:
#' `browseVignettes(package = "bioacoustics")`
#'
#' @useDynLib bioacoustics, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

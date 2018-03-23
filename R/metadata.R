#' Extract metadata
#'
#' @param x an object for which metadata will be extracted
#'
#' @export
#' @rdname metadata

metadata <- function(x)
{
  UseMethod("metadata", x)
}

#' @export
#' @rdname metadata
metadata.threshold_detection <- function(x)
{
  return(x$metadata)
}

#' @export
#' @rdname metadata
metadata.blob_detection <- function(x)
{
  return(x$metadata)
}

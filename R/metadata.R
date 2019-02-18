#' Extract metadata
#'
#' @param x an object for which metadata will be extracted
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @rdname metadata

metadata <- function(x, ...)
{
  UseMethod("metadata", x)
}

#' @param file_type type of file to read metadata from. Wav and Zero-Crossing files are currently supported.
#' @export
#' @rdname metadata
metadata.character <- function(x, file_type = c(file_type_guess(x), "wav", "zc"), ...)
{
  if (x == "wav")
  {
    return( guano_md(x) )
  }
  else if (x == "zc")
  {
    return( read_zc(x)$metadata )
  }
  else
  {
    stop("Not implemented")
  }
}

#' @export
#' @rdname metadata
metadata.blob_detection <- function(x, ...)
{
  return(x$metadata)
}

#' @export
#' @rdname metadata
metadata.threshold_detection <- function(x, ...)
{
  return(x$metadata)
}

#' Extract metadata from Zero-Crossing files
#'
#' @export
#' @rdname metadata
metadata.zc <- function(x, ...)
{
  return(x$metadata)
}

#' Extract metadata from a Wave object
#'
#' @export
#' @rdname metadata
metadata.Wave <- function(x, ...)
{
  md <- attr(x, "metadata")

  if(is.null(md))
    md <- list()

  md$file$sample_rate <- slot(x, "samp.rate")
  md$file$bit_depth <- slot(x, "bit")
  md
}


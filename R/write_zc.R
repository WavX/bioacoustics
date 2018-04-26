#' Write Zero-Crossing files
#'
#' Write Zero-Crossing files (.zc, .#)
#'
#' @param zc an object of class 'zc'.
#'
#' @param filename path or connection to write.
#'
#' @export
#'
#' @examples
#' data(zc)
#' filename <- tempfile()
#' write_zc(zc, filename = filename)
#'
#' @rdname write_zc
#'

write_zc <- function(zc, filename)
{
  con <- file(filename, open = "wb")
  on.exit(close(con))
  raw_data <- as.raw(zc$raw_data)
  raw_data[7:281] <- update_header.zc(zc)
  writeBin(raw_data, con = con, size = 1)
}


update_header.zc <- function(zc)
{
  lg <- list(8, 8, 40, 50, 16, 74, 79)

  if (zc$FTYPE == 132)
    zc$metadata$DATE <- ""

  unname(
    unlist(
      mapply(
        zc$metadata[c("TAPE", "DATE", "LOC", "SPECIES", "SPEC", "NOTE", "NOTE1")],
        lg = lg,
        SIMPLIFY = FALSE,
        FUN = function(w, lg)
        {
          if (is.null(w)) w <- ""
          charToRaw(stringr::str_pad(w, width = lg, side = "right"))
        }
      )
    )
  )
}


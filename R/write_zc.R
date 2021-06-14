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
  raw <- as.raw(zc$data$raw)
  raw[7:281] <- update_header.zc(zc$metadata$file)
  writeBin(raw, con = con, size = 1)
}


update_header.zc <- function(md_file)
{
  lg <- list(8, 8, 40, 50, 16, 74, 79)

  if (md_file$FTYPE == 132)
    md_file$DATE <- ""

  unname(
    unlist(
      mapply(
        md_file[c("TAPE", "DATE", "LOC", "SPECIES", "SPEC", "NOTE", "NOTE1")],
        lg = lg,
        SIMPLIFY = FALSE,
        FUN = function(w, lg)
        {
          if (is.null(w)) w <- ""
          charToRaw(substr(stringr::str_pad(w, width = lg, side = "right"), 1, lg))
        }
      )
    )
  )
}


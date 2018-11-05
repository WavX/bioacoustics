#' Convert MP3 to WAV
#'
#' Convert an MP3 file to a Wave file
#'
#' @param file path to a MP3 file.
#'
#' @param output_dir where to save the converted Wave file.
#' The Wave file is saved by default to the MP3 file location.
#'
#' @param delete delete the original MP3 file ?
#'
#' @export
#'
#' @rdname mp3_to_wav
#'

mp3_to_wav <- function(file, output_dir = dirname(file), delete = FALSE)
{
  mp3 <- tuneR::readMP3(file)

  extensible <- if (slot(mp3, "samp.rate") > 44100) TRUE else FALSE

  tuneR::writeWave(
    object = mp3,
    file = file.path(output_dir, basename(paste0(tools::file_path_sans_ext(file), ".wav"))),
    extensible = extensible
  )

  if (delete)
    file.remove(file)
}


#' Rotate 90° clockwise
#'
#' Rotate a matrix 90° clockwise
#'
#' @keywords internal
#'
#' @rdname rotate90
#'

rotate90 <- function(m) t(m)[,NROW(m):1]


#' Convert to dB
#'
#' Convert amplitude to decibel (dB) values
#'
#' @param x numeric. Vector of amplitude values (V1).
#'
#' @param ref numeric. Reference value (V0) to calculate the ratio (V1/V0).
#'
#' @keywords internal
#'
#' @rdname to_dB
#'

to_dB <- function(x, ref = 1) 20 * log10(x / ref)

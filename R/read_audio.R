#' Decode audio files
#'
#' Read audio files into a \link{Wave} object. WAV, WAC and MP3 files are
#' currently supported.
#'
#' @param filename a \link{Wave}, WAC or MP3 recording containing animal vocalizations.
#'
#' @param time_exp integer. Time expansion factor of the recording.
#' Set to 1 for real-time recording or above for time expanded recording. Default setting is 1.
#'
#' @param from optional. Numeric. Where to start reading the recording, in seconds (s).
#'
#' @param to optional. Numeric. Where to end reading the recording, in seconds (s).
#'
#' @return A \link{Wave} object.
#'
#' @export
#'
#' @importClassesFrom tuneR Wave
#' @importFrom tuneR readWave
#'
#' @examples
#'
#' \donttest{
#' filepath <- system.file("extdata", "recording.wav", package = "bioacoustics")
#' read_audio(filepath)
#' }
#'
#' @rdname read_audio

read_audio <- function(filename, time_exp = 1, from = NULL, to = NULL)
{
  file_format <- tolower(tools::file_ext(filename))

  switch (file_format,
          wav = read_wav,
          mp3 = read_mp3,
          wac = read_wac,
          stop("This file format ('.", file_format, "') is not supported.")
  ) -> read_fun

  if (!is.null(from) && !is.null(to) && from > to)
    stop("'from' > 'to'")

  read_fun(filename = filename, time_exp = time_exp, from = from, to = to)
}

#' Read MP3 files
#'
#' A thin wrapped around \link[tuneR]{readMP3} from the package tuneR.
#'
#' @param filename a MP3 file.
#'
#' @inheritParams read_audio
#'
#' @param ... currently not implemented.
#'
#' @return A \link[tuneR]{Wave} object.
#'
#' @export
#'
#' @importFrom methods slot slot<-
#'
#' @examples
#' \donttest{
#' filepath <- system.file("extdata", "recording.mp3", package = "bioacoustics")
#' read_mp3(filepath)
#' }
#'
#' @rdname read_mp3

read_mp3 <- function(filename, time_exp = 1, ...)
{
  w <- tuneR::readMP3(filename = filename)
  slot(w, "samp.rate") <- slot(w, "samp.rate") * time_exp
  attr(w, "filename") <- basename(filename)

  return(w)
}

#' Read WAC files from Wildlife Acoustics recorders
#'
#' Convert a Wildlife Acoustics' proprietary compressed WAC file into a \link[tuneR]{Wave} object
#'
#' @param filename a WAC file.
#'
#' @inheritParams read_audio
#'
#' @return A \link[tuneR]{Wave} object.
#'
#' @param write_wav optional folder path where WAV files will be written.
#'
#' @param ... currently not implemented.
#'
#' @export
#'
#' @import tools
#'
#' @examples
#' \donttest{
#' filepath <- system.file("extdata", "recording_20170716_230503.wac", package = "bioacoustics")
#' read_wac(filepath)
#' }
#'
#' @rdname read_wac

read_wac <- function(filename, time_exp = 1, write_wav = NULL, ...)
{
  basenm <- basename(filename)

  w <- read_wac_impl(
    normalizePath(filename, winslash = "/", mustWork = TRUE),
    basenm
  )

  w$sample_rate <- w$sample_rate * time_exp

  fpse <- tools::file_path_sans_ext(basenm)

  # Extract the date time info from the filename
  date_time <- strptime(
    substr( fpse, nchar(fpse) - 14, nchar(fpse) ),
    "%Y%m%d_%H%M%S"
  )

  if (is.na(date_time))
    stop("Could not extract date/time info from filename: '", filename, "'")

  if (length(w$trigger) > 0)
  {
    nm <- paste0(
      substr(
        fpse, 1, nchar(fpse) - 15
      ),
      unlist(
        lapply(
          w$trigger,
          function(tg)
          {
            sub(
              pattern = "\\.(?=[^\\.]*$)",
              replacement = "_",
              x = format(date_time + tg / w$sample_rate, "%Y%m%d_%H%M%OS3"),
              perl = TRUE
            )
          }
        )
      ),
      ".wav"
    )

    return(
      setNames(
        lapply(
          1:length(w$trigger),
          function(i)
          {
            with(
              w,
              {
                new("Wave",
                    stereo = length(right) > 0,
                    samp.rate = sample_rate,
                    bit = bits,
                    pcm = TRUE,
                    left = left[[i]],
                    right = if(length(right) > 0) right[[i]] else numeric()
                ) -> wv

                if (!is.null(write_wav))
                  writeWave(wv, file.path(write_wav, nm[[i]]))

                attr(wv, "filename") <- nm[[i]]
                wv
              }
            )
          }
        ),
        nm = nm
      )
    )
  }

  return(
    with(
      w,
      {
        nm <- paste0(
          substr(
            fpse, 1, nchar(fpse) - 15
          ),
          sub(
            pattern = "\\.(?=[^\\.]*$)",
            replacement = "_",
            x = format(date_time, "%Y%m%d_%H%M%OS3"),
            perl = TRUE
          ),
          ".wav"
        )
        new("Wave",
            stereo = length(right) > 0,
            samp.rate = sample_rate,
            bit = bits,
            pcm = TRUE,
            left = left[[1L]],
            right = if(length(right) > 0) right[[1L]] else numeric()
        ) -> wv

        if (!is.null(write_wav))
          writeWave(wv, file.path(write_wav, nm))

        attr(wv, "filename") <- nm
        wv
      }
    )
  )
}


#' Read WAV files
#'
#' A thin wrapped around \link[tuneR]{readWave} from the package tuneR.
#'
#' @param filename a WAV file.
#'
#' @inheritParams read_audio
#'
#' @return A \link[tuneR]{Wave} object.
#'
#' @export
#'
#' @importFrom methods slot slot<-
#'
#' @examples
#' \donttest{
#' filepath <- system.file("extdata", "recording.wav", package = "bioacoustics")
#' read_wav(filepath)
#' }
#'
#' @rdname read_wav

read_wav <- function(filename, time_exp = 1, from = NULL, to = NULL)
{
  tuneR::readWave(
    filename = filename,
    from = ifelse(is.null(from), 0, from * time_exp),
    to = ifelse(is.null(to), Inf, to * time_exp),
    units = "seconds"
  ) -> w
  slot(w, "samp.rate") <- slot(w, "samp.rate") * time_exp
  attr(w, "filename") <- basename(filename)

  return(w)
}

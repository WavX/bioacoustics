#' Amplitude threshold detector above Signal to Noise Ratio (SNR)
#'
#' This function is a modified version of the Bat Bioacoustics freeware developed by Christopher Scott (2012).
#' It combines several detection, filtering and audio feature extraction algorithms.
#'
#' @param wave either a path to a file, or a \link[tuneR]{Wave} object.
#'
#' Audio files will be automatically decoded internally using the function \link{read_audio}.
#'
#' @param threshold integer. Sensitivity of the audio event detection function (peak-picking algorithm) in dB.
#' A threshold value of 14 dB above SNR is recommended. Higher values increase the risk of leaving audio events undetected (false negative).
#' In a noisy recording (low SNR) this sensitivity threshold may be set at 12 dB,
#' but a value below 10 dB is not recommended. Default setting is 14 dB above SNR.
#'
#' @param channel character. Channel to keep for analysis in a stereo recording: 'left' or 'right'.
#' Do not need to be specified for mono recordings, recordings with more than two channels are not
#' yet supported. Default setting is 'left'.
#'
#' @inheritParams read_audio
#'
#' @param min_dur numeric. Minimum duration threshold in milliseconds (ms).
#' Extracted audio events shorter than this threshold are ignored. Default setting is 1.5 ms.
#'
#' @param max_dur numeric. Maximum duration threshold in milliseconds (ms).
#' Extracted audio events longer than this threshold are ignored. The default setting is 80 ms.
#'
#' @param min_TBE numeric. Minimum time window between two audio events in milliseconds (ms). If the time interval between two
#' successive audio events is shorter than this window, they are ignored. The default setting is 20 ms.
#'
#' @param max_TBE numeric. Maximum time window between two audio events in milliseconds (ms). If the time interval between two
#' successive audio events is longer than this window, they are ignored. The default setting is 1000 ms.
#'
#' @param EDG numeric. Exponential Decay Gain from 0 to 1. Sets the degree of temporal masking at the end of each audio event.
#' This filter avoids extracting noise or echoes at the end of the audio event. The default setting is 0.996.
#'
#' @param LPF integer. Low-Pass Filter (Hz). Frequencies above the cutoff are greatly attenuated.
#' Default is set internally at the Nyquist frequency of the recording.
#'
#' @param HPF integer. High-Pass Filter (Hz). Frequencies below the cutoff are greatly attenuated.
#' Default setting is 16000 Hz. A default of 1000 Hz is recommended for most bird vocalizations.
#'
#' @param FFT_size integer. Size of the Fast Fourrier Transform (FFT) window. Default setting is 256.
#'
#' @param FFT_overlap numeric. Percentage of overlap between two FFT windows (from 0 to 1). Default setting is 0.875.
#'
#' @param start_thr integer. Right to left amplitude threshold (dB) for audio event extraction, from the audio event centroid.
#' The last FFT where the amplitude level is equal or above this threshold is considered the start of the audio event.
#' Default setting is 40 dB. 20 dB is recommended for extracting bird vocalizations.
#'
#' @param end_thr integer. Left to right amplitude threshold (dB) for audio event extraction, from the audio event centroid.
#' The last FFT where the amplitude level is equal or above this threshold is considered the end of the audio event.
#' Default setting is 20 dB. 30 dB is recommended for extracting bird vocalizations.
#'
#' @param SNR_thr integer. SNR threshold (dB) at which the extraction of the audio event stops.
#' Default setting is 10 dB. 8 dB is recommended for bird vocalizations.
#'
#' @param angle_thr integer. Angle threshold (°) at which the audio event extraction stops.
#' Default setting is 40°. 125° is recommended for extracting bird vocalizations.
#'
#' @param duration_thr integer. Maximum duration threshold in milliseconds (ms) after which the monitoring of the background noise is resumed.
#' Default setting is 80 ms for bat echolocation calls. A higher threshold value is recommended for extracting bird vocalizations.
#'
#' @param NWS integer. Length of the time window used for background noise estimation in the recording (ms).
#' A longer window size is less sensitive to local variations in the background noise.
#' Default setting is 100 ms.
#'
#' @param KPE numeric. Set the Process Error parameter of the Kalman filter.
#' Default setting is 1e-05.
#'
#' @param KME numeric. Set the Measurement Error parameter of the Kalman filter.
#' Default setting is 1e-05.
#'
#' @param settings logical. \code{TRUE} or \code{FALSE}. Save on a list the parameters set with the threshold_detection function.
#' Default setting is \code{FALSE}.
#'
#' @param acoustic_feat logical. \code{TRUE} or \code{FALSE}. Extracts the acoustic and signal quality parameters from each audio event in a data frame.
#' The sequences of smoothed amplitude (dB) and frequency (Hz) bins of each audio event, temporal values (in ms)
#'  of the beginning and the end of each audio event are also extracted in separate lists. Default setting is \code{TRUE}.
#'
#' @param metadata logical. \code{TRUE} or \code{FALSE}. Extracts on a list the metadata embedded with the Wave file
#' GUANO metadata extraction is not -yet- implemented. Default setting is \code{FALSE}.
#'
#' @param spectro_dir character (path) or \code{NULL}. Generate an HTML page with the spectrograms numbered by order
#' of detection in the recording. Spectrograms are generated as individual .PNG files and stored in the
#' 'spectro_dir/spectrograms' subdirectory. The R working directory is used if \code{spectro_dir} is \code{NULL}.
#' \code{spectro_dir} is set to \code{NULL} by default.
#'
#' @param time_scale numeric. Time resolution of the spectrogram in milliseconds (ms) per pixel (px). Default setting is 0.1 ms for bat echolocation calls.
#' A default of 2 ms/px is recommended for most bird vocalizations.
#'
#' @param ticks either logical or numeric. If \code{TRUE} tickmarks are drawn on the (frequency)
#' y-axis and their positions are computed automatically. If numeric, sets the
#' lower and upper limits of the tickmarks and their interval (in Hz). Default setting is \code{TRUE}.
#'
#' @return an object of class 'bioacoustics_output'.
#'
#' @export
#'
#' @import tools
#' @importClassesFrom tuneR Wave
#' @importFrom grDevices dev.off
#' @importFrom methods slot slot<-
#'
#' @examples
#' data(myotis)
#' Output <- threshold_detection(myotis, time_exp = 10, HPF = 16000, LPF = 200000)
#' Output$data
#'
#' @rdname threshold_detection
#'

threshold_detection <- function(wave,
                                threshold = 14,
                                channel = 'left',
                                time_exp = 1,
                                min_dur = 1.5,
                                max_dur = 80,
                                min_TBE = 20,
                                max_TBE = 1000,
                                EDG = 0.996,
                                LPF,
                                HPF = 16000,
                                FFT_size = 256,
                                FFT_overlap = .875,
                                start_thr = 40,
                                end_thr = 20,
                                SNR_thr = 10,
                                angle_thr = 40,
                                duration_thr = 80,
                                NWS = 100,
                                KPE = 0.00001,
                                KME = 0.00001,
                                settings = FALSE,
                                acoustic_feat = TRUE,
                                metadata = FALSE,
                                spectro_dir = NULL,
                                time_scale = 0.1,
                                ticks = TRUE)

{
  if (is.character(wave))
    wave <- read_audio(wave)

  filepath <- attr(wave, 'filepath')

  filename <-
    if (!is.null(attr(wave, 'filename')))
      attr(wave, 'filename')
    else if (!is.null(filepath))
      basename(filepath)
    else
      NA_character_

  sample_rate <- slot(wave, 'samp.rate') * time_exp
  bit_depth <- slot(wave, 'bit')

  if(missing(LPF))
    LPF <- sample_rate / 2
  else
    LPF <- min(LPF, sample_rate / 2)

  event_list <- threshold_detection_impl(
    audio_samples = slot(wave, channel <- ifelse(slot(wave, 'stereo'), channel, 'left')),
    sample_rate = sample_rate,
    threshold = threshold,
    min_d = min_dur,
    max_d = max_dur,
    min_TBE = min_TBE,
    max_TBE = max_TBE,
    EDG = EDG,
    LPF = LPF,
    HPF = HPF,
    FFT_size = FFT_size,
    FFT_overlap = FFT_overlap,
    dur_t = duration_thr,
    snr_t = SNR_thr,
    angl_t = angle_thr,
    start_t = start_thr,
    end_t = end_thr,
    NWS = NWS,
    KPE = KPE,
    KME = KME
  )

  if (length(event_list) == 0)
  {
    message(
      "No audio events found",
      if (!is.na(filename)) paste0(" for file '", basename(filename), "'")
    )

    output <- list(data = NULL)
  }
  else
  {
    if ( !is.null(spectro_dir) )
    {
      if ( !dir.exists(file.path(spectro_dir, 'spectrograms')) )
        stopifnot(dir.create(file.path(spectro_dir, 'spectrograms'), recursive = TRUE)) # Stop if directory creation fails

      spectro_dir <- tools::file_path_as_absolute(spectro_dir)

      if (is.logical(ticks))
      {
        if (ticks)
        {
          ticks_at <- seq(0, FFT_size, length.out = 12) / FFT_size
        }
      }
      else if (is.numeric(ticks))
      {
        ticks_at <- seq(ticks[1L], ticks[2L], ticks[3L]) / (sample_rate / 2)
        ticks <- TRUE
      }

      bare_name <- if (!is.na(filename)) tools::file_path_sans_ext(filename)
      html_file <- paste0('spectrograms--', bare_name, format(Sys.time(), '--%y%m%d--%H%M%S'), '.html')

      tags <- htmltools::tags

      fft_step <- floor(FFT_size * (1 - FFT_overlap))

      html_doc <- tags$html(
        tags$head(
          tags$title(html_file)
        ),
        tags$body(
          htmltools::tagList(
            lapply(
              1:length(event_list[['event_start']]),
              function(i)
              {
                offset <- .5 / (1 - FFT_overlap)

                start <- max(event_list[['event_start']][[i]] - offset * fft_step, 1)
                end <- min(length(slot(wave, channel)), event_list[['event_end']][[i]] + offset * fft_step)

                offset <- floor( c(event_list[['event_start']][[i]] - start, end - event_list[['event_end']][[i]]) / fft_step )

                LPF <- min(LPF, sample_rate / 2)
                LPF_bin <- min(ceiling(LPF * FFT_size / sample_rate), FFT_size / 2 - 1)

                spec <- .fspec_impl(
                  audio_samples = slot(wave, channel)[start:end], FFT_size,
                  FFT_overlap, 'blackman7', HPF_bin = 0, LPF_bin = LPF_bin,
                  FLL_bin = 0, FUL_bin = LPF_bin, rotate = TRUE
                )
                png_file <- file.path('spectrograms', paste0(format(Sys.time(), '%y%m%d--%H%M%S--'), bare_name, '--', i, '.png'))
                png(file.path(spectro_dir, png_file), width = (((1 - FFT_overlap) * FFT_size) / sample_rate * 1000) * nrow(spec) / time_scale)
                par(mar = rep_len(0L,4L), oma = c(.1,.1,.1,.1))

                .spectro(data = to_dB(spec), colors = gray.colors(25, 1, 0))
                X <- (offset[1L] + 1):(nrow(spec) - offset[2L]) / nrow(spec)
                # lines(I(event_list[['freq_track']][[i]] * 2 / sample_rate) ~ X, col = 'red', lwd = 2)
                lines(I(event_list[['freq_track']][[i]] / LPF) ~ X, col = 'red', lwd = 2)
                lgd <- legend('topright', legend = NA, inset = 0, box.col = NA)
                text(x = lgd$rect$left + min(lgd$rect$w, 1) / 2L, y = lgd$text$y, labels = i, adj = .5)
                if (ticks)
                {
                  axis(2, at = ticks_at, tcl = .4, col = NA, col.ticks = 'black', labels = FALSE)
                }
                box()
                dev.off()
                tags$img(src = png_file)
              }
            )
          )
        )
      )

      htmltools::save_html(html = html_doc, file = file.path(spectro_dir, html_file))
      utils::browseURL(file.path(spectro_dir, html_file)) # Open html page on favorite browser
    }

    # Convert amp_track to dB
    event_list[['amp_track']] <- lapply(event_list[['amp_track']], to_dB, ref = 2 ^ (bit_depth - 1))

    if (!acoustic_feat)
    {
      event_list[['amp_track']] <- NULL
      event_list[['freq_track']] <- NULL
      event_list$event_data <- event_list$event_data[, c('starting_time', 'duration')]
    }

    # Add info about filename if available
    event_list$event_data <- cbind(
      data.frame(
        filename = basename(filename),
        stringsAsFactors = FALSE
      ),
      event_list$event_data
    )

    # Convert event_start and event_end to vector
    event_list$event_start <- unlist(event_list$event_start)
    event_list$event_end <- unlist(event_list$event_end)

    output <- list(data = event_list)
  }

  if (metadata)
  {
    if (!is.null(attr(wave, "metadata")))
    {
      output$metadata <- metadata(wave)
    }
    else
    {
      if (!is.null(filepath))
      {
        output$metadata$file$sample_rate <- sample_rate
        output$metadata$file$bit_depth <- bit_depth

        if (length(guano_md <- guano_md(filepath)) > 0)
          output$metadata$guano <- guano_md
      }
    }
  }

  if (settings)
  {
    output$metadata$settings_threshold_detection <- list(
      threshold = threshold,
      channel = channel,
      min_dur = min_dur,
      max_dur = max_dur,
      min_TBE = min_TBE,
      max_TBE = max_TBE,
      EDG = EDG,
      LPF = LPF,
      HPF = HPF,
      FFT_size = FFT_size,
      FFT_overlap = FFT_overlap,
      start_thr = start_thr,
      end_thr = end_thr,
      SNR_thr = SNR_thr,
      angle_thr = angle_thr,
      duration_thr = duration_thr,
      NWS = NWS,
      KPE = KPE,
      KME = KME,
      stringsAsFactors = FALSE
    )
  }

  class(output) <- "threshold_detection"

  return(output)
}

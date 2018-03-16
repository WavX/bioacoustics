#' Blob detection of a region into a spectrogram representation of the recording
#'
#' This function is a modified version of the Bat classify software developed by Christopher Scott (2014).
#' It combines several detection, filtering and audio feature extraction processes.
#'
#' @inheritParams threshold_detection
#'
#'@param min_area integer. Minimum area threshold in number of pixels.
#' Extracted segments with an area shorter than this threshold are discarded.
#' Default setting is 40 pixels.
#'
#' @param blur integer. Gaussian smoothing function for blurring the spectrogram of the audio event to reduce image noise.
#' Default setting is 2.
#'
#' @param contrast_boost integer.Edge contrast enhancement filter of the spectrogram of the audio event to improve its apparent sharpness.
#'Default setting is 20.
#'
#' @param bg_substract integer. Foreground extraction with a mean filter applied on the spectrogram of the audio even for image denoising.
#' Default setting is 20.
#'
#' @export
#'
#' @importClassesFrom tuneR Wave
#' @importFrom graphics axis box legend lines par plot text
#' @importFrom grDevices png
#' @importFrom methods is slot
#' @importFrom stats setNames
#'
#' @examples
#' data(myotis)
#' Output <- blob_detection(myotis, time_exp = 10, contrast_boost = 30, bg_substract = 30)
#' Output$event_data
#'
#' @rdname blob_detection
#'

blob_detection <- function(wave,
                           channel = 'left',
                           time_exp = 1,
                           min_dur = 1.5,
                           max_dur = 80,
                           min_area = 40,
                           TBE = 20,
                           EDG = .9,
                           LPF,
                           HPF = 16000,
                           FFT_size = 256,
                           FFT_overlap = 0.875,
                           blur = 2,
                           bg_substract = 20,
                           contrast_boost = 20,
                           settings = FALSE,
                           acoustic_feat = TRUE,
                           metadata = FALSE,
                           spectro_dir = NULL,
                           time_scale = 0.1,
                           ticks = TRUE)
{
  if (!is(wave, 'Wave'))
    wave <- read_audio(wave)

  filename <- attr(wave, 'filename')

  sample_rate <- slot(wave, 'samp.rate') * time_exp
  n_bits <- slot(wave, 'bit')

  if(missing(LPF))
    LPF <- sample_rate / 2

  if(missing(LPF))
    LPF <- sample_rate / 2
  else
    LPF <- min(LPF, sample_rate / 2)

  blobs <- blob_detection_impl(audio_samples = slot(wave, channel <- ifelse(slot(wave, 'stereo'), channel, 'left')),
                               sample_rate = sample_rate,
                               LPF = LPF,
                               HPF = HPF,
                               TBE = TBE,
                               EDG = EDG,
                               FFT_size = FFT_size,
                               FFT_overlap = FFT_overlap,
                               min_d = min_dur,
                               max_d = max_dur,
                               area = min_area,
                               blur_f = blur,
                               bg_substract = bg_substract,
                               boost = contrast_boost)

  if (length(blobs) == 0)
  {
    message(
      "No audio events found",
      if (!is.null(filename)) paste0(" for file '", basename(filename), "'")
    )
  }
  else
  {
    audio_events <- list( event_data = blobs[[1L]] )

    # features <- lapply(blobs[[2]], function(blob)
    # {
    #   blob <- log10(20.0 * blob + 1) # log compress
    #   height <- nrow(blob)
    #   power_spectrum <- vector('Numeric', height)
    #   histogram <- vector('Numeric', height * 2L)
    #   temporal_envelope <- colSums(blob)
    #   prev_peak <- 0
    #
    #   for(i in 1:NCOL(blob))
    #   {
    #     if(temporal_envelope[i] > .Machine$double.eps)
    #     {
    #       # MaskSpectrum
    #       wm <- which.max(blob[,i])
    #       if(wm != height)
    #       {
    #         distance <- height - wm +1L
    #         if(distance > 8L)
    #         {
    #           blob[(wm + 8L):height, i] <- 0
    #         }
    #
    #         distance <- wm - 1L
    #         if(distance > 8L)
    #         {
    #           blob[1:(wm - 9L), i] <- 0
    #         }
    #       }
    #     }
    #
    #     cur_peak <- ifelse(sum(blob[2:(height - 1L),i]) > 0,
    #                        round(sum(blob[2:(height - 1L),i] * (1:(height - 2L))) / sum(blob[2:(height - 1L),i])), ## round(centroid)
    #                        0)
    #
    #     if (cur_peak > 0 && prev_peak > 0)
    #     {
    #       bin <- (cur_peak - prev_peak) + height + 1L
    #       histogram[bin] <- histogram[bin] + blob[cur_peak + 1L,i]
    #       histogram[bin - 1L] <- histogram[bin - 1L] + blob[cur_peak,i]
    #       histogram[bin + 1L] <- histogram[bin + 1L] + blob[cur_peak + 2L,i]
    #     }
    #     prev_peak <- cur_peak
    #   }
    #
    #   area <- sum(blob[2:(height - 1L),] > 0.00001)
    #   power_spectrum[2:(height - 1L)] <- rowSums(blob[2:(height - 1L),])
    #
    #   ## // frequency
    #     ## Moments
    #     f_x <- power_spectrum/max(sum(power_spectrum), .Machine$double.eps)
    #     f_centroid <- drop(f_x %*% 1:height)
    #     f_delta <- (1:height - f_centroid)
    #     f_bandwith <- sqrt(sum(f_delta ^ 2L * f_x))
    #     f_skew <- ifelse(f_bandwith > .Machine$double.eps, sum((f_delta ^ 3L) * f_x) / (f_bandwith ^ 3L), 0)
    #     f_kurtosis <- ifelse(f_bandwith > .Machine$double.eps, sum((f_delta ^ 4L) * f_x) / (f_bandwith ^ 4L), 3) - 3
    #
    #     ## Q-value (centroid frequency divided by full width at half-maximum)
    #     Q <- f_centroid / ifelse(f_bandwith > 1, f_bandwith, 1)
    #
    #     ## GiniImpurity
    #     f_gini <- 1 - sum((power_spectrum / sum(power_spectrum)) ^ 2)
    #
    #     ## Quantiles
    #     partial_sum <- cumsum(power_spectrum)
    #     quantiles <- unlist(lapply(c(0.025, 0.25, 0.5, 0.75, 0.975), function(q) setNames(match(TRUE, partial_sum >= sum(power_spectrum) * q) - 1L, nm = paste0('Quant_',q*100,'%'))))
    #
    #   ## // temporal
    #     ## Moments
    #     t_x <- temporal_envelope/max(sum(temporal_envelope), .Machine$double.eps)
    #     t_centroid <- drop(t_x %*% seq_along(temporal_envelope))
    #     t_delta <- (seq_along(temporal_envelope) - t_centroid)
    #     t_bandwith <- sqrt(sum(t_delta ^ 2L * t_x))
    #     t_skew <- ifelse(t_bandwith > .Machine$double.eps, sum((t_delta ^ 3L) * t_x) / (t_bandwith ^ 3L), 0)
    #     t_kurtosis <- ifelse(t_bandwith > .Machine$double.eps, sum((t_delta ^ 4L) * t_x) / (t_bandwith ^ 4L), 3) - 3
    #
    #     ## GiniImpurity
    #     t_gini <- 1 - sum((temporal_envelope / sum(temporal_envelope)) ^ 2)
    #
    #   ## // gradient
    #     ## Moments
    #     g_x <- histogram/max(sum(histogram), .Machine$double.eps)
    #     g_centroid <- drop(g_x %*% seq_along(histogram))
    #     g_delta <- (seq_along(histogram) - g_centroid)
    #     g_bandwith <- sqrt(sum(g_delta ^ 2L * g_x))
    #     g_skew <- ifelse(g_bandwith > .Machine$double.eps, sum((g_delta ^ 3L) * g_x) / (g_bandwith ^ 3L), 0)
    #     g_kurtosis <- ifelse(g_bandwith > .Machine$double.eps, sum((g_delta ^ 4L) * g_x) / (g_bandwith ^ 4L), 3) - 3
    #
    #     ## GiniImpurity
    #     g_gini <- 1 - sum((histogram / sum(histogram)) ^ 2)
    #
    #   c(area = area,
    #     freq_centroid = f_centroid,
    #     freq_bandwith = f_bandwith,
    #     freq_skew = f_skew,
    #     freq_kurtosis = f_kurtosis,
    #     Q = Q,
    #     freq_GiniImpurity = f_gini,
    #     quantiles,
    #     '95%_CI' = quantiles[5L] - quantiles[1L],
    #     IQR = quantiles[4L] - quantiles[2L],
    #     temporal_centroid = t_centroid,
    #     temporal_bandwith = t_bandwith,
    #     temporal_skew = t_skew,
    #     temporal_kurtosis = t_kurtosis,
    #     temporal_GiniImpurity = t_gini,
    #     grad_centroid = g_centroid,
    #     grad_bandwith = g_bandwith,
    #     grad_skew = g_skew,
    #     grad_kurtosis = g_kurtosis,
    #     grad_GiniImpurity = g_gini)
    # })

    # if (acoustic_feat)
    #   audio_events$event_data <- cbind(audio_events$event_data, as.data.frame(do.call('rbind', features)))

    if (!acoustic_feat)
      audio_events$event_data <- audio_events$event_data[, c('starting_time', 'duration')]

    if ( !is.null(spectro_dir) )
    {
      if ( !dir.exists(file.path(spectro_dir, 'spectrograms')) )
        stopifnot(dir.create(file.path(spectro_dir, 'spectrograms'), recursive = TRUE)) # Stop if directory creation fails

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

      bare_name <- if (!is.null(filename)) tools::file_path_sans_ext(filename)
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
              1:length(blobs[[2]]),
              function(i)
              {
                png_file <- file.path('spectrograms', paste0(format(Sys.time(), '%y%m%d--%H%M%S--'), bare_name, '--', i, '.png'))
                png(file.path(spectro_dir, png_file), width = (((1 - FFT_overlap) * FFT_size) / sample_rate * 1000) * ncol(blobs[[2]][[i]]) / time_scale)
                par(mar = rep_len(0L,4L), oma = c(.1,.1,.1,.1))
                .spectro(data = to_dB(rotate90(blobs[[2]][[i]]), ref = 2^(n_bits-1)), colors = gray.colors(25, 1, 0))
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

    if (metadata)
    {
      audio_events$metadata <- list(
        sample_rate = sample_rate,
        n_bits = n_bits
      )
    }

    if (settings)
    {
      audio_events$settings$blob_detection <- data.frame(
        channel = channel,
        min_dur = min_dur,
        max_dur = max_dur,
        min_area = min_area,
        LPF = LPF,
        HPF = HPF,
        EDG = EDG,
        FFT_size = FFT_size,
        FFT_overlap = FFT_overlap,
        blur = blur,
        bg_substract = bg_substract,
        EDG = EDG,
        contrast_boost = contrast_boost,
        stringsAsFactors = FALSE
      )
    }

    # Add info about filename if x is a file path
    if (!is.null(filename)) audio_events$event_data <- cbind(data.frame(filename = filename), audio_events$event_data)

    return(audio_events)
  }
}

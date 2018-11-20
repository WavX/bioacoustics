#' Generate spectrograms
#'
#' This function returns the spectrographic representation of a time wave in decibels (dB) using the Fast Fourier transform (FFT).
#'
#' @param wave a \link[tuneR]{Wave} object.
#'
#' @param channel character. Channel to keep for analysis in a stereo recording: "left" or "right". Default setting is left.
#'
#' @param FFT_size integer. Size of the Fast Fourrier Transform (FFT) window. Default setting is 256.
#'
#' @param FFT_overlap numeric. Percentage of overlap between two FFT windows (from 0 to 1). Default setting is 0.875.
#'
#' @param FFT_win character. Specify the type of FFT window: "hann", "blackman4", or "blackman7".
#' Default setting is "hann".
#'
#' @param LPF integer. Low-Pass Filter (Hz). Frequencies above the cutoff are greatly attenuated.
#' Default setting is the Nyquist frequency of the recording.
#'
#' @param HPF integer. High-Pass Filter (Hz). Frequencies below the cutoff are greatly attenuated.
#' Default setting is 0 Hz.
#'
#' @param tlim numeric. Specify the time limits on the X-axis in seconds (s).
#' Default setting is \code{NULL}, i.e no time limits.
#'
#' @param flim numeric. Specify the frequency limits on the Y-axis in Hz. Default
#' setting is \code{NULL}, i.e. frequency limits are equal to \code{c(0, LPF)}.
#'
#' @param rotate logical. Should the matrix be rotated 90° counter clockwise ?
#' Default setting is \code{FALSE}.
#'
#' @return A matrix of decibel (dB) values in the time / frequency domain.
#'
#' @examples
#' data(myotis)
#' image(fspec(myotis, tlim = c(1, 2), rotate = TRUE))
#'
#' @export
#'
#' @rdname fspec

fspec <- function(wave,
                  channel = "left",
                  FFT_size = 256,
                  FFT_overlap = .875,
                  FFT_win = "hann",
                  LPF,
                  HPF = 0,
                  tlim = NULL,
                  flim = NULL,
                  rotate = FALSE)
{
  sample_rate <- slot(wave, "samp.rate")
  audio_samples <- slot(wave, channel)

  if (!is.null(tlim))
  {
    if (length(tlim) != 2)
      stop("'tlim' should be of length 2")

    from <- max(1, floor(tlim[1L] * sample_rate))
    to <- min(length(audio_samples), ceiling(tlim[2L] * sample_rate))
    audio_samples <- audio_samples[from:to]
  }

  if(missing(LPF))
    LPF <- sample_rate / 2
  else
    LPF <- min(LPF, sample_rate / 2)

  HPF_bin <- max(floor(HPF * FFT_size / sample_rate), 0)
  LPF_bin <- min(ceiling(LPF * FFT_size / sample_rate), FFT_size / 2 - 1)

  if (is.null(flim))
  {
    FLL_bin = 0
    FUL_bin = LPF_bin
  }
  else
  {
    if (length(flim) != 2)
      stop("'flim' should be of length 2")

    if (flim[2L] > sample_rate / 2)
    {
      flim[2L] <- sample_rate / 2
      warning("'flim[2]' was above the Nyquist, reset to the Nyquist")
    }

    FLL_bin = max(floor(flim[1L] * FFT_size / sample_rate), 0);

    if (flim[1L] > HPF)
    {
      HPF_bin = FLL_bin;
    }

    FUL_bin = min(ceiling(flim[2L] * FFT_size / sample_rate), FFT_size / 2 - 1);

    if (flim[2L] < LPF)
    {
      LPF_bin = FUL_bin;
    }
  }

  bit_depth <- slot(wave, "bit")

  to_dB(
    .fspec_impl(
      audio_samples, FFT_size, FFT_overlap, FFT_win,
      HPF_bin, LPF_bin, FLL_bin, FUL_bin, rotate = rotate
    ),
    ref = 2^(bit_depth - 1) # dBFS scale
  )
}

#' @importFrom graphics image plot
#'

.spectro <- function(data, colors)
{
  X <- seq(0, 1, length.out = nrow(data))
  Y <- seq(0, 1, length.out = ncol(data))

  image(x = X, y = Y, z = data, col = colors, useRaster = TRUE, axes = FALSE)
}

#' Plot a spectrogram
#'
#'
#' @inheritParams fspec
#'
#' @param ticks_y numeric. Whether tickmarks should be drawn on the frequency Y-axis or not.
#' The lower and upper bounds of the tickmarks and their intervals (in Hz) has to be specified.
#' Default setting is \code{NULL}.
#'
#' @param col set the colors for the amplitude scale (dB) of the spectrogram.
#'
#' @export
#'
#' @importFrom grDevices gray.colors
#'
#' @examples
#' data(myotis)
#' spectro(myotis, tlim = c(1, 2))
#'
#' @rdname spectro

spectro <- function(wave,
                    channel = "left",
                    FFT_size = 256,
                    FFT_overlap = .875,
                    FFT_win = "hann",
                    LPF,
                    HPF = 0,
                    tlim = NULL,
                    flim = NULL,
                    ticks_y = NULL,
                    col = gray.colors(25, 1, 0))
{
  op <- par(mar = rep_len(0L, 4L), oma = rep_len(0L, 4L))
  on.exit(par(op))

  sample_rate <- slot(wave, "samp.rate")

  if(missing(LPF))
    LPF <- sample_rate / 2
  else
    LPF <- min(LPF, sample_rate / 2)

  spec <- fspec(wave = wave,
                channel = channel,
                FFT_size = FFT_size,
                FFT_overlap = FFT_overlap,
                FFT_win = "hann",
                LPF = LPF,
                HPF = HPF,
                tlim = tlim,
                flim = flim,
                rotate = TRUE)

  .spectro(spec, col)

  tcl <- .4

  if (!is.null(ticks_y))
  {
    if (length(ticks_y) != 3)
      stop("'ticks_y' should be of length 3")

    axis(2, at = seq(ticks_y[1L], ticks_y[2L], ticks_y[3L]) / (sample_rate / 2), tcl = tcl, col = NA, col.ticks = "black", labels = FALSE)
  }

  box()
}





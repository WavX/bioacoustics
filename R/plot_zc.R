#' Generate spectrogram for Zero-Crossing files
#'
#' Generate spectrogram for Zero-Crossing files.
#'
#' @param x an object of class 'zc'.
#'
#' @param LPF numeric. Low-Pass Filter (Hz). Frequencies above the cutoff are
#' greatly attenuated. Default is set to 125000 Hz.
#'
#' @param HPF numeric. High-Pass Filter (Hz). Frequencies below the cutoff are
#' greatly attenuated. Default setting is 16000 Hz.
#'
#' @param tlim numeric. Time limits of the plot in seconds (s). Default setting
#' is set to \code{c(0, Inf)}.
#'
#' @param flim numeric. Frequency limits of plot in Hz. Default setting is set
#' to \code{c(HPF, LPF)}
#'
#' @param ybar should horizontal scale bars be plotted. Default is \code{TRUE}.
#'
#' @param ybar.lty line type of the horizontal scale bars.
#'
#' @param ybar.col color of the horizontal scale bars.
#'
#' @param dot.size dot size.
#'
#' @param dot.col dot color.
#'
#' @param ... not currently implemented.
#'
#' @export
#'
#' @importFrom graphics abline plot
#'
#' @examples
#' data(zc)
#' plot_zc(zc)

plot_zc <- function(x, LPF = 125000, HPF = 16000, tlim = c(0, Inf),
                    flim = c(HPF, LPF), ybar = TRUE, ybar.lty = 2,
                    ybar.col = "gray", dot.size = .3, dot.col = "red",  ...)
{
  x$data$freq_data[x$data$freq_data >= LPF] <- NA
  x$data$freq_data[x$data$freq_data <= HPF] <- NA

  tlim[2L] <- min(tlim[2L], max(x$data$time_data) / 1e06)

  plot(x = x$data$time_data / 1e06, y = x$data$freq_data, xlim = tlim, ylim = flim,
       xlab = "Time (Seconds)", ylab = "Frequency (Hz)",
       main = x$metadata$SPECIES, col = dot.col, cex = dot.size, pch = 16)

  if (ybar)
  {
    # Find out y-tick positions, but for now...
    nLines <- (flim[2L] + 10000) %/% 10000
    abline(h = (1:nLines - 1) * 10000, lty = ybar.lty, col = ybar.col)
  }
}

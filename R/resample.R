#' Resample
#'
#' Resample a Wave object to a given sampling frequency.
#'
#' @param wave a \link[tuneR]{Wave} object.
#' @param to target frequency in Hz
#' @return a \link[tuneR]{Wave} object.
#'
#' @export
#'
#' @importFrom methods slot
#' @importClassesFrom tuneR Wave
#'
#' @examples
#' data(myotis)
#' myotis_192 <- resample(myotis, to = 192000)
#' spectro(myotis_192, tlim = c(1, 1.5))
#'
#' @rdname resample
#'

resample <- function(wave, to)
{
  stopifnot(is(wave, "Wave"))
  from <- slot(wave, "samp.rate")
  r <- to / from

  slot(wave, "left") <- resample_impl(slot(wave, "left"), r)
  if ( slot(wave, "stereo") ) slot(wave, "right") <- resample_impl(slot(wave, "right"), r)

  slot(wave, "samp.rate") <- to
  return(wave)
}

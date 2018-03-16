#-------------------------------------------------------------------------------
# Copyright (C) 2013 Peter D. Wilson (peter@peterwilson.id.au)
# Copyright (C) 2017-2018 WavX, inc. (www.wavx.ca)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See http://www.gnu.org/licenses/ for a copy of the license.
#
# Originally written in C by Chris Corben (chris@hoarybat.com)
# Available here: http://users.lmi.net/corben/loader.txt
# More here: http://users.lmi.net/corben/fileform.htm
# Ported in R by Peter Wilson (peter@peterwilson.id.au)
# Available here: http://www.peterwilson.id.au/Rcode/AnabatTools.R
# Further improvements and additions by WavX inc. (info@wavx.ca)
#-------------------------------------------------------------------------------

#' Read Zero-Crossing files
#'
#' Read Zero-Crossing files (.zc, .#) from various bat recorders
#'
#' @param filename a Zero-Crossing file.
#'
#' @return an object of class 'zc'.
#'
#' @export
#'
#' @import moments
#' @import stringr
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' zc <- read_zc("filename")
#' }
#'
#' @rdname read_zc
#'

read_zc <- function(filename)
{
  zc <- list(raw_data = readBin(filename, what = "integer", n = 16384, size = 1, signed = FALSE))

  nBytes <- length(zc$raw_data)

  zc$FTYPE <- FTYPE <- zc$raw_data[4]

  if (FTYPE < 0x81 || FTYPE > 0x84)
    stop("File type must be betwen 129 and 132.")

  # Set Pointer to parameter table:
  p <- zc$raw_data[1] + zc$raw_data[2] * 256 + 1L

  if(p != 0x11a + 1)
    stop("File is corrupted!")

  # Read text header
  zc$metadata <- read_header.zc(zc$raw_data)

  # Read settings
  zc$metadata <- c(zc$metadata, read_settings.zc(p, zc$raw_data))

  if (zc$metadata$RES1 > 60000 || zc$metadata$RES1 < 10000) stop("File is corrupted!")
  if (zc$raw_data[p + 4] > 64 || zc$raw_data[p + 4] < 1) stop("File is corrupted!")

  # Status byte:
  # if 0 dot is out of range
  # if 1 dot switched off
  # if 2 dot normal
  # if 3 dot is maindot (has been selected as part of call body)

  zc$S <- rep_len(NA_integer_, 16384)

  if (FTYPE == 129)
  {
    # Set pointer to data block
    p <- zc$raw_data[p] + zc$raw_data[p + 1] * 256 + 1L

    # Read data from old 129 file
    zc <- read_zc.129(zc, p, zc$metadata, zc$raw_data)
  }
  else
  {
    if (FTYPE == 132)
    {
      zc$metadata$DATE <- read_date_time.zc(p, zc$raw_data)
      zc$metadata$GPS <- read_gps_data.zc(p, zc$raw_data)
      zc$metadata$ALTITUDE <- read_altitude_data.zc(p, zc$raw_data)
    }
    else
    {
      zc$metadata$DATE <- NULL # Time data stored in the filename
    }

    # Set pointer to data block
    p <- zc$raw_data[p] + zc$raw_data[p + 1] * 256 + 1L

    # Test for the existence of GUANO metadata
    if (p != 0x150 + 1 && FTYPE == 0x84)
    {
      # Parse GUANO data
      raw <- rawToChar(as.raw(zc$raw_data[(0x150+1):(p - 1)]))
      splitted <- strsplit(raw, "[|]")[[1]]

      tmp <- list()

      for (x in splitted[-1])
      {
        tmp <- c(
          tmp,
          if (grepl(x, pattern = "[:]"))
          {
            tmp2 <- strsplit(x, "[\n]")[[1]]
            tmp2 <- tmp2[tmp2 != "Anabat"]

            tmp3 <- lapply(tmp2, function(x) strsplit(x, "[:]")[[1]])
            setNames(unlist(lapply(tmp3, "[[", 2)), unlist(lapply(tmp3, "[[", 1)))
          }
          else
            setNames(list(NULL), x))
      }

      names(tmp)[1] <- paste0(splitted[1L], "|", names(tmp)[1L])

      zc$metadata$GUANO <- as.list(tmp)
    }

    # Read data from file types 130, 131 & 132
    zc <- read_zc.130(zc, p, zc$metadata, FTYPE, zc$raw_data)
  }

  RES1 <- 25000

  zc <- calc_freq.zc(zc, zc$metadata) # Add frequency results to zc
  # zc$S[zc$N] <- -2
  zc$N <- zc$N - 1
  class(zc) <- "zc"

  return(zc)
}


read_zc.129 <- function(zc, p, metadata, raw_data)
{
  if(p != 0x120 + 1L) stop("File is corrupted!")

  RES1 <- metadata$RES1

  t <- 1L
  prevtime_interval <- s <- time <- 0L
  time_data <- integer(16384)

  nBytes <- length(raw_data)

  t <- t + 1L # C/C++ indexing to R indexing

  while (p <= nBytes)
  {
    if (raw_data[p] < 0x80)
    {
      # If bit8 = 0 (same as raw_data < 0x80 or 128):
      #   only one byte is used. The remaining 7 bits form a signed, two
      #   complement number (between -64 and 63) which is added to the value of
      #   the previous time interval to derive the value of the current time interval.

      time_interval <- raw_data[p]
      if (time_interval > 63) time_interval <- time_interval - 128
      prevtime_interval <- prevtime_interval + time_interval

      time_interval <- if (RES1 != 25000) floor(zc$metadata$time_factor * prevtime_interval + .5) else prevtime_interval
      time_data[t] <- time <- time + time_interval
      t <- t + 1
      p <- p + 1
    }
    else
    {
      if (raw_data[p] > 0xf8)
      {
        # If bits 8 to 4 are all 1 (same as raw_data > 0xf8 or 248):
        # 	the first 3 bits tell us how many of the following data points should
        #   be turned off, even if they represent valid signal times. This is
        #   because these points were edited off before the file was saved, but
        #   they are needed to correctly calculate the positions of following points.

        s <- raw_data[p] - 0xf8
        zc$S[t + 1:s] <- 1
        p <- p + 1
      }
      else
      {
        # Bits 4 to 7 form an unsigned value n, and bits 1 to 3 form a value H.
        # The 8 bits of the next data byte form a value L. The resulting time interval
        # is: (H * 256 + L) shifted left n times.

        n <- raw_data[p]
        n <- bitwShiftR(n, 3) # Integer division
        n <- bitwAnd(n, 0x0f)
        time_interval <- bitwAnd(raw_data[p], 7) * 256 + raw_data[p + 1]
        if (n > 0) time_interval <- bitwShiftL(time_interval, n)
        prevtime_interval <- time_interval
        if (RES1 != 25000) time_interval <- floor(zc$metadata$time_factor * prevtime_interval + .5)
        time_data[t] <- time <- time_interval + time
        t <- t + 1
        p <- p + 2
      }
    }
  }

  zc$time_data <- time_data
  zc$N <- t
  zc
}


read_zc.130 <- function(zc, p, metadata, FTYPE, raw_data)
{
  time <- time_interval <- prevtime_interval <- s <- 0L
  n <- 1 # Number of data points
  time_data <- integer(16384)

  n <- n + 1 # C indexing to R indexing

  RES1 <- metadata$RES1
  if (p != 0x120 + 1 && FTYPE < 0x84) stop("File is corrupted!")

  nBytes <- length(raw_data) # File length in bytes

  while (p <= nBytes && n <= 16384)
  {
    if (raw_data[p] < 0x80)
    {
      # If bit8 = 0 (same as raw_data < 0x80 or 128):
      #   only one byte is used. The remaining 7 bits form a signed, two
      #   complement number (between -64 and 63) which is added to the value of
      #   the previous time interval to derive the value of the current time interval.

      time_interval <- raw_data[p]
      if (time_interval > 63) time_interval <- time_interval - 128
      prevtime_interval <- time_interval <- prevtime_interval + time_interval
      if (RES1 != 25000) time_interval <- floor(metadata$time_factor * prevtime_interval + 0.5)
      time_data[n] <- time <- time + time_interval
      n <- n + 1
      p <- p + 1
    }
    else
    {
      if (raw_data[p] >= 0xe0) # Indicates status change
      {
        # If bits 8 to 4 are all 1 (same as raw_data > 0xf8 or 248):
        # 	the first 3 bits tell us how many of the following data points should
        #   be turned off, even if they represent valid signal times. This is
        #   because these points were edited off before the file was saved, but
        #   they are needed to correctly calculate the positions of following points.

        if (FTYPE > 130) # File types 131 and 132
        {
          if (p + 1 >= nBytes) break
          c <- bitwAnd(raw_data[p], 3)
          s <- raw_data[p + 1]
          if ((n + s - 1) > 16383) s <- 16384 - n # limit index to arrays
          zc$S[n + 1:s] <- c
          p <- p + 2

        }
        else # File type 130
        {
          s <- raw_data[p] - 0xe0
          if ((n + s - 1) > 16383) s <- 16384 - n # limit index to arrays
          zc$S[n + 1:s - 1] <- 1
          p <- p + 1
        }
      }
      else
      {
        # Bits 6 and 7 indicate the case:
        #   0: time interval encoded in one byte
        #  32: time interval encoded in two bytes
        #  64: time interval encoded in three bytes

        case <- bitwAnd(raw_data[p], 0x60)

        if (case == 0)
        {
          zc$S[n] <- 2

          if (p + 1 >= nBytes) break
          prevtime_interval <- time_interval <- bitwShiftL(bitwAnd(raw_data[p], 0x1f), 8) + raw_data[p + 1]
          if (RES1 != 25000) time_interval <- floor(metadata$time_factor * prevtime_interval + 0.5)
          time_data[n] <- time <- time + time_interval
          n <- n + 1
          p <- p + 2
        }
        else if (case == 0x20) # 32
        {
          if (p + 2 >= nBytes) break
          prevtime_interval <- time_interval <- bitwShiftL(bitwAnd(raw_data[p], 0x1f), 16) + bitwShiftL(raw_data[p + 1], 8) + raw_data[p + 2]
          if (RES1 != 25000) time_interval <- floor(metadata$time_factor * prevtime_interval + 0.5)
          time_data[n] <- time <- time + time_interval
          n <- n + 1
          p <- p + 3
        }
        else if (case == 0x40) # 64
        {
          if (p + 3 >= nBytes) break
          prevtime_interval <- time_interval <- bitwShiftL(bitwAnd(raw_data[p], 0x1f), 24) + bitwShiftL(raw_data[p + 1], 16) + bitwShiftL(raw_data[p + 2], 8) + raw_data[p + 3]
          if (RES1 != 25000) time_interval <- floor(metadata$time_factor * prevtime_interval + 0.5)
          time_data[n] <- time <- time + time_interval
          n <- n + 1
          p <- p + 4
        }
      }
    }
  }

  zc$time_data <- time_data
  zc$N <- n
  class(zc) <- "zc"
  zc
}


read_altitude_data.zc <- function(p, raw_data)
{
  as.integer(rawToChar(as.raw(raw_data[p + 50:53])))
}

read_date_time.zc <- function(p, raw_data)
{
  # Date and time at start of file
  # 0x120 word year (eg 2001)
  # 0x122 byte month (Jan=1, Dec=12)
  # 0x123 byte day (1st=1 etc)
  # 0x124 byte hour (0-23)
  # 0x125 byte minute (0-59)
  # 0x126 byte second (0-59)
  # 0x127 byte hundredsth of a second (0-99)
  # 0x128 word microseconds (0-9999)

  YEAR <- raw_data[p + 6] + 256 * raw_data[p + 7]
  MON <- stringr::str_pad(raw_data[p + 8], width = 2, pad = 0)
  DAY <- stringr::str_pad(raw_data[p + 9], width = 2, pad = 0)
  HOURS <- stringr::str_pad(raw_data[p + 10], width = 2, pad = 0)
  MINS <- stringr::str_pad(raw_data[p + 11], width = 2, pad = 0)
  SECS <- stringr::str_pad(raw_data[p + 12], width = 2, pad = 0)
  HUNDS <- raw_data[p + 13]
  MICROS <- raw_data[p + 14] + 256 * raw_data[p + 15]
  paste0(YEAR, MON, DAY, "_", HOURS, ":", MINS, ":", SECS, ".", str_pad(HUNDS * 100 + MICROS, width = 6, pad = 0))
}

read_gps_data.zc <- function(p, raw_data)
{
  # 0x0130-0x014f 32 bytes GPS position data
  # 0x0130-0x0139 10 bytes specifying datum
  # 0x013a if numeric, specifies UTM data
  # 0x013a-0x013c 3 bytes for Zone (eg "10S")
  # 0x013d char(32)
  # 0x013e-0x0143 6 bytes for Easting in metres
  # 0x0144 char(32)
  # 0x0145-0x014b 7 bytes for Northing in metres
  # 0x013a if 'N' or 'S', specifies Degree data
  # 0x013b-0x013c 2 bytes for latitude Degrees
  # 0x013d-0x0141 5 bytes for decimal Degrees
  # 0x0142 char(32)
  # 0x0143 'W' or 'E'
  # 0x0144-0x0146 3 bytes for longitude Degrees
  # 0x0147-0x014b 5 bytes for decimal Degrees
  # 0x014c-0x014f 4 bytes for Altitude in metres (-999 to 9999)

  DATUM <- stringr::str_trim(rawToChar(as.raw(raw_data[p + 22:31])))
  UTM_OR_DEGREE <- rawToChar(as.raw(raw_data[p + 32]))

  if (UTM_OR_DEGREE %in% c('N', 'S'))
  {
    LAT_DEG <- as.integer(rawToChar(as.raw(raw_data[p + 33:34])))
    LAT_DEC_DEG <- rawToChar(as.raw(raw_data[p + 35:39]))
    LAT <- paste0(as.integer(LAT_DEG), ".", LAT_DEC_DEG, " ", UTM_OR_DEGREE)
    W_OR_E <- rawToChar(as.raw(raw_data[p + 41]))
    LONG_DEG <- as.integer(rawToChar(as.raw(raw_data[p + 42:44])))
    LONG_DEC_DEG <- rawToChar(as.raw(raw_data[p + 45:49]))
    LONG <- paste0(LONG_DEG, ".", LONG_DEC_DEG, " ", W_OR_E)
    GPS <- paste0(DATUM, " ", LAT, ", ", LONG)

  }
  else
  {
    ZONE <- rawToChar(as.raw(raw_data[p + 32:34]))
    EASTING <- rawToChar(as.raw(raw_data[p + 36:41]))
    NORTHING <- rawToChar(as.raw(raw_data[p + 43:49]))
    GPS <- paste0(DATUM, " ", ZONE, " ", EASTING, ", ", NORTHING)
  }

  GPS
}


read_header.zc <- function(raw_data)
{
  raw_data <- as.raw(raw_data[7:281])

  # At some point in the past it appears that there was a bug in the Anabat
  # software so that null or zero bytes would be written into the text header
  # fields instead of ASCII 32 for a space. It only shows up in a few files for
  # example in the NSW Bat Call Library data set. This patch should fix this
  # problem. Patch contributed by Peter Wilson (peter@peterwilson.id.au)
  raw_data[raw_data == 0] <- charToRaw(" ")

  TAPE <- stringr::str_trim(rawToChar(raw_data[1:8]))
  DATE <- stringr::str_trim(rawToChar(raw_data[9:16]))
  LOC <- stringr::str_trim(rawToChar(raw_data[17:56]))
  SPECIES <- stringr::str_trim(rawToChar(raw_data[57:106]))
  SPEC <- stringr::str_trim(rawToChar(raw_data[107:122]))
  NOTE <- stringr::str_trim(rawToChar(raw_data[123:196]))
  NOTE1 <- stringr::str_trim(rawToChar(raw_data[197:275]))

  list(TAPE = TAPE, DATE = DATE, LOC = LOC, SPECIES = SPECIES, SPEC = SPEC, NOTE = NOTE, NOTE1 = NOTE1)
}

read_settings.zc <- function(p, raw_data)
{
  RES1 <- raw_data[p + 2] + 256 * raw_data[p + 3] # Number of counts (derived from ZCAIM) which represent a time interval
                                                  # of 25ms = 25000 unless recalibration is in effect
  time_factor <- RES1 / 25000
  DIVRAT <- raw_data[p + 4] # Division ratio

  VRES <- raw_data[p + 5] # This value is used to determine what value of SCALE is to be used.
                          # The value (VRES AND 0x70) / 16 is used to lookup a table which gives SCALE, ie:
                          #
                          # (VRES AND 70h)/16    SCALE
                          #         0               10
                          #         1               25
                          #         2               50
                          #         3              100
                          #         4              250
                          #         5              500
                          #         6             1000
                          #         7             2500

  list(RES1 = RES1, DIVRAT = DIVRAT, VRES = VRES, time_factor = time_factor)
}

# Calculates the frequency data from the time data
# filling F with data derived from DIVRAT and time_data

calc_freq.zc <- function(zc, metadata)
{
  DIVRAT <- metadata$DIVRAT # Division ratio

  time_data <- zc$time_data
  N <- zc$N

  freq_data <- rep_len(NA, length(time_data)) # Frequency data (Hz)

  zc$S[1:2] <- if (zc$FTYPE == 129) 0 else 0:1

  Tmin <- ceiling(DIVRAT * 4)  # corresponds to 4 kHz
  Tmax <- floor(DIVRAT * 250)  # corresponds to 250 kHz
  if (Tmin < 48) Tmin <- 48
  if (Tmax > 12859) Tmax <- 12859

  time_diff <- diff(time_data[1:N], lag = 2)
  freq_data[3:N][time_diff >= Tmin & time_diff <= Tmax] <- (DIVRAT * 1e6) %/% time_diff[time_diff >= Tmin & time_diff <= Tmax]
  # zc$S[3:N][td >= Tmin & td <= Tmax] <- 2

  zc$freq_data <- freq_data
  zc
}


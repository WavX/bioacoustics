#-------------------------------------------------------------------------------
# Copyright (C) 2017 David A. Riggs (driggs@myotisoft.com)
# Copyright (C) 2019 WavX, inc. (www.wavx.ca)
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
# Originally written by David A. Riggs (driggs@myotisoft.com)
# Available here: https://github.com/riggsd/guano-r
# Further improvements and additions by WavX inc. (info@wavx.ca)
#-------------------------------------------------------------------------------

#' Read GUANO metadata in audio file
#'
#' @param file Path to a wav file
#' @return list of named metadata fields
#'
guano_md <- function(file)
{
  file_checks(file)

  f <- file(file, "rb")
  on.exit(close(f))

  RIFF <- readChar(f, 4)

  if (length(RIFF) == 0 || RIFF != "RIFF")
  {
    stop("This does not seem to be a valid wav file.")
  }

  riff.size <- readBin(f, integer(), size = 4, endian = "little")
  WAVE <- readChar(f, 4)  # "WAVE"

  if (length(WAVE) == 0 || WAVE != "WAVE")
  {
    stop("This does not seem to be a valid wav file.")
  }

  read.subchunk <- function()
  {
    id <- tryCatch(
      readChar(f, 4),
      error = function(e)
      {
        return("")
      }
    )

    if (length(id) == 0 || id == "")
      return(NULL)

    size <- readBin(f, integer(), size = 4, endian = "little")
    list(id = id, size = size)
  }

  skip.subchunk <- function(chunk)
  {
    #print(sprintf("Skipping subchunk '%s' ...", chunk$id))
    pos <- seek(f, NA)
    seek(f, pos + chunk$size)
  }

  # Maps metadata keys to a data type coercion function
  data_types <- list(
    `Filter HP` = as.double,
    `Filter LP` = as.double,
    Humidity = as.double,
    Length = as.double,
    `Loc Accuracy` = as.integer,
    `Loc Elevation` = as.double,
    Note = function(val) gsub("\\\\n", "\n", val),
    Samplerate = as.integer,
    #`Species Auto ID`=?, `Species Manual ID`=?,  # TODO: comma separated
    #Tags=?,  # TODO: comma separated
    TE = function(val) if (is.na(val) || is.null(val) || val == "") 1 else as.integer(val),
    `Temperature Ext` = as.double,
    `Temperature Int` = as.double,
    Timestamp = .parse.timestamp
  )


  md <- list()

  while (!is.null(chunk <- read.subchunk()))
  {
    if (chunk$id != "guan")
    {
      skip.subchunk(chunk)
      next
    }

    md.txt <- readChar(f, chunk$size)
    #Encoding(md.txt) <- "UTF-8"  # FIXME: this still isn't setting the encoding to UTF-8

    for (line in strsplit(md.txt, "\n")[[1]])
    {
      line <- trimws(line)

      if (line == "")
        next

      toks <- strsplit(sub(":", "\n", line), "\n")
      key <- trimws(toks[[1]][1])
      val <- trimws(toks[[1]][2])

      if (is.na(key) || is.null(key) || key == "")
        next

      if (!is.null(data_types[[key]]))
        val <- data_types[[key]](val)

      md[[key]] <- val
    }

    if ("Loc Position" %in% names(md))
    {
      coords <- lapply(strsplit(md[["Loc Position"]], " "), as.double)[[1]]
      md[["Loc Position Lat"]] <- coords[1]
      md[["Loc Position Lon"]] <- coords[2]
      md[["Loc Position"]] <- NULL
    }
  }

  return(md)
}

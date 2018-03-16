//------------------------------------------------------------------------------
//  Copyright (C) 2012 Chris Scott (fbscds@gmail.com)
//  Copyright (C) 2017-2018 WavX, inc. (www.wavx.ca)
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with This program.  If not, see <https://www.gnu.org/licenses/>.
//------------------------------------------------------------------------------

#ifndef BB_EXTRACT_H
#define BB_EXTRACT_H

#include <Rcpp.h>
#include "fft.h"
#include "bb_audio_event.h"

void extract_impl (Audio_Event &audio_event,
                   const size_t &sample_rate,
                   const size_t &fft_size,
                   const size_t &step,
                   Rcpp::List &out,
                   const size_t &index,
                   const double KPE,
                   const double KME);

#endif // BB_EXTRACT_H

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
//  along with This program. If not, see <https://www.gnu.org/licenses/>.
//------------------------------------------------------------------------------

#ifndef BB_AUDIO_EVENT_H
#define BB_AUDIO_EVENT_H

#include <Rcpp.h>
#include <vector>

class Audio_Event
{
public:
  size_t duration = 0;
  int end, start;
  double noise = 0, amp_peak = 0, signal = 0;
  Rcpp::NumericVector amp_track, freq_track, harmonic_amp_track;
  std::vector<double> power_spectrum;
};

#endif // BB_AUDIO_EVENT_H

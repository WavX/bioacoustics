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

#ifndef BB_DETECT_H
#define BB_DETECT_H

#include "fft.h"

// std::vector< std::vector< std::vector<double> > > detect_impl (const std::vector<int> &audio_samples,
void detect_impl (const std::vector<int> &audio_samples,
                  std::deque<int> &peak_locations,
                  std::deque< std::vector<double> > &background_noises,
                  const size_t &sample_rate,
                  const size_t &threshold,
                  FFT &fft,
                  const size_t &LPF,
                  const size_t &HPF,
                  const double &dur_t,
                  const double &EDG,
                  const size_t &noise_estim_win_size);

#endif // BB_DETECT_H

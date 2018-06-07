//------------------------------------------------------------------------------
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

#ifndef SPECTRO_H
#define SPECTRO_H

#include <string>
#include <vector>
#include <Rcpp.h>
#include "fft.h"

// std::vector< std::vector<double> > spectro(const std::vector<int> &audio_samples,
//                                            size_t fft_size,
//                                            size_t fft_step,
//                                            std::string win,
//                                            const size_t& HPF,
//                                            const size_t& LPF);
Rcpp::NumericMatrix fspec_impl(const std::vector<int>& audio_samples,
                               const size_t& fft_size,
                               const double& fft_overlap,
                               const std::string& win,
                               const size_t& HPF,
                               const size_t& LPF,
                               const size_t& FLL_bin,
                               const size_t& FUL_bin,
                               const bool& rotate);

#endif // SPECTRO_H

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
//  along with This program.  If not, see <https://www.gnu.org/licenses/>.
//------------------------------------------------------------------------------

#include <cmath>
#include <string>
#include <vector>
#include <Rcpp.h>
#include "fft.h"

// [[Rcpp::export(.fspec_impl)]]
Rcpp::NumericMatrix fspec_impl(const std::vector<int>& audio_samples,
                               const size_t& fft_size,
                               const double& fft_overlap,
                               const std::string& win,
                               const size_t& HPF_bin,
                               const size_t& LPF_bin,
                               const size_t& FLL_bin,
                               const size_t& FUL_bin,
                               const bool& rotate)
{
  FFT::WIN_TYPE win_type = fft_win_str_to_enum(win);
  FFT fft(fft_size, win_type);
  size_t seek_step = std::max((double)fft_size * (1 - fft_overlap), (double)1);

  int n_rows = FUL_bin - FLL_bin + 1;
  int k = FLL_bin - 1;
  int n_cols = std::ceil((double)audio_samples.size() / (double)seek_step);

  if (rotate)
  {
    Rcpp::NumericMatrix mat(n_cols, n_rows);

    for (size_t seek = 0, j = 0; seek < audio_samples.size(); seek+=seek_step, j++)
    {
      fft.impl(seek, audio_samples);

      for (size_t i = LPF_bin + 1; i-- > HPF_bin;)
      {
        double x = fft.magnitude[i];
        mat(j, i - HPF_bin) = x;
      }
    }

    return mat;
  }
  else
  {
    Rcpp::NumericMatrix mat(n_rows, n_cols);

    for (size_t seek = 0, j = 0; seek < audio_samples.size(); seek+=seek_step, j++)
    {
      fft.impl(seek, audio_samples);

      for (size_t i = LPF_bin + 1; i-- > HPF_bin;)
      {
        double x = fft.magnitude[i];
        mat(n_rows - i + k, j) = x;
      }
    }

    return mat;
  }
}



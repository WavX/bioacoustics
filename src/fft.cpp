//------------------------------------------------------------------------------
//  Copyright (C) 2017-2019 WavX, inc. (www.wavx.ca)
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

#include <algorithm>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>
#include "fft.h"

FFT::FFT() {}

FFT::FFT(size_t fft_sz, FFT::WIN_TYPE win_type)
{
  this->fft_size = fft_sz;

  set_window(win_type);
  set_plan(fft_sz);
}

void FFT::set_plan(const size_t &fft_sz)
{
  original.resize(fft_sz, 0);
  transformed.resize(fft_sz, 0);
  magnitude.resize(fft_sz / 2, 0);
  plan = fftw_plan_r2r_1d(fft_sz, &original[0], &transformed[0], FFTW_R2HC, FFTW_PATIENT);
}

void FFT::set_window(const WIN_TYPE& win_type)
{
  window.resize(fft_size, 0);
  z = M_PI / (fft_size-1);

  switch(win_type)
  {
  case WIN_TYPE::BLACKMAN_HARRIS_4 :
    blackman_harris_4 (fft_size);
    break;
  case WIN_TYPE::BLACKMAN_HARRIS_7 :
    blackman_harris_7 (fft_size);
    break;
  case WIN_TYPE::HANN :
    hann (fft_size);
    break;
  }

  double sum { std::accumulate(window.begin(), window.end(), 0.0) };
  normalise = 1. / sum;
}

// FFT windows
FFT::WIN_TYPE fft_win_str_to_enum(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  FFT::WIN_TYPE win_type;

  if (s == "hann")
  {
    win_type = FFT::WIN_TYPE::HANN;
  }
  else if (s == "blackman4")
  {
    win_type = FFT::WIN_TYPE::BLACKMAN_HARRIS_4;
  }
  else if (s == "blackman7")
  {
    win_type = FFT::WIN_TYPE::BLACKMAN_HARRIS_7;
  }
  else Rcpp::stop("This type of window is not implemented.");

  return win_type;
}

void FFT::blackman_harris_4 (size_t fft_sz)
{
  for (size_t i = 0; i < fft_sz; i++)
  {
    window[i] = 0.35875 - 0.48829 * std::cos(2*z*i) + \
      0.14128 * std::cos(4*z*i) -                   \
      0.01168 * std::cos(6*z*i);
  }
}


void FFT::blackman_harris_7 (size_t fft_sz)
{
  // doi:10.1109/icassp.2001.940309
  for (size_t i = 0; i < fft_sz; i++)
  {
    window[i] = 0.2712203606 - 0.4334446123 * std::cos(2*z*i) +      \
      0.21800412 * std::cos(4*z*i) - 0.0657853433 * std::cos(6*z*i) +   \
      0.0107618673 * std::cos(8*z*i) - 0.0007700127 * std::cos(10*z*i) + \
      0.00001368088 * std::cos(12*z*i);
  }
}

void FFT::hann (size_t fft_sz)
{
  for (size_t i = 0; i < fft_sz; i++)
  {
    window[i] = 0.5 * (1 - std::cos(2*z*i));
  }
}


FFT::~FFT() { fftw_destroy_plan( plan ); }

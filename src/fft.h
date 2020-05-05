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

#ifndef FFT_H
#define FFT_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <stddef.h>
#include <vector>
#include <fftw3.h>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

class FFT
{
public:
  enum class WIN_TYPE {
    BLACKMAN_HARRIS_4 = 0,
    BLACKMAN_HARRIS_7 = 1,
    HANN = 2
  };

  FFT();
  FFT(size_t fft_sz, WIN_TYPE win_type);
  ~FFT();

  std::vector<double> magnitude, original, transformed;

  inline void impl(std::size_t seek, const std::vector<int> &samples);
  void set_plan(const size_t &fft_sz);
  void set_window(const FFT::WIN_TYPE& win_type);

  std::size_t fft_size;

private:
  double normalise, z;
  std::vector<double> window;
  fftw_plan plan;

  void blackman_harris_4 (const size_t fft_sz);
  void blackman_harris_7 (const size_t fft_sz);
  void hann (const size_t fft_sz);
};


inline
  void FFT::impl(size_t seek, const std::vector<int> &samples)
  {
    size_t N = samples.size();

    std::fill(original.begin(), original.end(), 0.0);
    std::fill(transformed.begin(), transformed.end(), 0.0);

    for (size_t i = 0; i < fft_size; i++, seek++)
    {
      if (seek < N)
      {
        original[i] = samples[seek] * window[i];
      }
    }

    fftw_execute(plan);
    size_t sk = fft_size;

    for (size_t i = 0; i < fft_size / 2; i++)
    {
      magnitude[i] = std::abs(std::complex<double>(transformed[i], transformed[--sk])) * normalise;
    }
  }

FFT::WIN_TYPE fft_win_str_to_enum(std::string s);

#endif // FFT_H

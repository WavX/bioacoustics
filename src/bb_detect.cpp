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

#include <cfloat>
#include <cmath>
#include <deque>
#include <vector>
#include "fft.h"
#include "bb_tools.h"
#include <Rcpp.h>

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
                  const size_t &noise_estim_win_size)
{

  double freq_res = (double)sample_rate / fft.fft_size;
  int n_samples = audio_samples.size();
  size_t step_size = fft.fft_size / 2;
  size_t frames = 1 + (n_samples - fft.fft_size) / step_size;
  int seek = -step_size;

  double EXPT = -DBL_MAX;
  double sum_squares = 0;

  while(sum_squares < 0.00001 && seek < n_samples)
  {
    for (int j = 0; j < (int)fft.fft_size; j++)
    {
      if ((seek + j) >= 0 && (seek + j) < n_samples)
      {
        sum_squares += std::pow(audio_samples.at(seek + j), 2);
      }
    }
    seek += step_size;
  }

  fft.impl(seek, audio_samples);
  std::vector<double> power_spectrum = fft.magnitude;

  band_pass_filter(power_spectrum, LPF, HPF, freq_res);

  std::vector<double> prev = power_spectrum, background_noise_filtered_prev = power_spectrum, background_noise_filtered_new = power_spectrum;

  size_t win_size = ((float)noise_estim_win_size * (float)sample_rate / 1000) / (float)step_size;
  std::vector< std::vector<double> > bg_noiz_f_win (fft.fft_size / 2, std::vector<double> (win_size));

  size_t bg_noiz_f_insert_at = 0, bg_noiz_f_win_size = 0;

  bool triggered = false;
  double const step_ms = (double)1000 * (double)step_size / (double)sample_rate;
  double duration = 0;
  double local_SNR = ZEROLOG, max_SNR = ZEROLOG;
  size_t peak_location = 0;

  for (size_t i = 1; i < frames; i++)
  {
    if (seek > n_samples) continue;

    sum_squares = 0;
    for (size_t j = 0; j < fft.fft_size; j++)
    {
      if ((seek + (int)j) >= 0 && (seek + (int)j) < n_samples)
      {
        sum_squares += std::pow(audio_samples.at(seek + j), 2);
      }
    }

    if (sum_squares < 0.00001)
    {
      seek += step_size;
      continue;
    }

    fft.impl(seek, audio_samples);
    power_spectrum = fft.magnitude;
    seek += step_size;

    std::vector<double> filtered_spectrum = power_spectrum;
    band_pass_filter(filtered_spectrum, LPF, HPF, freq_res);

    for (size_t i = 0; i < prev.size(); i++)
    {
      double y = filtered_spectrum[i];
      filtered_spectrum[i] = (filtered_spectrum[i] + prev[i]) / (double)2;
      prev[i] = y;
    }

    if (!triggered || duration > dur_t)
    {
      // wrap around if end of buffer is reached
      if (bg_noiz_f_insert_at == win_size) bg_noiz_f_insert_at = 0;

      if (bg_noiz_f_win_size < win_size)
      {
        bg_noiz_f_win_size++;
        for (size_t i = 0; i < bg_noiz_f_win.size(); i++)
        {
          background_noise_filtered_new[i] =   \
            background_noise_filtered_new[i] + \
            (filtered_spectrum[i] - background_noise_filtered_new[i]) / (double)bg_noiz_f_win_size;
        }
      }
      else
      {
        for (size_t i = 0; i < bg_noiz_f_win.size(); i++)
        {
          background_noise_filtered_new[i] =                             \
            background_noise_filtered_new[i] -                           \
            bg_noiz_f_win[i][bg_noiz_f_insert_at] / bg_noiz_f_win_size + \
            filtered_spectrum[i] / bg_noiz_f_win_size;
        }
      }

      for (size_t i = 0; i < bg_noiz_f_win.size(); i++)
      {
        bg_noiz_f_win[i][bg_noiz_f_insert_at] = filtered_spectrum[i];
      }
      bg_noiz_f_insert_at++;
    }

    local_SNR = SNR(filtered_spectrum, background_noise_filtered_prev);

    EXPT = std::max(local_SNR, EDG * EXPT + (1 - EDG) * local_SNR);

    if (local_SNR >= threshold && local_SNR >= EXPT && !triggered)
    {
      triggered = true;
    }

    if (triggered)
    {
      duration += step_ms;
      if (local_SNR > max_SNR)
      {
        peak_location = std::max(seek - 2 * (int)step_size, 0);
        max_SNR = local_SNR;
      }

      if (local_SNR < threshold)
      {

        peak_locations.push_back(peak_location);
        background_noises.push_back(background_noise_filtered_prev);

        // reset variables
        duration = 0;
        triggered = false;
        max_SNR = ZEROLOG;
      }
    }
    background_noise_filtered_prev = background_noise_filtered_new;
  }
}

//------------------------------------------------------------------------------
//  Copyright (C) 2012 Chris Scott (fbscds@gmail.com)
//  Copyright (C) 2017-2018 WavX inc. (www.wavx.ca)
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

#ifndef BB_ANALYSE_H
#define BB_ANALYSE_H

#include "fft.h"
#include "bb_audio_event.h"
#include "bb_kalman.h"

class Analyse
{
public:
  Analyse(const std::vector<int> &audio_samples,
          const size_t &LPF,
          const size_t &HPF,
          const size_t &sample_rate,
          const int &step,
          const double &start_t,
          const double &end_t,
          const double &angl_t,
          const double &snr_t);

  Audio_Event impl(const size_t &fft_size,
                   const int &peak_location,
                   std::vector<double> &background_noise,
                   const double KPE,
                   const double KME);

  FFT fft;

private:
  Kalman kalman;

  const std::vector<int> &audio_samples;
  const size_t &LPF, &HPF, &sample_rate;
  const int &step;
  const double &start_t, &end_t, &angl_t, &snr_t;

  size_t bin_fundamental, bin_harmonic;
  double energy = 0, freq_res, smoothing_gain = 0.25;
  double bin_centroid;
  std::vector<double> power_spectrum, mask;

  void analyse_frame (const int &seek, double &noise, double &signal, const std::vector<double> &background_noise);
  void backward_analyse (Audio_Event &audio_event, int &seek, const size_t &peak_location, const std::vector<double> &background_noise, double &noise, double &signal);
  void forward_analyse (Audio_Event &audio_event, int &seek, const std::vector<double> &background_noise, double &noise, double &signal);
  bool is_start(Audio_Event &audio_event, const double &angl_diff, const double &amp_diff, const double &SNR, const int &seek);
  bool is_end(Audio_Event &audio_event, const double &angl_diff, const double &amp_diff, const double &SNR, const int &seek);
  void store_back (Audio_Event &audio_event, const double &noise, const double &signal);
  void store_front (Audio_Event &audio_event, const double &noise, const double &signal);

};

#endif // BB_ANALYSE_H

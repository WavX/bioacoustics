//------------------------------------------------------------------------------
//  Copyright (C) 2012 Chris Scott (fbscds@gmail.com)
//  Copyright (C) 2017 WavX, inc. (www.wavx.ca)
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
#include <deque>
#include <vector>
#include "fft.h"
#include "bb_analyse.h"
#include "bb_audio_event.h"
#include "bb_tools.h"
#include "bb_kalman.h"

Analyse::Analyse(const std::vector<int> &audio_samples,
                 const size_t &LPF,
                 const size_t &HPF,
                 const size_t &sample_rate,
                 const int &step,
                 const double &start_t,
                 const double &end_t,
                 const double &angl_t,
                 const double &snr_t) :
  audio_samples(audio_samples),
  LPF(LPF),
  HPF(HPF),
  sample_rate(sample_rate),
  step(step),
  start_t(start_t),
  end_t(end_t),
  angl_t(angl_t),
  snr_t(snr_t)
  {}

Audio_Event Analyse::impl(const size_t &fft_size,
                          const int &peak_location,
                          std::vector<double> &background_noise,
                          const double KPE,
                          const double KME)
{
  double noise = 0, signal = 0;
  int seek = peak_location;
  fft.fft_size = fft_size;
  fft.set_plan(fft_size);
  fft.set_window(FFT::WIN_TYPE::BLACKMAN_HARRIS_7);
  freq_res = sample_rate / fft.fft_size;

  smooth_spectrum(background_noise, smoothing_gain);
  analyse_frame(seek, noise, signal, background_noise);
  seek += step;
  kalman = Kalman(KPE, KME, bin_centroid);

  Audio_Event audio_event;
  audio_event.power_spectrum.resize(fft.fft_size / 2, 0);
  store_back(audio_event, noise, signal);

  forward_analyse(audio_event, seek, background_noise, noise, signal);
  backward_analyse(audio_event, seek, peak_location, background_noise, noise, signal);

  return audio_event;
}

void Analyse::analyse_frame (const int &seek,
                             double &noise,
                             double &signal,
                             const std::vector<double> &background_noise)
{
  noise = 0;
  signal = 0;
  fft.impl(seek, audio_samples);
  power_spectrum = fft.magnitude;
  band_pass_filter(power_spectrum, LPF, HPF, freq_res);

  smooth_spectrum(power_spectrum, smoothing_gain);

  for (size_t i = 0; i < background_noise.size(); i++)
  {
    if (power_spectrum[i] > background_noise[i])
    {
      signal += power_spectrum[i];
      noise += background_noise[i];
    }
  }

  for (size_t i = 0; i < power_spectrum.size(); i++)
  {
    power_spectrum[i] -= background_noise[i];
    power_spectrum[i] = (power_spectrum[i] + std::abs(power_spectrum[i])) / (double)2;
  }

  size_t peak = std::distance(power_spectrum.begin(), std::max_element(power_spectrum.begin(), power_spectrum.end()));
  mask = gaussian_mask(peak, power_spectrum);
  bin_fundamental = peak;
  energy = std::accumulate(mask.begin(), mask.end(), (double)0);

  // tracked frequency - weighted average
  bin_centroid = centroid(mask);

  // Harmonic octave above tracked
  bin_harmonic = std::round(2 * bin_centroid);
}


void Analyse::backward_analyse (Audio_Event &audio_event,
                                int &seek,
                                const size_t &peak_location,
                                const std::vector<double> &background_noise,
                                double &noise,
                                double &signal)
{
  seek = peak_location - step;

  kalman.p_error_prev = 1;
  kalman.p_state_prev = bin_centroid;
  kalman.data.clear();

  size_t win_size = 5;
  std::vector<double> snr_win (win_size);

  size_t snr_win_insert_at = 0, snr_win_size = 0;

  while (seek >= 0)
  {
    analyse_frame(seek, noise, signal, background_noise);
    seek -= step;
    kalman.impl(bin_centroid);
    kalman.data.push_back(kalman.p_state_prev);

    double AD = ang_diff(kalman.data);
    double amp_diff = to_dB(audio_event.amp_peak) - to_dB(energy);
    double SNR = to_dB(signal / no_zero(noise));

    if (snr_win_insert_at == win_size) snr_win_insert_at = 0;

    snr_win[snr_win_insert_at] = SNR;

    snr_win_insert_at++;
    if (snr_win_size < win_size) snr_win_size++;

    double sum = 0;
    for (size_t i = 0; i < snr_win_size; i++)
    {
      sum += snr_win[i];
    }
    SNR = sum / (double)snr_win_size;

    if(is_start(audio_event, AD, amp_diff, SNR, seek)) break;
    store_front(audio_event, noise, signal);
  }

}


void Analyse::forward_analyse (Audio_Event &audio_event,
                               int &seek,
                               const std::vector<double> &background_noise,
                               double &noise,
                               double &signal)
{
  size_t frames = (audio_samples.size() - seek) / step;

  size_t win_size = 5;
  std::vector<double> snr_win (win_size);
  size_t snr_win_insert_at = 0, snr_win_size = 0;

  for (size_t i = 0; i < frames; i++)
  {
    analyse_frame(seek, noise, signal, background_noise);
    seek += step;
    kalman.impl(bin_centroid);
    kalman.data.push_back(kalman.p_state_prev);
    double AD = ang_diff(kalman.data);
    double amp_diff = to_dB(audio_event.amp_peak) - to_dB(energy);
    double SNR = to_dB(signal / no_zero(noise));

    if (snr_win_insert_at == win_size) snr_win_insert_at = 0;

    snr_win[snr_win_insert_at] = SNR;

    snr_win_insert_at++;
    if (snr_win_size < win_size) snr_win_size++;

    double sum = 0;
    for (size_t i = 0; i < snr_win_size; i++)
    {
      sum += snr_win[i];
    }
    SNR = sum / (double)snr_win_size;

    if(i > 1 && is_end(audio_event, AD, amp_diff, SNR, seek)) break;
    store_back(audio_event, noise, signal);
  }
}


bool Analyse::is_end (Audio_Event &audio_event,
                      const double &angl_diff,
                      const double &amp_diff,
                      const double &SNR,
                      const int &seek)
{
  if (amp_diff > end_t || angl_diff > angl_t || SNR < snr_t)
  {
    audio_event.end = seek - step - 1;
    return true;
  }
  return false;
}


bool Analyse::is_start (Audio_Event &audio_event,
                        const double &angl_diff,
                        const double &amp_diff,
                        const double &SNR,
                        const int &seek)
{
  if (amp_diff > start_t || angl_diff > angl_t || SNR < snr_t)
  {
    audio_event.start = seek + 2 * step;
    return true;
  }
  return false;
}


void Analyse::store_back (Audio_Event &audio_event,
                          const double &noise,
                          const double &signal)
{
  audio_event.signal += signal;
  audio_event.noise += noise;
  audio_event.amp_track.push_back(power_spectrum[bin_fundamental]);
  audio_event.freq_track.push_back(bin_centroid);
  audio_event.amp_peak = std::max(audio_event.amp_peak, energy);
  audio_event.duration++;
  std::transform(mask.begin(), mask.end(), audio_event.power_spectrum.begin(), audio_event.power_spectrum.begin(), std::plus<double>());
  audio_event.harmonic_amp_track.push_back(power_spectrum[std::min(bin_harmonic, fft.fft_size/2 - 1)]);
}


void Analyse::store_front (Audio_Event &audio_event,
                  const double &noise,
                  const double &signal)
{
  audio_event.signal += signal;
  audio_event.noise += noise;
  audio_event.amp_track.push_front(power_spectrum[bin_fundamental]);
  audio_event.freq_track.push_front(bin_centroid);
  audio_event.amp_peak = std::max(audio_event.amp_peak, energy);
  audio_event.duration++;
  std::transform(mask.begin(), mask.end(), audio_event.power_spectrum.begin(), audio_event.power_spectrum.begin(), std::plus<double>());
  audio_event.harmonic_amp_track.push_front(power_spectrum[std::min(bin_harmonic, fft.fft_size/2 - 1)]);
}



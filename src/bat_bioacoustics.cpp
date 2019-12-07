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

#include <algorithm>
#include <deque>
#include <vector>
#include <Rcpp.h>
#include "fft.h"
#include "bb_audio_event.h"
#include "bb_detect.h"
#include "bb_analyse.h"
#include "bb_extract.h"

// [[Rcpp::export]]
Rcpp::List threshold_detection_impl(
    const std::vector<int>& audio_samples,
    size_t sample_rate,
    size_t threshold,
    double min_d,
    double max_d,
    double min_TBE,
    double max_TBE,
    double EDG,
    size_t LPF,
    size_t HPF,
    double dur_t,
    double snr_t,
    double angl_t,
    size_t FFT_size,
    double FFT_overlap,
    double start_t,
    double end_t,
    const size_t NWS,
    double KPE,
    double KME
)
{
  std::deque<int> peak_locations;
  std::deque< std::vector<double> > background_noises;

  FFT fft(FFT_size, FFT::WIN_TYPE::BLACKMAN_HARRIS_7);
  fft.set_plan(FFT_size);
  detect_impl(audio_samples, peak_locations, background_noises, sample_rate, threshold, fft, LPF, HPF, dur_t, EDG, NWS);


  if (peak_locations.empty()) return(Rcpp::List());

  size_t duration = 0;
  int end = 0;
  std::vector<Audio_Event> audio_events;

  int step = fft.fft_size * (1 - FFT_overlap);

  std::vector<size_t> to_rm;

  for (size_t i = 0; i < peak_locations.size(); i++)
  {
    if (peak_locations[i] < end) continue;

    Analyse analyse(audio_samples, LPF, HPF, sample_rate, step, start_t, end_t, angl_t, snr_t);
    Audio_Event audio_event = analyse.impl(FFT_size, peak_locations[i], background_noises[i], KPE, KME);

    double d = 1000 * ((double)FFT_size * (1 - FFT_overlap)) / (double)sample_rate * audio_event.duration;

    if(d < min_d || d > max_d){
      to_rm.push_back(i);
      continue;
    }

    if(audio_event.start < (int)end)
    {
      if(audio_event.duration > duration)
      {
        if(!audio_events.empty())
        {
          audio_events.pop_back();
          to_rm.push_back(i-1);
        }
      }
    }

    duration = audio_event.duration;
    end = audio_event.end;
    audio_events.push_back(audio_event);
  }

  // Check for cut signals
  size_t n_events = audio_events.size();
  std::vector<bool> del;
  del.resize(n_events, false);

  if(n_events >= 2)
  {
    for (size_t i = 1; i < n_events; i++)
    {
      if (
          ((audio_events[i].start - audio_events[i - 1].end) / (double)sample_rate * 1000) < min_TBE ||
          ((audio_events[i].start - audio_events[i - 1].end) / (double)sample_rate * 1000) > max_TBE
         )
      {
        del[i - 1] = true;
        del[i] = true;
      }
    }

    for (int i = n_events-1; i >= 0; i--)
    {
      if (del[i])
      {
        audio_events.erase(std::next(audio_events.begin(), i));
        del.erase(std::next(del.begin(), i));
      }
    }
  }

  Rcpp::CharacterVector data_names = {
    "starting_time", "duration", "freq_max_amp", "freq_max", "freq_min",
    "bandwidth", "freq_start", "freq_center", "freq_end", "freq_knee", "fc",
    "freq_bw_knee_fc", "bin_max_amp", "pc_freq_max_amp", "pc_freq_max",
    "pc_freq_min", "pc_knee", "temp_bw_knee_fc", "slope", "kalman_slope",
    "curve_neg", "curve_pos_start", "curve_pos_end", "mid_offset", "snr",
    "hd", "smoothness"
  };
  size_t n_var = data_names.size();
  n_events = audio_events.size();

  if (n_events == 0)
    return Rcpp::List();

  Rcpp::List event_data(n_var), amp_track(n_events), freq_track(n_events), event_start(n_events), event_end(n_events);
  event_data.attr("names") = data_names;
  event_data.attr("class") = "data.frame";

  event_data[0] = Rcpp::StringVector(n_events);

  Rcpp::NumericVector vec(n_events);

  for (size_t i = 1; i < n_var; ++i)
  {
    event_data[i] = clone(vec);
  }


  Rcpp::StringVector row_names(n_events);

  for(size_t i = 0; i < n_events; i++)
  {
    row_names(i) = std::to_string(i + 1);
  }

  event_data.attr("row.names") = row_names;

  Rcpp::CharacterVector out_names = { "event_data", "amp_track", "freq_track", "event_start", "event_end" };
  Rcpp::List out(out_names.size());
  out.attr("names") = out_names;
  out["event_data"] = event_data;
  out["amp_track"] = amp_track;
  out["freq_track"] = freq_track;
  out["event_start"] = event_start;
  out["event_end"] = event_end;

  for(size_t i = 0; i < n_events; i++)
  {
    extract_impl(audio_events[i], sample_rate, FFT_size, step, out, i, KPE, KME);
  }

  return(out);

}



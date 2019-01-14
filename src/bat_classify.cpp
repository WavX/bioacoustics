//------------------------------------------------------------------------------
//  Copyright (C) 2014 Chris Scott (fbscds@gmail.com)
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
#include <cmath>
#include <queue>
#include <vector>
#include <Rcpp.h>
#include "fft.h"
#include "spectro.h"
#include "tools.h"
#include "bc_img_processing.h"
#include "bc_blob_finder.h"
#include "bc_features.h"

bool sort_blobs (std::pair<int,Rcpp::List> lhs, std::pair<int,Rcpp::List> rhs);

// [[Rcpp::export]]
Rcpp::List blob_detection_impl(const std::vector<int>& audio_samples,
                               size_t sample_rate,
                               size_t FFT_size,
                               double FFT_overlap,
                               double min_TBE,
                               double max_TBE,
                               size_t HPF,
                               size_t LPF,
                               double min_d,
                               double max_d,
                               size_t area,
                               double blur_f,
                               double bg_substract,
                               double EDG,
                               double boost)
{
  size_t HPF_bin = (double)HPF * (double)FFT_size / (double)sample_rate;
  size_t LPF_bin = std::ceil((double)LPF * (double)FFT_size / (double)sample_rate) - 1;

  Rcpp::NumericMatrix spectro = \
    fspec_impl(audio_samples, FFT_size, FFT_overlap, "blackman4", HPF_bin, LPF_bin, HPF_bin, LPF_bin, false);

  int width  = spectro.ncol();
  int height = spectro.nrow();

  blur(spectro, blur_f);
  background_subtract(spectro, bg_substract);
  post_mask(spectro, EDG);
  contrast_boost(spectro, boost);

  Rcpp::NumericMatrix label (height, width);

  auto blob_map = blob_finder(spectro, label);

  std::vector<std::pair<int, Rcpp::List> > blobs(blob_map.begin(), blob_map.end());
  std::sort(blobs.begin(), blobs.end(), sort_blobs);


  size_t FFT_step = (double)FFT_size * (1 - FFT_overlap);

  std::vector<std::string> start_times;
  std::vector<double> durations;
  std::vector<int> min_offsets, max_offsets;

  start_times.reserve(blobs.size());
  durations.reserve(blobs.size());
  min_offsets.reserve(blobs.size());
  max_offsets.reserve(blobs.size());

  std::vector <Rcpp::NumericMatrix> pixel_map;

  Rcpp::CharacterVector data_names = {
    "starting_time", "duration", "area", "freq_centroid", "freq_bandwith",
    "freq_skew", "freq_kurtosis", "q", "freq_gini", "quant_2.5", "quant_25",
    "quant_50", "quant_75", "quant_97.5", "freq_bw_95_ci", "freq_bw_75_ci",
    "temp_centroid", "temp_bandwith", "temp_skew", "temp_kurtosis", "temp_gini",
    "grad_centroid", "grad_bandwith", "grad_skew", "grad_kurtosis", "grad_gini"
  };

  size_t n_var = data_names.size();

  Rcpp::List event_data (n_var);
  event_data.attr("names") = data_names;
  event_data.attr("class") = "data.frame";

  for (auto& blob: blobs)
  {
    if (Rcpp::as<size_t> (blob.second["area"]) > area)
    {
      Rcpp::NumericMatrix segment = mask(spectro, label, blob);

      double duration = ((double)segment.ncol() * FFT_step / (double)sample_rate) * 1000;

      if (duration < min_d || duration > max_d)
        continue;

      int min_offset = Rcpp::as<int> (blob.second["min_offset"]);
      int max_offset = Rcpp::as<int> (blob.second["max_offset"]);

      std::string start_time = s2dhmsms((double)min_offset * FFT_step / ((double)height * sample_rate));

      start_times.push_back(start_time);
      durations.push_back(duration);
      min_offsets.push_back(min_offset);
      max_offsets.push_back(max_offset);
      pixel_map.push_back(clone(segment));
    }
  }

  size_t n_events = min_offsets.size();
  std::vector<bool> del;
  del.resize(n_events, false);

  if(n_events >= 2)
  {
    for (size_t i = 1; i < n_events; i++)
    {
      if ( (((double)(min_offsets[i] - max_offsets[i - 1]) * FFT_step / ((double)height * sample_rate)) * 1000) < min_TBE ||
           (((double)(min_offsets[i] - max_offsets[i - 1]) * FFT_step / ((double)height * sample_rate)) * 1000) > max_TBE
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
        n_events--;
        start_times.erase(std::next(start_times.begin(), i));
        durations.erase(std::next(durations.begin(), i));
        min_offsets.erase(std::next(min_offsets.begin(), i));
        max_offsets.erase(std::next(max_offsets.begin(), i));
        pixel_map.erase(std::next(pixel_map.begin(), i));
      }
    }
  }

  if (start_times.empty())
    return Rcpp::List();

  event_data[0] = Rcpp::StringVector(n_events);

  event_data["starting_time"] = start_times;
  event_data["duration"] = durations;

  Rcpp::NumericVector vec(n_events);

  for (size_t i = 2; i < n_var; i++)
  {
    event_data[i] = clone(vec);
  }

  for (size_t i = 0; i < pixel_map.size(); i++)
  {
    calc_features(event_data, pixel_map[i], i, sample_rate, FFT_size);
  }

  Rcpp::StringVector row_names(n_events);

  for(size_t i = 0; i < n_events; i++)
  {
    row_names(i) = std::to_string(i + 1);
  }

  event_data.attr("row.names") = row_names;

  Rcpp::CharacterVector out_names = { "event_data", "pixel_map" };

  Rcpp::List out (2);
  out.attr("names") = out_names;
  out["event_data"] = event_data;
  out["pixel_map"] = pixel_map;

  return out;
}

bool sort_blobs (std::pair<int, Rcpp::List> lhs, std::pair<int, Rcpp::List> rhs)
{
  return Rcpp::as<int>(lhs.second["min_offset"]) < Rcpp::as<double>(rhs.second["min_offset"]);
}



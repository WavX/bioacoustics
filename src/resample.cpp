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

#include <vector>
#include <Rcpp.h>
#include <soxr-lsr.h>

// [[Rcpp::export]]
std::vector<float> resample_impl(std::vector<float> &audio_samples, double ratio)
{
  int converter = SRC_SINC_BEST_QUALITY;

  std::vector<float> transformed = audio_samples;
  int output_frames = (audio_samples.size() * ratio);
  transformed.resize(output_frames);

  SRC_DATA src_data;
  src_data.data_in = audio_samples.data();
  src_data.input_frames = audio_samples.size();
  src_data.data_out = transformed.data();
  src_data.output_frames = output_frames;
  src_data.src_ratio = ratio;
  src_simple(&src_data, converter, 1);
  return transformed;
}


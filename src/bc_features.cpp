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

#include <cmath>
#include <vector>
#include <Rcpp.h>
#include "bc_tools.h"

void calc_features(Rcpp::List& event_data, Rcpp::NumericMatrix& segment, const size_t& index, const size_t& sample_rate, const size_t& fft_size)
{
  int height = segment.nrow();
  int width = segment.ncol();

  std::vector<double> power_spectrum(height, 0);
  std::vector<double> temporal_envelope(width, 0);
  std::vector<double> histogram(2 * height, 0);

  size_t area = 0;
  int prev_peak = 0;

  for (int x = 0; x < width; x++)
  {
    std::vector<double> spectrum(height, 0);
    double sum = 0;

    for (int y = 0; y < height; y++)
    {
      spectrum[y] = segment(height - y-1,x);
      sum += spectrum[y];
    }

    temporal_envelope[x] = sum;
    int cur_peak = 0;

    if (sum > std::numeric_limits<double>::epsilon())
    {
      mask_spectrum(std::begin(spectrum), std::end(spectrum), 8);
      sum = 0;
      double centroid = 0;

      for (int y = 1; y < height-1; y++)
      {
        power_spectrum[y] += spectrum[y];
        area += spectrum[y] > 0.00001 ? 1 : 0;
        sum += spectrum[y];
        centroid += spectrum[y] * y;
      }

      centroid /= sum;
      cur_peak = std::round(centroid);
    }

    if (cur_peak > 0 && prev_peak > 0)
    {
      int bin = (cur_peak-prev_peak) + height;

      histogram[bin] += spectrum[cur_peak];
      histogram[bin-1] += spectrum[cur_peak-1];
      histogram[bin+1] += spectrum[cur_peak+1];
    }
    prev_peak = cur_peak;
  }

  // Segment features:
  Rcpp::as<Rcpp::NumericVector>(event_data["area"])[index] = area;

    double bin2freq = (double)sample_rate / (double)fft_size;

    // frequency
    moments(std::begin(power_spectrum), std::end(power_spectrum), event_data, std::string("freq"), index, bin2freq);

    double freq = Rcpp::as<Rcpp::NumericVector>(event_data["freq_centroid"])[index] / bin2freq;
    double bw = Rcpp::as<Rcpp::NumericVector>(event_data["freq_bandwith"])[index] / bin2freq;
    bw = bw > 1.0 ? bw : 1.0;
    double q = freq / bw;
    Rcpp::as<Rcpp::NumericVector>(event_data["q"])[index] = q * bin2freq;
    Rcpp::as<Rcpp::NumericVector>(event_data["freq_gini"])[index] = gini_impurity(std::begin(power_spectrum), std::end(power_spectrum));
    quantiles(std::begin(power_spectrum), std::end(power_spectrum), event_data, index, bin2freq);

    // temporal
    moments(std::begin(temporal_envelope), std::end(temporal_envelope), event_data, std::string("temp"), index, bin2freq);
    Rcpp::as<Rcpp::NumericVector>(event_data["temp_gini"])[index] = gini_impurity(std::begin(temporal_envelope), std::end(temporal_envelope));

    // gradient
    moments(std::begin(histogram), std::end(histogram), event_data, std::string("grad"), index, bin2freq);
    Rcpp::as<Rcpp::NumericVector>(event_data["grad_gini"])[index] = gini_impurity(std::begin(histogram), std::end(histogram));
}

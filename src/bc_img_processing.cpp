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
#include "bc_img_processing.h"

void background_subtract(Rcpp::NumericMatrix& pixels, const double& C)
{
  int height = pixels.nrow();
  int width = pixels.ncol();

  std::vector<double> mean (height, 0);

  for (int x = 1, size = 1; x < width; x++, size++)
  {
    for (int y = 0; y < height; y++)
    {
      double z = pixels(y, x) - C * mean[y];
      pixels(y, x) = (z + std::abs(z)) * 0.5;
      mean[y] = mean[y] + (z - mean[y]) / size;
    }
  }
}

void blur(Rcpp::NumericMatrix& pixels, const double& f)
{
  int height = pixels.nrow();
  int width = pixels.ncol();

  for (int x = 0; x < width; x++)
  {
    Rcpp::NumericVector spectrum = pixels.column(x);

    for (int y = 1; y < height - 1; y++)
    {
      pixels(y, x) = f * spectrum[y] + spectrum[y - 1] + spectrum[y + 1];
    }
  }
}

void contrast_boost(Rcpp::NumericMatrix& pixels, const double &boost)
{
  int height = pixels.nrow();
  int width = pixels.ncol();

  double scale = boost / (height - 7);
  Rcpp::NumericVector spectrum (height);

  for (int x = 0; x < width; x++)
  {
    pixels(0, x) = pixels(height - 1, x) = 0;
    pixels(1, x) = pixels(height - 2, x) = 0;
    pixels(2, x) = pixels(height - 3, x) = 0;
    spectrum = pixels.column(x);

    double sum { std::accumulate(spectrum.begin(), spectrum.end(), 0.0) };
    for (int y = 3; y < height - 3; y++)
    {
      pixels(y, x) -= (sum - spectrum[y] - \
        spectrum[y - 1] - spectrum[y + 1] - \
        spectrum[y + 2] - spectrum[y - 2] - \
        spectrum[y + 3] - spectrum[y - 3]) * scale;
      pixels(y, x) = pixels(y, x) < 0.0 ? 0.0 : pixels(y, x);
    }
  }
}

void post_mask(Rcpp::NumericMatrix& pixels, const double& EDG)
{
  int height = pixels.nrow();
  int width = pixels.ncol();

  double beta = 1 - EDG;
  std::vector<double> threshold(height, 0);

  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      double z = pixels(y,x);
      double decayed = EDG * threshold[y] + beta * z;
      threshold[y] = z > decayed ? z : decayed;
      pixels(y,x) = (2.0 * z) < threshold[y] ? 0.0 : z;
    }
  }
}

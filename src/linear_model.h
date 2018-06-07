//------------------------------------------------------------------------------
//  Copyright (c) 2012 Chris Scott (fbscds@gmail.com)
//  Copyright (c) 2017-2018 WavX inc. (www.wavx.ca)
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

#ifndef LINEAR_MODEL_H
#define LINEAR_MODEL_H

#include <cmath>
#include <numeric>
#include <vector>

inline
  std::vector<double> linear_model(const std::vector<int> &x, const std::vector<double> &y)
  {
    const auto n = x.size();
    const auto sX = std::accumulate(x.begin(), x.end(), 0);
    const auto sY = std::accumulate(y.begin(), y.end(), (double)0);
    const auto sXX = std::inner_product(x.begin(), x.end(), x.begin(), 0);
    const auto sXY = std::inner_product(x.begin(), x.end(), y.begin(), (double)0);
    const auto mX = (double)sX / n;
    const auto mY = sY / (double)n;
    double b = ((double)n * sXY - sX * sY) / (n * (double)sXX - std::pow(sX, 2));
    double a = mY - b * mX;

    std::vector<double> fit(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
      fit[i] = a + b * x[i];
    }

    return fit;
  }

#endif // LINEAR_MODEL_H

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

#ifndef BB_TOOLS_H
#define BB_TOOLS_H

#include <algorithm>
#include <cmath>
#include <vector>

const int ZEROLOG = -120;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double no_zero (const double &x);
inline double to_dB (const double &x);

inline
  double ang_diff (const std::vector<double> &x)
  {
    if (x.size() < 3) return 0;

    double r = x[x.size() - 1];
    double m = x[x.size() - 2];
    double l = x[x.size() - 3];
    double z = std::abs((atan2(r - m, 1) - atan2(m - l, 1)) * 180 / M_PI);
    return z;
  }


inline
void band_pass_filter(std::vector<double> &x, const size_t &LPF, const size_t &HPF, const double &freq_res)
{
  double low = LPF / freq_res;
  double high = HPF / freq_res;

  for (size_t i = 0; i < x.size(); i++)
  {
    if (low < i) {
      x[i] = 0.;
    }

    if (high > i) {
      x[i] = 0.;
    }
  }
}


inline
  double centroid (const std::vector<double> &x)
  {
    if(x.empty()) return 0;

    double index_sum = .0;
    double magnitude_sum  = .0;
    for(size_t i = 0; i < x.size(); i++)
    {
      index_sum += i * x[i];
      magnitude_sum += x[i];
    }
    return (magnitude_sum == 0) ? 0 : index_sum / magnitude_sum;
  }


template <typename T>
inline
  double clamp(T x, T min, T max) {
    return x < min ? min : x > max ? max : x;
  }


inline
  double cubic_interp (const std::vector<double> &x, const size_t i, const double &f)
  {
    return std::pow(f, 3) * (- x[i - 1] - 3 * x[i + 1] + x[i + 2] + 3 * x[i]) / 6 + \
      std::pow(f, 2) * ((x[i - 1] + x[i + 1]) / 2 - x[i]) +                         \
      f * (x[i + 1] + (-2 * x[i - 1] - (x[i + 2] + 3 * x[i])) / 6) + x[i];
  }


inline
  std::vector<double> gaussian_mask (const size_t &peak, const std::vector<double> &x)
  {
    std::vector<double> result (x.size(), 0);
    std::vector<double> kernel = {1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1};
    const double c = 4096;
      std::transform(kernel.begin(), kernel.end(), kernel.begin(), std::bind(std::divides<double>(), std::placeholders::_1, c));
    int i = peak - kernel.size() / 2;
    for (size_t j = 0; j < kernel.size(); j++) {
      i++;
      if (i >= 0 && i <  (int)x.size())
      {
        result[i] += x[i] * kernel[j];
      }
    }
    return result;
  }


inline
  double linear_interp (const double &x, const double a, const double b)
  {
    double r = a + x * (b - a);
    return r;
  }


inline
  double no_zero (const double &x)
  {
    return std::max(x, (double)0.000001);
  }


inline
  double quad_interp (std::vector<double> &x, const size_t &i)
  {
    if(i >= x.size() || i < 3)
    {
      return 0;
    }
    else
    {
      return i + (x[i + 1] - x[i - 1]) / ((double)2 * ((double)2 * x[i] - x[i - 1] - x[i + 1]));
    }
  }


inline
  double sample_at(const std::vector<double> &x, double y)
  {
    y = clamp(y, 0.0, 1.0);
    double i = y * (double)(x.size() - 1);
    double f = std::modf (i, &i);

    if (i > 0 && i < (x.size() - 2)) {
      double r = cubic_interp(x, (size_t)i, f);
      return r;
    }
    else {
      double r = ((i + 1) >= x.size() ? x[i] : linear_interp(f, x[i], x[i + 1]));
      return r;
    }
  }

inline
  void smooth_spectrum (std::vector<double> &x, const double y)
  {
    for(int i = (int)(x.size() - 2); i >= 0; i--)
    {
      x[i] = ((double)1 - y) * x[i] + y * x[i+1];
    }

    for(size_t i = 1; i < x.size(); i++)
    {
      x[i] = ((double)1 - y) * x[i] + y * x[i - 1];
    }
  }


inline
  double smoothness (const std::vector<double> &x)
  {
    double y = 0;

    for (size_t i = 1; i < ( x.size() - 1 ); i++)
    {
      y += std::abs(x[i + 1] - x[i] - x[i] + x[i - 1]);
    }
    return y;
  }


inline
  double SNR (const std::vector<double> &x, const std::vector<double> &y)
  {
    double signal = 0, noise = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
      if (x[i] > y[i])
      {
        signal += x[i];
        noise += y[i];
      }
    }
    return(to_dB(signal / no_zero(noise)));
  }


inline
  double to_dB(const double &x)
  {
    double y = 20 * std::log10(no_zero(x));
    return y;
  }


inline
  std::vector<double> up_sample(const std::vector<double> &x, const size_t &size)
  {
    std::vector<double> y (size);
    for (size_t i = 0; i < size; i++)
    {
      y[i] = sample_at(x, (double)i / (double)(size - 1));
    }
    return y;
  }

#endif // BB_TOOLS_H

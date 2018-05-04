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
//  along with This program.  If not, see <https://www.gnu.org/licenses/>.
//------------------------------------------------------------------------------

#ifndef BB_KALMAN_H
#define BB_KALMAN_H

#include <vector>

class Kalman
{
public:
  Kalman();
  Kalman(double Q, double R, double p_state_prev);

  double p_state_prev, p_error_prev;
  std::vector<double> data;

  void impl(const double &p_state_cur);
  double Q;// = 0.00001; // process
  double R;// = 0.000001; // measurement
private:

};

inline
  Kalman::Kalman(double Q, double R, double p_state_prev)
  {
    p_error_prev = 1;
    this->Q = Q;
    this->R = R;
    this->p_state_prev = p_state_prev;
  }

inline
  Kalman::Kalman()
  {
    p_error_prev = 1;
    Q = 0.00001;
    R = 0.0001;
  }

inline
  void Kalman::impl(const double &p_state_cur)
  {
    double p_error_cur = p_error_prev + Q;
    double K = p_error_cur / (p_error_cur + R);

    p_error_prev = (1 - K) * p_error_cur;
    p_state_prev = p_state_prev + K * (p_state_cur - p_state_prev);
  }


#endif // BB_KALMAN_H

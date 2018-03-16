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
  Kalman(double Q, double R, double _X);

  double _X, _P;
  std::vector<double> data;

  void impl(const double &x);
  double Q;// = 0.00001; // process
  double R;// = 0.000001; // measurement
private:

};

inline
  Kalman::Kalman(double Q, double R, double _X)
  {
    _P = 1;
    this->Q = Q;
    this->R = R;
    this->_X = _X;
  }

inline
  Kalman::Kalman()
  {
    _P = 1;
    Q = 0.00001;
    R = 0.0001;
  }

inline
  void Kalman::impl(const double &x)
  {
    double P = _P + Q;
    double K = P / (P + R);

    _P = (1 - K) * P;
    _X = _X + K * (x - _X);
  }


#endif // BB_KALMAN_H

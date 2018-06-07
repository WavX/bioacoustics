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

#ifndef IMG_PROCESSING_H
#define IMG_PROCESSING_H

#include <Rcpp.h>

void background_subtract(Rcpp::NumericMatrix& pixels, const double& C);
void blur(Rcpp::NumericMatrix& pixels, const double& f);
void contrast_boost(Rcpp::NumericMatrix& pixels, const double& boost);
void post_mask(Rcpp::NumericMatrix& pixels, const double& EDG);

#endif // IMG_PROCESSING_H

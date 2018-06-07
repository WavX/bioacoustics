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

#ifndef BC_BLOB_FINDER_H
#define BC_BLOB_FINDER_H

#include <unordered_map>
#include <Rcpp.h>

const int UNSET  =  0;
const int MARKED = -1;

std::unordered_map<int, Rcpp::List> blob_extract(Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label);
std::unordered_map<int, Rcpp::List> blob_finder(Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label);
void contour_tracing(Rcpp::NumericMatrix& mat, int offset_x, int offset_y, int blob_id, bool external_ring, std::unordered_map<int, Rcpp::List>& blob_map, Rcpp::NumericMatrix& label);
Rcpp::NumericMatrix mask(const Rcpp::NumericMatrix& mat, const Rcpp::NumericMatrix& label, std::pair<int, Rcpp::List> blob);
std::vector<int> tracer(const Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label, int offset, int p);

#endif // BC_BLOB_FINDER_H

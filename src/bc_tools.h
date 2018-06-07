//------------------------------------------------------------------------------
//  Copyright 2011-2014 Chris Scott (fbscds@gmail.com)
//  Copyright 2017-2018 WavX inc. (www.wavx.ca)
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

#ifndef BC_TOOLS_H
#define BC_TOOLS_H

#include <cmath>
#include <Rcpp.h>

template <class InputIterator>
inline
  double gini_impurity(InputIterator first, InputIterator last)
  {
    double sum = std::accumulate(first, last, 0.0);
    double gini = 1;
    if (sum > std::numeric_limits<double>::epsilon())
    {
      while (first != last)
      {
        double p = *first / sum;
        gini -= p * p;
        ++first;
      }
    }
    return gini;
  }


template <class InputIterator>
inline
  void mask_spectrum(InputIterator first, InputIterator last, int width)
  {
    auto max_it = std::max_element(first, last);

    if (max_it == last)
      return;

    auto dist = std::distance(max_it, last);

    if (dist > width)
    {
      auto it = max_it;
      std::advance(it, width);

      while (++it != last)
      {
        *it = 0;
      }
    }

    dist = std::distance(first, max_it);

    if (dist > width)
    {
      auto it = max_it;
      std::advance(it, -width);

      while (--it != first)
      {
        *it = 0;
      }
    }
  }


template <class InputIterator>
inline
  void moments(InputIterator first, InputIterator last, Rcpp::List& event_data, const std::string& which, const size_t& index, const double& bin2freq)
  {
    std::vector<double> where(std::distance(first, last));
    std::iota(std::begin(where), std::end(where), 1);

    std::vector<double> x (first, last);
    // x.assign(first, last);
    auto sum = std::accumulate(first, last, 0.0);
    auto epsilon = std::numeric_limits<double>::epsilon();
    sum = sum > epsilon ? sum : epsilon;

    std::for_each(std::begin(x), std::end(x), [sum] (double& i) { i /= sum; } );

    //	Rprintf("%5f\n", sum);
    //	return(x);
    double centroid = std::inner_product(std::begin(x), std::end(x), std::begin(where), 0.0);

    double bandwidth = 0, skew = 0, kurtosis = 0;

    for (size_t i = 0; i < x.size(); i++)
    {
      double delta = where[i] - centroid;
      double tmp = delta * delta * x[i];
      bandwidth += tmp;
      tmp *= delta;
      skew += tmp;
      tmp *= delta;
      kurtosis += tmp;
    }

    bandwidth = std::sqrt(bandwidth);
    skew = (bandwidth > epsilon) ? (skew / std::pow(bandwidth, 3.0)) : 0.0;
    kurtosis = (bandwidth > epsilon) ? (kurtosis / std::pow(bandwidth, 4.0)) : 3.0;
    kurtosis -= 3.0;

    Rcpp::as<Rcpp::NumericVector>(event_data[which + "_centroid"])[index] = centroid * bin2freq;
    Rcpp::as<Rcpp::NumericVector>(event_data[which + "_bandwith"])[index] = bandwidth * bin2freq;
    Rcpp::as<Rcpp::NumericVector>(event_data[which + "_skew"])[index] = skew * bin2freq;
    Rcpp::as<Rcpp::NumericVector>(event_data[which + "_kurtosis"])[index] = kurtosis * bin2freq;
  }


template <class InputIterator>
inline
  void quantiles(InputIterator first, InputIterator last, Rcpp::List& event_data, const size_t& index, const double& bin2freq)
  {
    std::vector<double> tmp(std::distance(first, last));
    std::partial_sum(first, last, std::begin(tmp));

    auto sum = tmp.back();
    std::vector<double> quantiles = {0.025, 0.25, 0.5, 0.75, 0.975};
    std::vector<double> freq;

    for (auto& q: quantiles)
    {
      double threshold = q * sum;
      auto it = std::find_if(std::begin(tmp), std::end(tmp),
                             [threshold] (double i) { return i >= threshold; } );
      freq.push_back(std::distance(std::begin(tmp), it) * bin2freq);
    }

    Rcpp::as<Rcpp::NumericVector>(event_data["quant_2.5"])[index] = freq[0];
    Rcpp::as<Rcpp::NumericVector>(event_data["quant_25"])[index] = freq[1];
    Rcpp::as<Rcpp::NumericVector>(event_data["quant_50"])[index] = freq[2];
    Rcpp::as<Rcpp::NumericVector>(event_data["quant_75"])[index] = freq[3];
    Rcpp::as<Rcpp::NumericVector>(event_data["quant_97.5"])[index] = freq[4];

    Rcpp::as<Rcpp::NumericVector>(event_data["freq_bw_95_ci"])[index] = freq[4] - freq[0];
    Rcpp::as<Rcpp::NumericVector>(event_data["freq_bw_75_ci"])[index] = freq[3] - freq[1];
  }


#endif // BC_TOOLS_H

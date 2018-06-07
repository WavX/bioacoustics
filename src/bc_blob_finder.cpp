//------------------------------------------------------------------------------
//  Copyright (C) 2014 Chris Scott (fbscds@gmail.com)
//  Copyright (C) 2017-2018 WavX, inc. (www.wavx.ca)
//
//  Algorithm (Chang et al. 2003)
//  "A linear-time component labeling algorithm using contour tracing technique"
//
//  Code adapted from
//  https://github.com/bramp/Connected-component-labelling/blob/master/connected-component-labelling.js
//  Andrew Brampton (2011)
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

#include <unordered_map>
#include <Rcpp.h>
#include "bc_blob_finder.h"

std::unordered_map<int, Rcpp::List> blob_extract(Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label)
{
  std::unordered_map<int, Rcpp::List> blob_map;
  int count = 1;
  int height = mat.nrow();
  int width  = mat.ncol();

  int y = 1; // We start at 1 to avoid looking above the image
  do {
    int x = 0;
    do {
      // We skip white pixels or previous labeled pixels
      if (mat(y,x) < 0.00001) continue;

      // Step 1 - P not labelled, and above pixel is white
      if (mat(y-1,x) < 0.00001 && label(y,x) == UNSET)
      {
        // P must be external contour
        contour_tracing(mat, x, y, count, true, blob_map, label);
        count++;
      }

      // Step 2 - Below pixel is white, and unmarked
      if (mat(y+1,x) < 0.00001 && label(y+1,x) == UNSET)
      {
        // Use previous pixel label, unless this is already labelled
        int n = label(y-1, x);
        if (label(y,x) != UNSET) n = label(y,x);

        // P must be a internal contour
        contour_tracing(mat, x, y, n, false, blob_map, label);
      }

      // Step 3 - Not dealt with in previous two steps
      if (label(y,x) == UNSET)
      {
        int n = label(y-1,x);
        // Assign P the value of N
        label(y,x) = n;

        if (blob_map.find(n) == blob_map.end())
        {
          double magnitude = 0;
          size_t area = 0;
          blob_map[n] = Rcpp::List::create(Rcpp::Named("magnitude") = magnitude,
                                           Rcpp::Named("area") = area,
                                           Rcpp::Named("min_offset") = std::numeric_limits<int>::max(),
                                           Rcpp::Named("max_offset") = 0,
                                           Rcpp::Named("blob_id") = n);
        }

        // Update blob_map
        blob_map[n]["magnitude"] = mat(y, x) + Rcpp::as<double> (blob_map[n]["magnitude"]);
        blob_map[n]["area"] = Rcpp::as<size_t> (blob_map[n]["area"]) + 1;

        int offset = x * height + y;

        if (offset < Rcpp::as<int>(blob_map[n]["min_offset"]))
        {
          blob_map[n]["min_offset"] = offset;
        }
        if (offset > Rcpp::as<int>(blob_map[n]["max_offset"]))
        {
          blob_map[n]["max_offset"] = offset;
        }
      }
    } while (++x < width);
  } while (++y < height-1); // We end one before the end to to avoid looking below the spectrogram
  return blob_map;
}

std::unordered_map<int, Rcpp::List> blob_finder(Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label)
{
  int height = mat.nrow();
  int width = mat.ncol();

  // We change the border to be white. We could add a pixel around
  // but we are lazy and want to do this in place.
  // Set the outer rows/cols to background
  for (int y = 0; y < height; y++)
  {
    mat(y, 0) = 0;
    mat(y, width - 1) = 0;
  }

  for (int x = 0; x < width; x++)
  {
    mat(0, x) = 0;
    mat(height - 1, x) = 0;
  }

  std::unordered_map<int, Rcpp::List> blob_map = blob_extract(mat, label);
  return blob_map;
}

void contour_tracing(Rcpp::NumericMatrix& mat,
                     int offset_x,
                     int offset_y,
                     int blob_id,
                     bool external_ring,
                     std::unordered_map<int, Rcpp::List>& blob_map,
                     Rcpp::NumericMatrix& label)
{
  int p = external_ring ? 7 : 3;

  int height = mat.nrow();
  int offset = offset_x * height + offset_y;

  // Find out our default next pos (from offset)
  std::vector<int> traced = tracer(mat, label, offset, p);
  int default_offset = traced[0];
  int q              = traced[1];

  if (blob_map.find(blob_id) == blob_map.end())
  {
    double magnitude = 0;
    size_t area = 0;
    blob_map[blob_id] = Rcpp::List::create(Rcpp::Named("magnitude") = magnitude,
                                           Rcpp::Named("area") = area,
                                           Rcpp::Named("min_offset") = std::numeric_limits<int>::max(),
                                           Rcpp::Named("max_offset") = 0,
                                           Rcpp::Named("blob_id") = blob_id);
  }

  label(offset_y, offset_x) = blob_id;

  // Update blob_map
  blob_map[blob_id]["magnitude"] = mat(offset_y, offset_x) + Rcpp::as<double> (blob_map[blob_id]["magnitude"]);
  blob_map[blob_id]["area"] = Rcpp::as<size_t>(blob_map[blob_id]["area"]) + 1;

  if (offset < Rcpp::as<int>(blob_map[blob_id]["min_offset"]))
  {
    blob_map[blob_id]["min_offset"] = offset;
  }
  if (offset > Rcpp::as<int>(blob_map[blob_id]["max_offset"]))
  {
    blob_map[blob_id]["max_offset"] = offset;
  }

  // Single pixel check
  if (default_offset == offset) return;

  int stored_offset = default_offset;
  int new_offset    = default_offset;

  while ( new_offset != offset || stored_offset != default_offset )
  {
    int stored_offset_x = stored_offset / height;
    int stored_offset_y = stored_offset % height;
    label(stored_offset_y, stored_offset_x) = blob_id;

    // Update blob_map
    blob_map[blob_id]["magnitude"] = mat(stored_offset_y, stored_offset_x) + Rcpp::as<double>(blob_map[blob_id]["magnitude"]);
    blob_map[blob_id]["area"] = Rcpp::as<size_t>(blob_map[blob_id]["area"]) + 1;

    if (stored_offset < Rcpp::as<int>(blob_map[blob_id]["min_offset"]))
    {
      blob_map[blob_id]["min_offset"] = stored_offset;
    }
    if (stored_offset > Rcpp::as<int>(blob_map[blob_id]["max_offset"]))
    {
      blob_map[blob_id]["max_offset"] = stored_offset;
    }

    new_offset = stored_offset;
    p = (q + 5) % 8;
    std::vector<int> traced = tracer(mat, label, new_offset, p);
    stored_offset = traced[0];
    q             = traced[1];
  }
}

Rcpp::NumericMatrix mask(const Rcpp::NumericMatrix& mat, const Rcpp::NumericMatrix& label, std::pair<int, Rcpp::List> blob)
{
  int height = mat.nrow();
  int min = Rcpp::as<int>(blob.second["min_offset"]) / height;
  int max = Rcpp::as<int>(blob.second["max_offset"]) / height;

  int n_col = max - min + 1;

  Rcpp::NumericMatrix segment (height, n_col);

  for (int x = min, xs = 0; x <= max; x++, xs++)
  {
    for (int y = 0; y < height; y++)
    {
      if (label(y,x) == blob.first)
      {
        segment(y, xs) = mat(y, x);
      }
    }
  }

  return segment;
}


std::vector<int> tracer(const Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& label, int offset, int p)
{
  std::vector<int> pos(8);
  int height = mat.nrow();
  int width = mat.ncol();

  pos[0] = 1;         // below
  pos[1] = height+1;  // below-right
  pos[2] = height;    // right
  pos[3] = height-1;  // above-right
  pos[4] = -1;        // above
  pos[5] = -height-1; // above-left
  pos[6] = -height;   // left
  pos[7] = -height+1; // below-left

  int size = height * width;

  for (int d = 0; d < 8; d++)
  {
    int q = (p + d) % 8;
    int new_offset = offset + pos[q];

    // Make sure we are inside the spectrogram
    if (new_offset < 0 || new_offset >= size) continue;

    int offset_x = new_offset / height;
    int offset_y = new_offset % height;

    if (mat(offset_y, offset_x) >= 0.00001)
    {
      return std::vector<int> {new_offset, q};
    }

    label(offset_y, offset_x) = MARKED;
  }

  // No move
  return std::vector<int> {offset, -1};
}

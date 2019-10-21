//------------------------------------------------------------------------------
// Copyright (C) 2014-2017 Wildlife Acoustics, inc.
// Copyright (C) 2017-2018 WavX, inc. (www.wavx.ca)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//------------------------------------------------------------------------------

#ifndef READ_WAC_H
#define READ_WAC_H

#include <string>
#include <cstring>
#include <deque>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Rcpp.h>

struct wac_s
{
  int            version;       // WAC file version number
  int            flags;         // WAC flags
  int            framesize;     // samples per channel per frame
  int            blocksize;     // frames per block
  int            seeksize;      // blocks per seek-table entry
  int            seekentries;   // number of seek-table entries
  int            channelcount;  // number of channels
  int            samplerate;    // sample rate
  unsigned long  sample_read = 0;   // number of samples already decoded
  unsigned long  samplecount;   // number of samples in file per channel

  FILE           *fp;           // input file descriptor
  int            frameindex;    // current frame index
  int            filebit_index; // index to current bit in word (0-15)
  unsigned short bitbuffer;     // remaining bits in bit buffer

  int it;

  // TRIGGER AWARE BEHAVIOUR BEGINS //
  std::deque< std::vector<signed short> > left; // data container (left channel) (accomodate for TRIGGER wac)
  std::deque< std::vector<signed short> > right; // data container (right channel) (accomodate for TRIGGER wac)
  std::vector<signed short> left_tmp;
  std::vector<signed short> right_tmp;

  std::vector<size_t> trigger;

  // TRIGGER AWARE BEHAVIOUR ENDS //

  // ORIGINAL BEHAVIOUR BEGINS //
  // std::vector<signed short> left;
  // std::vector<signed short> right;
  // ORIGINAL BEHAVIOUR ENDS //
};



int ReadBits(wac_s *w, int _bits, const std::string &filename);
unsigned short ReadWord(wac_s *w, const std::string &filename);
void FrameDecode(wac_s *w, const std::string &filename);

Rcpp::List read_wac_impl (const std::string filepath, const std::string filename);


#endif // READ_WAC_H

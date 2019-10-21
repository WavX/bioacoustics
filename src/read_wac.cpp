//------------------------------------------------------------------------------
// Copyright (C) 2014-2017 Wildlife Acoustics, inc.
// Further modified at WavX, inc. (www.wavx.ca) since 2017
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

// A WAC file has the following format.  Note multi-byte values are
// little-endian.
//
//   1. WAC HEADER (24 bytes)
//      0x00 - 0x03 = "WAac" - identifies this as a WAC file
//      0x04        = Version number (<= 4)
//      0x05        = Channel count (1 for mono, 2 for stereo)
//      0x06 - 0x07 = Frame size = # samples per channel per frame
//      0x08 - 0x09 = Block size = # frames per block
//      0x0a - 0x0b = Flags
//                      0x0f = mask of "lossy bits".  For "WAC0", this is 0.
//                             For "WAC1", this is 1, and so on, representing
//                             increasing levels of compression as the number
//                             of discarded least-significant bits.  WAC0 is
//                             lossless compression, WAC1 is equivalent to
//                             15-bit dynamic range, WAC2 is equivalent to
//                             14-bit dynamic range, and so on.
//                      0x10 = TRIGGER WAC file e.g. one or both channels
//                             have triggers.  What this means is that highly
//                             compressed zero-value frames may be inserted
//                             in the data stream representing unTRIGGER
//                             time between TRIGGER recordings.  Software
//                             capable of handling triggers can break a file
//                             into pieces discarding these zero-value frames.
//                             Note: This example code does not support
//                             this, but look in the code for "ZERO FRAME"
//                             to see where this would be used.
//                      0x20 = GPS data present - GPS data is interleaved with
//                             data in block headers.
//                      0x40 = TAG data present - TAG data is interleaved with
//                             data in block headers.  The TAG corresponds to
//                             an EM3/EM3+ button press to tag a recording.
//      0x0c - 0x0f = Sample rate (samples per second)
//      0x10 - 0x13 = Sample count (number of samples in WAC file per channel)
//      0x14 - 0x15 = Seek size (number of blocks per seek table entry)
//      0x16 - 0x17 = Seek entries (size of seek table in 32-bit words)
//
// 2. SEEK TABLE
//    The Seek Table contains (Seek entries) number of 4-byte (32-bit)
//    values representing the absolute offset into the WAC file corresponding
//    to each (Seek size) blocks.  The offset is measured in 16-bit words so
//    you would double these values to convert to a byte offset into the file.
//    The intention of the seek table is to make it easier to jump to a position
//    in the WAC file without needing to decompress all the data before that
//    position.  This code example does not use the seek table so we simply
//    skip over it.
//
// 3. BLOCKS OF FRAMES OF SAMPLES
//    Samples are grouped into frames (according to the frame size), and
//    frames are organized into blocks (according to the block size).
//    Additionally, blocks are organized into seek table entries as described
//    above according to the seek size.
//
//    Each block is aligned to a 16-bit boundary and consists of a block
//    header followed by block size frames. The format of the block header is
//    as follows:
//
//    0x00 - 0x03 = 0x00018000 = unique block header pattern
//    0x04 - 0x07 = block index (starting with zero and incrementing by one
//                  for each subsequent block used to keep things synchronized
//                  and detect file corruption.  This is also convenient for
//                  seeking to a particular block as the patterns here will
//                  not occur in the data stream.
//
//    Following the block header are a series of variable-length bit-fields
//    which do not necessarily line up on byte boundaries.  Refer to the
//    ReadBits() function for specifics relating to the encoding.
//
//    If (flags & 0x20), then GPS data is present in every seek size blocks
//    beginning with the first block at index zero.  The GPS data is encoded as
//    25-bits of signed latitude and 26-bits of signed longitude information.
//    (using 2's complement notation). The latitude and longitude values in
//    degrees can be determined by dividing these signed quantities by 100,000
//    with positive values corresponding to North latitude and West longitude.
//
//    If (flags & 0x40), then tag data is present in every block and is
//    represented by 4-bits.  For tagged recordings (e.g. from an EM3), the
//    tag values 1-4 correspond to the buttons 'A' through 'D', and a value 0
//    indicates that no tag is present.  While the tag button is pressed,
//    blocks will be written with the corresponding tag.
//
//    Note that the GPS and TAG values are not implemented in this code and are
//    simply skipped, but please see the comments in the code for more
//    information.
//
//    Following the block header and optional GPS or tag data are block size
//    frames of frame size samples for each channel.  For multi-channel
//    recordings, samples are interleaved.
//
//    Compression uses Golumb coding of the deltas between successive samples.
//    The number of bits used to represent the remainder is variable and
//    optimized for each frame and for each channel.  The quotient is
//    represented by alternating 1/0 bits ahead of the remainder.
//
//    The frame begins with a 4-bit value for each channel indicating the
//    number of bits used to represent the remainder.  Note that a zero value
//    indicates that the frame contains no content e.g. representing the
//    space inbetween TRIGGER recordings.
//
//    What follows are Golumb-encoded representations of deltas of interleaved
//    (by channel) samples.  Details can be found in FrameDecode().
//
//    NOTE: We have not yet added Wildlife Acoustics metadata to the WAC
//    format and may do so in the future, quite likely by appending a
//    "Wamd" chunk at the end of the file.

#include "read_wac.h"
#include "tools.h"

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List read_wac_impl (const std::string filepath,
                          const std::string filename)
{
  size_t sz;
  unsigned char hdr[24];

  // Initialize WAC state
  wac_s ws;
  ws.frameindex = 0;
  ws.it = 0;
  ws.fp = std::fopen(filepath.c_str(), "rb");

  if (ws.fp == NULL)
  {
    Rcpp::Rcerr << "Error reading " << filename << "\n";
  }

  // Parse WAC header and validate supported formats
  if ((sz = std::fread(hdr, 1, 24, ws.fp)) != 24)
  {
    Rcpp::Rcerr << filename << ": Unexpected eof \n";
  }

  // Verify "magic" header
  if (   hdr[0] != 'W'
      || hdr[1] != 'A'
      || hdr[2] != 'a'
      || hdr[3] != 'c'
  )
  {
    Rcpp::Rcerr << filename << " is not a WAC file \n";
  }

  // Check version
  ws.version = hdr[4];
  if (ws.version > 4)
  {
    Rcpp::Rcerr << filename << " version (" << std::to_string(ws.version) << ") not supported\n";
  }

  // Read channel count and frame size
  ws.channelcount = hdr[5];
  ws.framesize = hdr[6] | (hdr[7] << 8);

  // All Wildlife Acoustics WAC files have 512-byte (256-sample mono or
  // 128 sample stereo) frames.
  if (ws.channelcount * ws.framesize != 256)
  {
    Rcpp::Rcerr << filename << ": Unsupported block size (" << std::to_string(ws.channelcount * ws.framesize) << ")\n";
  }

  // All Wildlife Acoustics WAC files have 1 or 2 channels
  if (ws.channelcount > 2)
  {
    Rcpp::Rcerr << filename << ": Unsupported channel count (" << std::to_string(ws.channelcount) << ")\n";
  }

  // Read flags
  ws.flags = hdr[10] | (hdr[11] << 8);

  // Parse additional fields from the WAC header
  ws.blocksize = hdr[8] | (hdr[9] << 8);
  ws.samplerate = hdr[12] | (hdr[13] << 8) | (hdr[14] << 16) | (hdr[15] << 24);
  ws.samplecount = hdr[16] | (hdr[17] << 8) | (hdr[18] << 16) | (hdr[19] << 24);
  ws.seeksize = hdr[20] | (hdr[21] << 8);
  ws.seekentries = hdr[22] | (hdr[23] << 8);

  // TRIGGER AWARE BEHAVIOUR BEGINS //
  if (ws.flags & 0x10)
  {
    Rcpp::Rcout << "WAC with trigger(s), converting individual segment(s) to WAV file(s)." << std::endl;
  }
  // TRIGGER AWARE BEHAVIOUR ENDS //

  // Use the seek table to find first data block
  unsigned short *b;
  b = (unsigned short*) std::malloc (sizeof(unsigned short)*2);
  if ((sz = std::fread(b, 2, 1, ws.fp)) != 1)
  {
    Rcpp::Rcerr << filename << ": Unexpected eof\n";
  }
  *b *= 2;
  // Position the cursor at the first data block
  std::fseek(ws.fp, *b, SEEK_SET);
  free(b);

  // Read frames of data and store samples in memory
  while (ws.samplecount > 0)
  {
    FrameDecode(&ws, filename);
    ws.samplecount -= ws.framesize;
    ws.sample_read += ws.framesize;
  }

  std::fclose(ws.fp);

  return Rcpp::List::create(Rcpp::Named("version") = ws.version,
                            Rcpp::Named("sample_rate") = ws.samplerate,
                            Rcpp::Named("bits") = 16,
                            Rcpp::Named("left") = ws.left,
                            Rcpp::Named("right") = ws.right,
                            Rcpp::Named("trigger") = ws.trigger);
}


// Read bits:
//
// This function reads the specified number of bits from the file (1-16) and
// returns the signed value.
//
// We buffer 16-bits at a time in bitbuffer and keep track of the current
// filebit_index number of bits remaining in the buffer.
//
int ReadBits(wac_s *w, int bits, const std::string &filename)
{
  unsigned long x = 0;
  std::size_t sz;

  // While we need bits...
  while (bits > 0)
  {
    // If starting a new 16-bit word, read the next 16 bits
    if (w->filebit_index == 0)
    {
      if ((sz = fread(&w->bitbuffer, 2, 1, w->fp)) != 1)
      {
        Rcpp::Rcerr << filename << ": Unexpected eof\n";
      }

      w->filebit_index = 16;
    }

    // If all the bits we need are in the current word, extract the bits,
    // update state, and break out of the loop.
    if (bits < w->filebit_index)
    {
      x <<= bits;
      x |= w->bitbuffer >> (16 - bits);
      w->bitbuffer <<= bits;
      w->filebit_index -= bits;
      break;
    }
    else
    {
      // Otherwise extract the bits we have and continue (which will
      // then load the next 16-bits into the bitbuffer
      x <<= w->filebit_index; // no need ?
      x |= w->bitbuffer >> (16 - w->filebit_index);
      bits -= w->filebit_index;
      w->filebit_index = 0;
    }
  }
  return (int) x;
}

// ReadWord
//
// Skip to 16-bit boundary and read and return the next 16 bits.
//
unsigned short ReadWord(wac_s *w, const std::string &filename)
{
  w->filebit_index = 0;
  return ReadBits(w, 16, filename);
}

// FrameDecode
//
// Decode the next frame and write interleaved 16-bit samples to output file
//
void FrameDecode(wac_s *w, const std::string &filename)
{
  unsigned short s;
  int i;
  int ch;
  unsigned short code;
  short lastsample[2];
  int g[2];
  int lossybits = w->flags & 0x0f;

  // At start of block parse block header
  if (0 == (w->frameindex % w->blocksize))
  {
    // Verify that the block header is valid and as expected
    int block = w->frameindex / w->blocksize;

    if (   ReadWord(w, filename) != 0x8000
        || ReadWord(w, filename) != 0x0001
        || ReadWord(w, filename) != (block & 0xffff)
        || ReadWord(w, filename) != ((block >> 16) & 0xffff)
    )
    {
      Rcpp::Rcerr << filename << ": Bad block header\n";
    }

    // If GPS data present and the block number is modulo the blocks per
    // seek table entry, then load the latitude and longitude.
    if ((w->flags & 0x20) && 0 == (block % w->seeksize))
    {
      // int lathi = ReadBits(w,9);
      // int latlo = ReadBits(w,16);
      // int lonhi = ReadBits(w,10);
      // int lonlo = ReadBits(w,16);

      // For this example, we are not actually using the latitude and
      // longitude.  But these values are 100,000 times the latitude and
      // longitude in degrees with positive sign indicating north latitude
      // and west longitude.  These can be derived as follows:
      //
      //   float lat = ((lathi<<16)|latlo) / 100000.0;
      //   float lon = ((lonhi<<16)|lonlo) / 100000.0;
    }

    // If tag data present read it.  For this example, we aren't using this
    // data.  But for EM3 recordings, for example, this tag value would
    // be 1-4 corresponding to buttons A-D being depressed during this frame.
    if (w->flags & 0x40)
    {
      // int tag = ReadBits(w,4);
    }
  }
  // Advance frame
  w->frameindex++;

  // Read the per-channel Golumb remainder code size
  for (ch = 0; ch < w->channelcount; ch++)
  {
    lastsample[ch] = 0;
    g[ch] = ReadBits(w, 4, filename);
  }

  // Read samples for frame
  for (i = 0; i < w->framesize; i++)
  {
    // Interleave channels
    for (ch = 0; ch < w->channelcount; ch++)
    {
      // TRIGGER AWARE BEHAVIOUR BEGINS //
      std::deque< std::vector<signed short> > *channel;
      std::vector<signed short> *tmp;

      if (ch == 0)
      {
        channel = &w->left;
        tmp = &w->left_tmp;
      }
      else
      {
        channel = &w->right;
        tmp = &w->right_tmp;
      }
      // TRIGGER AWARE BEHAVIOUR ENDS //

      // ORIGINAL BEHAVIOUR BEGINS //
      // std::vector<signed short> *channel;
      // if (ch == 0) {
      //   channel = &w->left;
      // } else {
      //   channel = &w->right;
      // }
      // ORIGINAL BEHAVIOUR ENDS //


      int delta;
      int stopbit;

      // ZERO FRAME
      // Special case: For TRIGGER WAC files, the code size is set to
      // zero to indicate a zero-value frame.  This represents the
      // space between TRIGGER events in the WAC file.  A trigger-aware
      // parser would look for onset/offset of this condition to break up
      // a WAC file into individual triggers.  In this example, we are just
      // filling the space with zero sample values.  And we won't actually
      // get to this code because the presence of these special zero-value
      // frames also requires that the header flags has 0x10 set.  In
      // main() above, we exit with an error if this is the case.

      // TRIGGER AWARE BEHAVIOUR BEGINS //
      if (g[ch] == 0)
      {
        size_t s = tmp->size();

        if(s > 0)
        {
          channel->push_back(*tmp);
          tmp->clear();
          w->trigger.push_back(w->sample_read - s);
        }
        return;
      }
      // TRIGGER AWARE BEHAVIOUR ENDS //

      // ORIGINAL BEHAVIOUR BEGINS //
      // if (g[ch] == 0)
      // {
      //   channel->push_back(0);
      //   continue;
      // }
      // ORIGINAL BEHAVIOUR ENDS //

      // Read the remainder code from the specified number of bits
      code = ReadBits(w, g[ch], filename);

      // Read the quotient represented by alternating 1/0 pattern
      // following the remainder code
      stopbit = (code & 1) ^ 1;
      while (stopbit != ReadBits(w, 1, filename))
      {
        code += 1 << g[ch];
        stopbit ^= 1;
      }

      // Adjust for sign
      if (code & 1)
      {
        delta = -(code+1)/2;
      }
      else
      {
        delta = code/2;
      }

      // Compute sample value as delta from previous sample
      delta += lastsample[ch];
      lastsample[ch] = delta;

      // Restore dropped least-significant bits used in higher levels
      // of compression e.g. WAC1, WAC2, etc.
      delta <<= lossybits;

      // Write 16-bit sample to output file
      s = delta;

      // TRIGGER AWARE BEHAVIOUR BEGINS //
      if (w->flags & 0x10)
      {
        tmp->push_back(s);
      }
      else
      {
        if (channel->size() == 0)
        {
          channel->push_back(*tmp);
        }
        channel->at(0).push_back(s);
      }
      // TRIGGER AWARE BEHAVIOUR ENDS //

      // channel->push_back(s); // ORIGINAL BEHAVIOUR
    }
  }
}


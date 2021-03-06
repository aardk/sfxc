/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * Copyright (c) 2007 University of Amsterdam (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *            Damien Marchal <dmarchal@science.uva.nl>, 2007
 *
 *
 * This file is part of:
 *   - SFXC/SCARIe project.
 * This file contains:
 *   - Implementation of a channel extractor with look up table.
 */

#include "channel_extractor_5.h"
#include "utils.h"

#include <iostream>
#include <cstring>

template<int size_of_one_input_word_, int n_subbands_>
class Channel_extractor_5_impl : public Channel_extractor_interface {
public:
  Channel_extractor_5_impl() {}

  void initialise(const std::vector< std::vector<int> > &track_positions,
                  int IGNORED_size_of_one_input_word,
                  int input_sample_size, 
                  int bits_per_sample_) {

    input_sample_size_ = input_sample_size;
    //n_subbands = track_positions.size();
    fan_out = track_positions[0].size();
    samples_per_byte = 8/fan_out;
    bits_per_sample = bits_per_sample_;
    SFXC_ASSERT(fan_out <= 8);
    SFXC_ASSERT(8%fan_out == 0);
    SFXC_ASSERT(fan_out*samples_per_byte == 8);

    memset(lookup_table, 0, size_of_one_input_word_*256*n_subbands_);
    // Lookup table for the tracks
    for (int subband=0; subband < n_subbands_; subband++) {
      int low_bit = 1<<(8-fan_out);
      for (int track_nr=0; track_nr < fan_out; track_nr++) {
        int track = track_positions[subband][track_nr];
        int n = track/8;
        SFXC_ASSERT(n < size_of_one_input_word_);
        int bit = track % 8;
        for (int sample=0; sample<256; sample++) {
          if (((sample >> bit)& 1) != 0) {
            lookup_table[subband][n][sample] |= low_bit<<((track_nr+1)%bits_per_sample);
	    //lookup_table[subband][n][sample] |= (1<<(fan_out-1-track_nr));
          }
        }
        if((track_nr+1)%bits_per_sample==0)
          low_bit=low_bit << bits_per_sample;
      }
    }
  }

  void extract(unsigned char *in_data1,
               unsigned char **output_data) {
    do_task(input_sample_size_, in_data1, output_data);
  };


  void do_task(int n_input_samples,
               unsigned char in_data[],
               unsigned char** output_data) {

    for (int subband=0; subband<n_subbands_; subband++) {
      const uint8_t *in_pos = (const uint8_t *)in_data;
      unsigned char *out_pos = &(output_data[subband][0]);
      const uint8_t *table = &lookup_table[subband][0][0];
      memset((void *)out_pos, 0, n_input_samples*fan_out/8);

      for (int sample=0; sample<n_input_samples; sample+=samples_per_byte) {
        // samples_per_byte-1 times:
        for (int i=1; i<samples_per_byte; i++) {
          process_sample(in_pos, out_pos, table);
          *out_pos = (*out_pos >> fan_out);
        }

        process_sample(in_pos, out_pos, table);
        out_pos++;
      }
    }
  }

private:
  inline void process_sample(const uint8_t *&in_array,
                             unsigned char *output_data,
                             const uint8_t *table) {
    for (int n=0; n<size_of_one_input_word_; n++) {
      *output_data |= table[*in_array];
      table += 256;
      in_array++;
    }
  }
  // Lookup table for the channel extraction
  // At most 8 channels
  // lookup_table[index in Type word][value of the byte][output sample per channel]
  uint8_t lookup_table[n_subbands_][size_of_one_input_word_][256];

  // Fan out
  int fan_out, samples_per_byte, bits_per_sample;

  int input_sample_size_;
};

Channel_extractor_5::Channel_extractor_5() {
  name_ = "Channel_extractor_5";
  hidden_implementation_ = NULL;
}


template<int Tsize_of_word>
Channel_extractor_interface* create_number5_(int n_subbands) {
  switch (n_subbands) {
    case 1:
    return new Channel_extractor_5_impl<Tsize_of_word, 1>();
    break;
    case 2:
    return new Channel_extractor_5_impl<Tsize_of_word, 2>();
    break;
    case 4:
    return new Channel_extractor_5_impl<Tsize_of_word, 4>();
    break;
    case 8:
    return new Channel_extractor_5_impl<Tsize_of_word, 8>();
    break;
    case 16:
    return new Channel_extractor_5_impl<Tsize_of_word, 16>();
    break;
    case 32:
    return new Channel_extractor_5_impl<Tsize_of_word, 32>();
    break;
    default:
    SFXC_ASSERT_MSG(false,
                    "Unsupported number of channel for this channelizer");
  }
  return NULL;
}

Channel_extractor_interface* create_number5_(int size_of_one_input_word, int n_subbands) {
  switch (size_of_one_input_word) {
    case 1:
    return create_number5_<1>(n_subbands);
    break;
    case 2:
    return create_number5_<2>(n_subbands);
    break;
    case 4:
    return create_number5_<4>(n_subbands);
    break;
    case 8:
    return create_number5_<8>(n_subbands);
    break;
    default:
    std::cerr << " size_of_input_word: " << size_of_one_input_word << std::endl;
    SFXC_ASSERT_MSG(false,
                    "Unsupported size of input word for channelization");
  }
  return NULL;
}


void Channel_extractor_5::initialise(const std::vector< std::vector<int> > &track_positions_,
                                     int size_of_one_input_word_,
                                     int input_sample_size_, int bits_per_sample_) {
  track_positions = track_positions_;
  size_of_one_input_word = size_of_one_input_word_;
  input_sample_size = input_sample_size_;
  fan_out = track_positions[0].size();
  n_subbands = track_positions.size();

  hidden_implementation_ = create_number5_(size_of_one_input_word, n_subbands);
  if ( hidden_implementation_ == NULL ) {
    std::cout << "UNABLE TO CREATE EXTRACTOR5 " << std::endl;
    exit(0);
  }
  hidden_implementation_->initialise(track_positions, size_of_one_input_word_, input_sample_size, bits_per_sample_);
}

void Channel_extractor_5::extract(unsigned char *in_data1,
                                  unsigned char **output_data) {
  SFXC_ASSERT( hidden_implementation_ != NULL && "INVALID EXTRACTOR NULL POINTER" );
  hidden_implementation_->extract(in_data1, output_data);
}

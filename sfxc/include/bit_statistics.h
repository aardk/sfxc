/* Copyright (c) 2010 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Aard Keimpema, <keimpema@jive.nl>
 *
 */
#ifndef BIT_STATISTICS_H
#define BIT_STATISTICS_H
#include <vector>
#include <boost/shared_ptr.hpp>
#include "utils.h"

class bit_statistics;
typedef boost::shared_ptr<bit_statistics>   bit_statistics_ptr;

class bit_statistics{
public:
  bit_statistics();
  ~bit_statistics();
  void reset_statistics(int bits_per_sample_, uint64_t sample_rate_, uint64_t base_sample_rate_);
  void inc_counter(unsigned char word, bool);
  void inc_invalid(int n);
  int64_t *get_statistics();
  int64_t *get_tsys();
  int bits_per_sample;
  uint64_t sample_rate;
  uint64_t base_sample_rate;
private:
  int64_t nInvalid;
  std::vector<int64_t> data_counts_on;
  std::vector<int64_t> data_counts_off;
  std::vector<int64_t> statistics;
  std::vector<int64_t> tsys;
};

inline void 
bit_statistics::inc_counter(unsigned char word, bool on){
  if (on)
    data_counts_on[(unsigned int)word]++;
  else
    data_counts_off[(unsigned int)word]++;
}

inline void 
bit_statistics::inc_invalid(int n){
  nInvalid+=n;
}
#endif

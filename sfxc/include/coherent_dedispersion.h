#ifndef COHERENT_DEDISPERSION_H
#define COHERENT_DEDISPERSION_H
#include "utils.h"
#include "control_parameters.h"
#include "correlator_node_types.h"
#ifdef USE_DOUBLE
#include "sfxc_fft.h"
#else
#include "sfxc_fft_float.h"
#endif

class Coherent_dedispersion {
public:
  typedef Correlator_node_types::Delay_memory_pool   Delay_memory_pool;
  typedef Correlator_node_types::Delay_queue         Delay_queue;
  typedef Correlator_node_types::Delay_queue_ptr     Delay_queue_ptr;
  typedef Delay_queue::value_type                    Delay_queue_element;
  typedef Pulsar_parameters::Pulsar                  Pulsar;

  Coherent_dedispersion(int stream_nr_, SFXC_FFT &fft_, SFXC_FFT &fft_cor_, 
                      Memory_pool_vector_element<std::complex<FLOAT> > &filter_,
                      Memory_pool_vector_element<std::complex<FLOAT> > &dispersion_buffer_,
                      Memory_pool_vector_element<FLOAT> &zeropad_buffer_ );
  ~Coherent_dedispersion();
  void do_task();
  bool has_work();
  void set_parameters(const Correlation_parameters &parameters);
  void connect_to(Delay_queue_ptr buffer);
  void empty_output_queue();
  /// Get the output
  Delay_queue_ptr get_output_buffer();
private:
  void allocate_element(int nfft);
  void overlap_add();
  size_t fft_rot_size();
  size_t fft_cor_size();
  size_t fft_dedisp_size();
  uint64_t sample_rate();
private:
  int stream_nr;
  int stream_idx;
  int out_pos;
  int fft_to_skip;
  int current_fft, current_buffer;
  int number_ffts_in_slice;
  Correlation_parameters correlation_parameters;
  int output_stride;
  Memory_pool_vector_element<std::complex<FLOAT> > &filter, &dedispersion_buffer;
  Memory_pool_vector_element<FLOAT> time_buffer[2], &zeropad_buffer;
  Delay_queue_element cur_output;

  Delay_memory_pool output_memory_pool;
  Delay_queue_ptr input_queue;
  Delay_queue_ptr output_queue;
  Time start_time;
  SFXC_FFT  &fft, &fft_cor;
};

inline size_t Coherent_dedispersion::fft_dedisp_size() {
  return (2 * correlation_parameters.fft_size_dedispersion * sample_rate()) / 
         correlation_parameters.sample_rate;
}

inline size_t Coherent_dedispersion::fft_rot_size() {
  return (fft_cor_size() * sample_rate()) / correlation_parameters.sample_rate;
}

inline size_t Coherent_dedispersion::fft_cor_size() {
  return 2 * correlation_parameters.fft_size_correlation;
}

inline uint64_t Coherent_dedispersion::sample_rate() {
  return correlation_parameters.station_streams[stream_idx].sample_rate;
}
#endif

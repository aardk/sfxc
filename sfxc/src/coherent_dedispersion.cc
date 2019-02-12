#include "coherent_dedispersion.h"

Coherent_dedispersion::Coherent_dedispersion(int stream_nr_,
                         SFXC_FFT  &fft_, 
                         SFXC_FFT  &fft_cor_, 
                         Memory_pool_vector_element<std::complex<FLOAT> > &filter_,
                         Memory_pool_vector_element<std::complex<FLOAT> > &dedispersion_buffer_,
                         Memory_pool_vector_element<FLOAT> &zeropad_buffer_): 
      output_queue(Delay_queue_ptr(new Delay_queue())),
      stream_nr(stream_nr_), output_memory_pool(2, NO_RESIZE), 
      current_buffer(0), fft(fft_), fft_cor(fft_cor_), filter(filter_), 
      dedispersion_buffer(dedispersion_buffer_),
      zeropad_buffer(zeropad_buffer_) { 
}

Coherent_dedispersion::~Coherent_dedispersion() {
  // Clear variable in order not to confuse the reference counting
  cur_output = Delay_queue_element();
  while (output_queue->size() > 0)
    output_queue->pop();
}

void
Coherent_dedispersion::do_task() {
  Delay_queue_element input = input_queue->front_and_pop();
  Memory_pool_vector_element<std::complex<FLOAT> > &input_data = input->data;
  const int input_stride = input->stride;
  const int n_dedisp_fft = input_data.size() / input_stride;
  const int n_corr_fft = 
    std::min((size_t) (number_ffts_in_integration - current_fft),
             n_dedisp_fft * fft_dedisp_size() / fft_cor_size());
  // Allocate output buffer
  allocate_element(n_corr_fft);

  for (int i = 0; i < n_dedisp_fft; i++) {
    // Apply dedispersion
    SFXC_MUL_FC(&input_data[i * input_stride], 
                &filter[0], 
                &dedispersion_buffer[0], 
                fft_dedisp_size() / 2 + 1);
    fft.irfft(&dedispersion_buffer[0], 
              &time_buffer[current_buffer][0]);

    // Perform overlap-add
    overlap_add();
    current_buffer = 1 - current_buffer;
  }
  // Write output data
  if (out_pos > 0) {
    cur_output->data.resize(out_pos);
    output_queue->push(cur_output);
  }
}


void
Coherent_dedispersion::empty_output_queue() {
  while (output_queue->size() > 0)
    output_queue->pop();
}

void
Coherent_dedispersion::overlap_add() {
  Memory_pool_vector_element<std::complex<FLOAT> > &out = cur_output->data;
  //NB : fft_size_dedispersion >= fft_size_correlation
  const int step_size = fft_cor_size() / 4;
  const int total_steps_in_input = fft_dedisp_size() / step_size / 2;
  const int nstep_out = fft_cor_size() / 2 / step_size;

  int k = 0;
  for (int n = 0; n < total_steps_in_input; n++) {
    int i, j;
    if (n < total_steps_in_input / 2) {
      i = 3 * fft_dedisp_size() / 4 + n * step_size;
      j = fft_dedisp_size() / 4 + n * step_size;
    } else { 
      i = (n - total_steps_in_input / 2) * step_size;
      j = fft_dedisp_size() / 2 + (n - total_steps_in_input / 2) * step_size;
    }

    // Sum overlapping windows
    SFXC_ADD_F(&time_buffer[current_buffer][i],
               &time_buffer[1 - current_buffer][j],
               &zeropad_buffer[k * step_size],
               step_size);
    k++;
    if (k == nstep_out) {
      if (fft_to_skip > 0) {
        fft_to_skip--;
      } else {
        fft_cor.rfft(&zeropad_buffer[0], &out[out_pos]);
        out_pos += output_stride;
        current_fft += 1;
      }
      k = 0;
      if (current_fft == number_ffts_in_integration)
        break;
    }
  }
}

void
Coherent_dedispersion::allocate_element(int nfft) {
  cur_output = output_memory_pool.allocate();
  cur_output->data.resize(output_stride * nfft);
  cur_output->stride = output_stride;
  out_pos = 0;
}

bool
Coherent_dedispersion::has_work() {
  if (input_queue->empty())
    return false;
  if (output_memory_pool.empty())
    return false;
  return true;
}

void Coherent_dedispersion::connect_to(Delay_queue_ptr buffer) {
  input_queue = buffer;
}

Coherent_dedispersion::Delay_queue_ptr
Coherent_dedispersion::get_output_buffer() {
  SFXC_ASSERT(output_queue != Delay_queue_ptr());
  return output_queue;
}

void 
Coherent_dedispersion::set_parameters(const Correlation_parameters &parameters)
{
  stream_idx = 0;
  while ((stream_idx < parameters.station_streams.size()) &&
         (parameters.station_streams[stream_idx].station_stream != stream_nr))
    stream_idx++;
  if (stream_idx == parameters.station_streams.size()) {
    // Data stream is not participating in current time slice
    return;
  }

  correlation_parameters = parameters;
  output_stride =  fft_cor_size() / 2 + 4; // for allignment
  
  current_fft = 0;
  current_buffer = 0;
  fft_to_skip = (parameters.integration_start.diff(parameters.stream_start) * 
                 sample_rate()  + fft_dedisp_size() / 4) * 2 / fft_rot_size();
  number_ffts_in_integration =
    Control_parameters::nr_correlation_ffts_per_integration(
      (int) parameters.integration_time.get_time_usec(),
      parameters.sample_rate,
      parameters.fft_size_correlation);

  start_time = parameters.integration_start;
  // Initialize buffers
  time_buffer[0].resize(fft_dedisp_size());
  time_buffer[1].resize(fft_dedisp_size());
}

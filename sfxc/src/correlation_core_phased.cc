#include "correlation_core_phased.h"
#include "output_header.h"
#include <utils.h>

Correlation_core_phased::Correlation_core_phased()
{
}

Correlation_core_phased::~Correlation_core_phased()
{
}

void
Correlation_core_phased::do_task() {
  SFXC_ASSERT(has_work());

  if (current_fft % 1000 == 0) {
    PROGRESS_MSG("node " << node_nr_ << ", "
                 << current_fft << " of " << number_ffts_in_slice);
  }

  if (current_fft % number_ffts_in_slice == 0) {
    integration_initialise();
  }

  for (size_t i = 0; i < number_input_streams(); i++) {
    int stream = station_stream(i);
    input_elements[i] = &input_buffers[stream]->front()->data[0];
  }
  const int first_stream = station_stream(0);
  const int stride = input_buffers[0]->front()->stride;
  const int nbuffer = input_buffers[0]->front()->data.size() / stride;
  for (size_t buf_idx = 0; buf_idx < nbuffer * stride ; buf_idx += stride){
    // Process the data of the current fft
    integration_step(accumulation_buffers, buf_idx);
    current_fft++;
  }

  for (size_t i = 0; i < number_input_streams(); i++) {
    int stream = station_stream(i);
    input_buffers[stream]->pop();
  }

  if (current_fft == number_ffts_in_slice) {
    PROGRESS_MSG("node " << node_nr_ << ", "
                 << current_fft << " of " << number_ffts_in_slice);

    int source = sources[delay_tables[first_stream].get_source(0)];
    integration_write(accumulation_buffers, 0, source, 1);
    current_integration++;
  }
}

void
Correlation_core_phased::set_parameters(const Correlation_parameters &parameters,
                                        std::vector<Delay_table_akima> &delays,
                                        std::vector<std::vector<double> > &uvw,
					int node_nr)
{
  node_nr_ = node_nr;
  current_integration = 0;
  current_fft = 0;
  delay_tables = delays;
  uvw_table = uvw;

  correlation_parameters = parameters;

  create_baselines(parameters);
  if (input_elements.size() != number_input_streams()) {
    input_elements.resize(number_input_streams());
  }
  if (input_conj_buffers.size() != number_input_streams()) {
    input_conj_buffers.resize(number_input_streams());
    for(int i = 0; i < number_input_streams(); i++)
      input_conj_buffers[i].resize(fft_size() + 1);
  }
  n_flagged.resize(baselines.size());
}

void
Correlation_core_phased::create_baselines(const Correlation_parameters &parameters){
  number_ffts_in_slice = parameters.slice_size /
    parameters.fft_size_correlation;
  // One less because of the overlapping windows
  if ((parameters.window != SFXC_WINDOW_NONE) &&
      (parameters.window != SFXC_WINDOW_PFB))
    number_ffts_in_slice -= 1;

  number_ffts_in_sub_integration =
    Control_parameters::nr_correlation_ffts_per_integration(
      (int) parameters.sub_integration_time.get_time_usec(),
      parameters.sample_rate,
      parameters.fft_size_correlation);

  baselines.clear();
  baselines.resize(number_ffts_in_slice / number_ffts_in_sub_integration + 1);
}

void Correlation_core_phased::integration_initialise() {
  int num_sub_integrations =
    (number_ffts_in_slice / number_ffts_in_sub_integration + 1);
  if (accumulation_buffers.size() != num_sub_integrations) {
    accumulation_buffers.resize(num_sub_integrations);
    for (size_t i = 0; i < accumulation_buffers.size(); i++) {
      accumulation_buffers[i].resize(fft_size() + 1);
    }
  }

  for (size_t i = 0; i < accumulation_buffers.size(); i++) {
    size_t size = accumulation_buffers[i].size() * sizeof(std::complex<FLOAT>);
    memset(&accumulation_buffers[i][0], 0, size);
  }

  memset(&n_flagged[0], 0, sizeof(std::pair<int64_t,int64_t>)*n_flagged.size());
}

void Correlation_core_phased::integration_step(std::vector<Complex_buffer> &integration_buffer, int buf_idx) {
#ifndef DUMMY_CORRELATION
  int sub_integration = current_fft / number_ffts_in_sub_integration;
  SFXC_ASSERT(sub_integration < integration_buffer.size());
  for (size_t i = 0; i < number_input_streams(); i++) {
    SFXC_ADD_FC_I(&input_elements[i][buf_idx],
		  &integration_buffer[sub_integration][0], fft_size() + 1);
  }
#endif // DUMMY_CORRELATION
}

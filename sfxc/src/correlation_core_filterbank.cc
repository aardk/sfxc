#include "correlation_core_filterbank.h"
#include "output_header.h"
#include "bandpass.h"
#include "utils.h"
#include <set>

Correlation_core_filterbank::Correlation_core_filterbank(): 
     DM(0), no_intra_channel_dedispersion(false)
{
}

Correlation_core_filterbank::~Correlation_core_filterbank() {
}

void Correlation_core_filterbank::do_task() {
  SFXC_ASSERT(has_work());
  if (current_fft % 1000 == 0) {
    PROGRESS_MSG("node " << node_nr_ << ", "
                 << current_fft << " of " << number_ffts_in_integration);
  }

  if (current_fft == 0) {
    integration_initialise();
  }
  SFXC_ASSERT(current_fft < number_ffts_in_integration);

  for (size_t i=0; i < number_input_streams_in_use(); i++) {
    int j = streams_in_scan[i];
    input_elements[i] = &input_buffers[j]->front()->data[0];
    if (input_buffers[j]->front()->data.size() > input_conj_buffers[i].size())
      input_conj_buffers[i].resize(input_buffers[j]->front()->data.size());
  }
  const int first_stream = streams_in_scan[0];
  const int stride = input_buffers[first_stream]->front()->stride;
  const int nbuffer = std::min((size_t)number_ffts_in_integration - current_fft,
                               input_buffers[first_stream]->front()->data.size() / stride);
  // Process the data from the current fft buffer
  for (int buf_idx=0; (buf_idx < nbuffer) && (current_fft < number_ffts_in_integration); buf_idx++) {
    integration_step(accumulation_buffers, buf_idx * stride);
    current_fft++;
    if (current_fft == number_ffts_in_integration) {
      PROGRESS_MSG("node " << node_nr_ << ", "
                   << current_fft << " of " << number_ffts_in_integration);

      sub_integration();
      for(int i = 0 ; i < phase_centers.size(); i++){
        int source_nr;
        int pcenter = 0;
        source_nr = streams_in_scan[i];
        integration_write_headers(pcenter, source_nr);
        integration_write_subints(phase_centers[i]);
      }
      current_integration++;
    }else if(current_fft >= next_sub_integration*number_ffts_in_sub_integration){
      sub_integration();
      next_sub_integration++;
    } else if(no_intra_channel_dedispersion) {
      sub_integration();
    }
  }
  for (size_t i=0, nstreams=number_input_streams_in_use(); i<nstreams; i++){
    int stream = streams_in_scan[i];
    input_buffers[stream]->pop();
  }
}

void
Correlation_core_filterbank::set_parameters(
                                  const Correlation_parameters &parameters,
                                  std::vector<Delay_table_akima> &delays,
                                  std::vector<std::vector<double> > &uvw,
                                  int node_nr, double DM_, 
                                  bool no_intra_channel_dedispersion_)
{
  node_nr_ = node_nr;
  current_integration = 0;
  current_fft = 0;
  delay_tables = delays;
  uvw_table = uvw;
  DM = DM_;
  no_intra_channel_dedispersion = no_intra_channel_dedispersion_;
  correlation_parameters = parameters;

  create_baselines(parameters);
  // NB: we overide the value from create_baselines
  number_ffts_in_integration = parameters.slice_size / parameters.fft_size_correlation;
  if (input_elements.size() != number_input_streams_in_use()) {
    input_elements.resize(number_input_streams_in_use());
  }
  if (input_conj_buffers.size() != number_input_streams_in_use()) {
    input_conj_buffers.resize(number_input_streams_in_use());
  }
  // Read calibration tables
  if (old_fft_size != fft_size()){
    old_fft_size = fft_size();
    if (cltable_name != std::string())
      cltable.open_table(cltable_name, fft_size());
    if (bptable_name != std::string())
      bptable.open_table(bptable_name, fft_size(), true);
  }
  n_flagged.resize(baselines.size());
  get_input_streams();
  create_channel_offsets();
}

void Correlation_core_filterbank::integration_initialise() {
  number_output_products = ceil(correlation_parameters.integration_time / correlation_parameters.sub_integration_time);
  previous_fft = 0;

  int nstreams = streams_in_scan.size();
  if(phase_centers.size() != nstreams)
    phase_centers.resize(nstreams);

  for(int i = 0 ; i < phase_centers.size() ; i++){
    if (phase_centers[i].size() != number_output_products){
      phase_centers[i].resize(number_output_products);
      for(int j = 0 ; j < phase_centers[i].size() ; j++){
        phase_centers[i][j].resize(fft_size() + 1);
      }
    }
  }

  for(int i = 0 ; i < phase_centers.size() ; i++){
    for (int j = 0; j < phase_centers[i].size(); j++) {
      SFXC_ASSERT(phase_centers[i][j].size() == fft_size() + 1);
      size_t size = phase_centers[i][j].size() * sizeof(std::complex<FLOAT>);
      memset(&phase_centers[i][j][0], 0, size);
    }
  }

  if (accumulation_buffers.size() != baselines.size()) {
    accumulation_buffers.resize(baselines.size());
    for (size_t i=0; i<accumulation_buffers.size(); i++) {
      accumulation_buffers[i].resize(fft_size() + 1);
    }
  }

  SFXC_ASSERT(accumulation_buffers.size() == baselines.size());
  for (size_t i=0; i<accumulation_buffers.size(); i++) {
    SFXC_ASSERT(accumulation_buffers[i].size() == fft_size() + 1);
    size_t size = accumulation_buffers[i].size() * sizeof(std::complex<FLOAT>);
    memset(&accumulation_buffers[i][0], 0, size);
  }
  memset(&n_flagged[0], 0, sizeof(std::pair<int64_t,int64_t>)*n_flagged.size());
  next_sub_integration = 1;
}

void Correlation_core_filterbank::
integration_step(std::vector<Complex_buffer> &integration_buffer, 
                 int buf_idx) 
{
  for (size_t i = 0; i < number_input_streams_in_use(); i++) {
    // get the complex conjugates of the input
    SFXC_CONJ_FC(&input_elements[i][buf_idx], &input_conj_buffers[i][buf_idx], 
                 fft_size() + 1);
  }
  
  for (size_t i = 0; i < number_input_streams_in_use(); i++) {
    // Cross correlations
    std::pair<size_t,size_t> &stations = baselines[i];
    //SFXC_ASSERT(stations.first == stations.second);
    SFXC_ADD_PRODUCT_FC(/* in1 */ &input_elements[stations.first][buf_idx], 
                    /* in2 */ &input_conj_buffers[stations.second][buf_idx],
                    /* out */ &integration_buffer[i][0], fft_size() + 1);
  }
}

void Correlation_core_filterbank::integration_write_subints(std::vector<Complex_buffer> &integration_buffer) {
  SFXC_ASSERT(accumulation_buffers.size() == baselines.size());

  int polarisation = (correlation_parameters.polarisation == 'R')? 0 : 1;

  std::vector<float> float_buffer;
  float_buffer.resize(number_channels() + 1);
  Output_header_baseline hbaseline;

  size_t n = fft_size() / number_channels();
  for (size_t i = 0; i < integration_buffer.size(); i++) {
    for (size_t j = 0; j < number_channels(); j++) {
      float_buffer[j] = integration_buffer[i][j * n].real();
      for (size_t k = 1; k < n ; k++)
        float_buffer[j] += integration_buffer[i][j * n + k].real();
      float_buffer[j] /= n;
    }
    if (RANK_OF_NODE == -10)
      std::cerr << i << " : " << float_buffer[0] << ", " << float_buffer[1] << ", " << float_buffer[2] 
                              << integration_buffer[i][0] << ", " << integration_buffer[i][1] << ", " 
                              << integration_buffer[i][2]  <<"\n";
    float_buffer[number_channels()] = 
                            integration_buffer[i][number_channels()*n].real();
    hbaseline.weight = 1;                        // The number of good samples
    hbaseline.station_nr1 = 0;
    hbaseline.station_nr2 = 0;

    // Polarisation for the first station
    SFXC_ASSERT((polarisation == 0) || (polarisation == 1)); // (RCP: 0, LCP: 1)
    hbaseline.polarisation1 = polarisation;
    hbaseline.polarisation2 = polarisation;
    // Upper or lower sideband (LSB: 0, USB: 1)
    if (correlation_parameters.sideband=='U') {
      hbaseline.sideband = 1;
    } else {
      SFXC_ASSERT(correlation_parameters.sideband == 'L');
      hbaseline.sideband = 0;
    }
    // The number of the channel in the vex-file,
    hbaseline.frequency_nr = (unsigned char)correlation_parameters.frequency_nr;
    // sorted increasingly
    // 1 byte left:
    hbaseline.empty = ' ';

    int nWrite = sizeof(hbaseline);
    writer->put_bytes(nWrite, (char *)&hbaseline);
    writer->put_bytes((number_channels() + 1) * sizeof(float),
                      ((char*)&float_buffer[0]));
  }
}

void 
Correlation_core_filterbank::sub_integration(){
  Time tfft(0., correlation_parameters.sample_rate); 
  tfft.inc_samples(fft_size());
  //const Time tmid = correlation_parameters.stream_start + 
  //                  tfft*(current_fft + 0.5); 
  //calibrate(accumulation_buffers, tmid); 

  Time tstart = correlation_parameters.stream_start - 
                 correlation_parameters.integration_start;
  double tmid = (tstart +
                 tfft*(previous_fft+(current_fft-previous_fft)/2.)).get_time_usec();
  double sub_integration_time = correlation_parameters.sub_integration_time.get_time_usec();
  const int n_fft = fft_size() + 1;
  const int n_phase_centers = phase_centers.size();
  const int n_station = number_input_streams_in_use();
  const int n_baseline = accumulation_buffers.size();
  if((RANK_OF_NODE == -5) && (next_sub_integration >= 15624))
                         std::cout << "nbaseline = " << n_baseline
                                   << ", subint = " << next_sub_integration-1
                                   << " / " << phase_centers[0].size()
                                   << "fft = " << current_fft
                                   << ", nfft_per_sub="<< number_ffts_in_sub_integration
                                   << ", tmid = " << (int64_t)tmid
                                   << "\n";
  for(int i = 0 ; i < n_station ; i++){
    std::pair<size_t,size_t> &inputs = baselines[i];
    int station1 = streams_in_scan[inputs.first];
    int station2 = streams_in_scan[inputs.second];

    for (int j = 0; j <= fft_size(); j++) {
      int bin = floor((tmid + offsets[j]) / sub_integration_time);
      if ((bin >= 0) and (bin < number_output_products))
          phase_centers[i][bin][j] += accumulation_buffers[i][j];
    }
  }
  // Clear the accumulation buffers
  for (size_t i=0; i<accumulation_buffers.size(); i++) {
    SFXC_ASSERT(accumulation_buffers[i].size() == n_fft);
    size_t size = accumulation_buffers[i].size() * sizeof(std::complex<FLOAT>);
    memset(&accumulation_buffers[i][0], 0, size);
  }
  previous_fft = current_fft;
}

void
Correlation_core_filterbank::create_channel_offsets() {
  offsets.resize(fft_size());
  if (!no_intra_channel_dedispersion) {
    for (int i = 0; i < offsets.size(); i++)
      offsets[i] = 0;
    return;
  }
  int m = fft_size() / number_channels();
  // Find the time offsets between frequency components
  int sb = correlation_parameters.sideband == 'L' ? -1 : 1;
  double base_freq = correlation_parameters.channel_freq*1e-6; // [MHz]
  double fmid = (correlation_parameters.channel_freq + sb * correlation_parameters.bandwidth / 2) * 1e-6; // [MHz]
  for (int i = 0; i <= number_channels(); i++) {
    double f_i = base_freq + sb * i * 1e-6 * correlation_parameters.bandwidth / number_channels();
    for (int j = 0; j < m; j++) {
      int k = i * m + j;
      if (k <= fft_size()) {
        offsets[k] = 4.148808e9 * DM * ( 1. / (f_i * f_i) - 1. / (fmid*fmid));
        std::cerr.precision(16);
        if (RANK_OF_NODE == -5)
          std::cerr << "offsets[" << k << "] = " << offsets[k] << " ;f_" <<i << " = " << f_i << ", fmid = " << fmid<< "\n";
      }
    }
  }
}

#include "correlation_core_voltages.h"
#include "output_header.h"
#include "bandpass.h"
#include "utils.h"
#include <set>

Correlation_core_voltages::Correlation_core_voltages() : ffts_to_skip(0)
{
}

Correlation_core_voltages::~Correlation_core_voltages() {
}

void Correlation_core_voltages::connect_to(size_t stream, Channel_queue_ptr buffer) {
  if (stream >= input_buffers.size()) {
    input_buffers.resize(stream+1);
  }
  input_buffers[stream] = buffer;
}

bool Correlation_core_voltages::has_work() {
  for (size_t i=0, nstreams=number_input_streams_in_use(); i < nstreams; i++) {
    int stream = streams_in_scan[i];
    if (input_buffers[stream]->empty())
      return false;
  }
  return true;
}

void Correlation_core_voltages::do_task() {
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
  }
  const int first_stream = streams_in_scan[0];
  const int nbuffer_received = input_buffers[first_stream]->front()->data.size() / fft_size();
  if (ffts_to_skip >= nbuffer_received) {
    ffts_to_skip -= nbuffer_received;
    for (size_t i=0, nstreams=number_input_streams_in_use(); i<nstreams; i++){
      int stream = streams_in_scan[i];
      input_buffers[stream]->pop();
    }
    return;
  }
  const int first = ffts_to_skip;
  ffts_to_skip = 0;
  const int nbuffer = std::min(number_ffts_in_integration - current_fft,
                               nbuffer_received - first);
  // Process the data from the current fft buffer
  integration_step(current_fft * fft_size(), first * fft_size(), nbuffer * fft_size());
  current_fft += nbuffer;
  if (current_fft == number_ffts_in_integration) {
    PROGRESS_MSG("node " << node_nr_ << ", "
                 << current_fft << " of " << number_ffts_in_integration);

    for(int i = 0 ; i < station_voltage.size(); i++){
      int pcenter = 0;
      int source_nr = streams_in_scan[i];
      integration_write_headers(pcenter, source_nr);
      integration_write_subints(station_voltage[i]);
    }
    current_integration++;
  }
  for (size_t i=0, nstreams=number_input_streams_in_use(); i<nstreams; i++){
    int stream = streams_in_scan[i];
    input_buffers[stream]->pop();
  }
}

void
Correlation_core_voltages::set_parameters(
                                  const Correlation_parameters &parameters,
                                  std::vector<std::vector<double> > &uvw,                                  
                                  int node_nr)
{
  node_nr_ = node_nr;
  current_integration = 0;
  current_fft = 0;
  correlation_parameters = parameters;
  uvw_table = uvw;

  ffts_to_skip = (int)round((parameters.integration_start - parameters.stream_start).get_time_usec() * 
                            (parameters.sample_rate * 1e-6)) / parameters.fft_size_correlation;
  //ffts_to_skip = 0;

  create_baselines(parameters);
  if (input_elements.size() != number_input_streams_in_use()) {
    input_elements.resize(number_input_streams_in_use());
  }
  // Read calibration tables
  // FIXME: Calibration doesn't work in this mode
  if (old_fft_size != fft_size()){
    old_fft_size = fft_size();
    if (cltable_name != std::string())
      cltable.open_table(cltable_name, fft_size());
    if (bptable_name != std::string())
      bptable.open_table(bptable_name, fft_size(), true);
  }
  n_flagged.resize(baselines.size());
  get_input_streams();
}

void
Correlation_core_voltages::create_baselines(const Correlation_parameters &parameters) {
  number_ffts_in_integration =
    Control_parameters::nr_ffts_to_output_node(
      parameters.integration_time,
      parameters.sample_rate,
      parameters.fft_size_correlation);
  if(parameters.window != SFXC_WINDOW_NONE)
    number_ffts_in_integration -= 1;

  baselines.clear();
  // Autos
  for (size_t sn = 0 ; sn < number_input_streams_in_use(); sn++) {
    baselines.push_back(std::pair<int,int>(sn,sn));
  }
}

void Correlation_core_voltages::integration_initialise() {
  // FIXME this assumes no mixed single and double polarizitions
  if (correlation_parameters.cross_polarize)
    number_output_products = 2; 
  else
    number_output_products = 1;
  previous_fft = 0;

  int nstreams = streams_in_scan.size();
  int nstations = (correlation_parameters.cross_polarize) ? nstreams / 2 : nstreams;
  if(station_voltage.size() != nstations)
    station_voltage.resize(nstations);

  int nprod = number_ffts_in_integration * fft_size() * (baselines.size() / nstations);
  int npol = correlation_parameters.cross_polarize ? 2 : 1; 
  for(int i = 0 ; i < station_voltage.size() ; i++) {
    station_voltage[i].resize(npol);
    for (int j = 0; j < station_voltage[i].size(); j++) {
      if (station_voltage[i][j].size() != nprod) { 
        station_voltage[i][j].resize(nprod);
      }
    }
  }

  for(int i = 0 ; i < station_voltage.size() ; i++){
    for (int j = 0; j < station_voltage[i].size(); j++) {
      size_t size = station_voltage[i][j].size() * sizeof(float);
      memset(&station_voltage[i][j][0], 0, size);
    }
  }

  memset(&n_flagged[0], 0, sizeof(std::pair<int64_t,int64_t>)*n_flagged.size());
}

void Correlation_core_voltages::
integration_step(int ostart, int istart, int nsamples) 
{
  const int nstation = station_voltage.size();
  for (size_t i = 0; i < baselines.size(); i++) {
    // Cross correlations
    std::pair<size_t,size_t> &inputs = baselines[i];
    int station = inputs.first % nstation;
    int ipol = (inputs.first < nstation) ? 0 : 1; 

    memcpy(&station_voltage[station][ipol][ostart],
           &input_elements[inputs.first][istart], 
           nsamples * sizeof(FLOAT));
  }
}

void Correlation_core_voltages::integration_write_subints(std::vector< std::vector<float> > &voltage) {
  Output_header_baseline hbaseline;

  int pol_ref = (correlation_parameters.polarisation == 'R')? 0 : 1;
  int npol = (correlation_parameters.cross_polarize) ? 3 : 1;
  int pol1[4] = {pol_ref, 1-pol_ref, 0};
  int pol2[4] = {pol_ref, 1-pol_ref, 1};
  for (int ipol = 0; ipol < npol; ipol++) {
    hbaseline.weight = 1;                        // The number of good samples
    hbaseline.station_nr1 = 0;
    hbaseline.station_nr2 = 0;

    // Polarisation for the first station
    hbaseline.polarisation1 = pol1[ipol];
    hbaseline.polarisation2 = pol2[ipol];
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
    writer->put_bytes(voltage[ipol].size() * sizeof(float),
                      ((char*)&voltage[ipol][0]));
  }
}
